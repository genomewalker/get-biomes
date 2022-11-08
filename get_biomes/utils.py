import argparse
import sys
import gzip
import os
import logging
import pandas as pd
from multiprocessing import Pool
from functools import partial
from os import devnull
import tqdm
from get_biomes import __version__
import time
from itertools import chain
import io
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import requests
from contextlib import contextmanager, redirect_stderr, redirect_stdout
import json

log = logging.getLogger("my_logger")
log.setLevel(logging.INFO)
timestr = time.strftime("%Y%m%d-%H%M%S")


API_BASE = "https://www.ebi.ac.uk/metagenomics/api/v1"
API_BASE_ENA = "https://www.ebi.ac.uk/ena/portal/api/"


def is_debug():
    return logging.getLogger("my_logger").getEffectiveLevel() == logging.DEBUG


def check_values(val, minval, maxval, parser, var):
    value = float(val)
    if value < minval or value > maxval:
        parser.error(
            "argument %s: Invalid value %s. Range has to be between %s and %s!"
            % (
                var,
                value,
                minval,
                maxval,
            )
        )
    return value


# From: https://note.nkmk.me/en/python-check-int-float/
def is_integer(n):
    try:
        float(n)
    except ValueError:
        return False
    else:
        return float(n).is_integer()


# function to check if the input value has K, M or G suffix in it
def check_suffix(val, parser, var):
    units = ["K", "M", "G"]
    unit = val[-1]
    value = int(val[:-1])

    if is_integer(value) & (unit in units) & (value > 0):
        return val
    else:
        parser.error(
            "argument %s: Invalid value %s. Memory has to be an integer larger than 0 with the following suffix K, M or G"
            % (var, val)
        )


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {
        "gz": (b"\x1f", b"\x8b", b"\x08"),
        "bz2": (b"\x42", b"\x5a", b"\x68"),
        "zip": (b"\x50", b"\x4b", b"\x03", b"\x04"),
    }
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, "rb")
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = "plain"
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == "bz2":
        sys.exit("Error: cannot use bzip2 format - use gzip instead")
        sys.exit("Error: cannot use zip format - use gzip instead")
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == "gz":
        return gzip.open
    else:  # plain text
        return open


# From: https://stackoverflow.com/a/11541450
def is_valid_file(parser, arg, var):
    if not os.path.exists(arg):
        parser.error("argument %s: The file %s does not exist!" % (var, arg))
    else:
        return arg


filters = [
    "experiment_type",
    "biome_name",
    "lineage",
    "geo_loc_name",
    "latitude_gte",
    "latitude_lte",
    "longitude_gte",
    "longitude_lte",
    "species",
    "instrument_model",
    "instrument_platform",
    "metadata_key",
    "metadata_value_gte",
    "metadata_value_lte",
    "metadata_value",
    "environment_material",
    "environment_feature",
    "study_accession",
    "include",
]
# From https://stackoverflow.com/a/59617044/15704171
def convert_list_to_str(lst):
    n = len(lst)
    if not n:
        return ""
    if n == 1:
        return lst[0]
    return ", ".join(lst[:-1]) + f" or {lst[-1]}"


def is_valid_filter(parser, arg, var):
    arg = json.loads(arg)
    # check if the dictionary keys are in the mdmg header list
    for key in arg.keys():
        if key not in filters:
            parser.error(
                f"argument {var}: Invalid value {key}.\n"
                f"Valid values are: {convert_list_to_str(filters)}"
            )
    return arg


defaults = {
    "threads": 1,
    "outfile": "output.tsv",
}

help_msg = {
    "biomes": "A txt file containing MGnify biomes. Ex: root:Environmental:Aquatic:Marine",
    "threads": "Number of threads to use",
    "mg_filter": f"Key-value pairs to filter the MGnify metadata. Valid values are: {convert_list_to_str(filters)}",
    "outfile": "Output file name",
    "help": "Help message",
    "debug": "Print debug messages",
    "version": "Print program version",
}


def get_arguments(argv=None):
    parser = argparse.ArgumentParser(
        description="A simple tool to get biomes related data from MGnify and ENA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
    )

    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")

    required.add_argument(
        "-b",
        "--biomes",
        dest="biomes",
        type=lambda x: is_valid_file(parser, x, "--biomes"),
        help=help_msg["biomes"],
        required=True,
    )
    optional.add_argument(
        "-f",
        "--filter",
        type=lambda x: is_valid_filter(parser, x, "--filter"),
        dest="mg_filter",
        default={
            "experiment_type": "metagenomic",
            "instrument_platform": "illumina",
        },
        help=help_msg["mg_filter"],
        required=False,
    )
    optional.add_argument(
        "--output",
        dest="outfile",
        help=help_msg["outfile"],
        required=False,
        default=defaults["outfile"],
    )
    optional.add_argument(
        "-t",
        "--threads",
        type=lambda x: int(
            check_values(x, minval=1, maxval=1000, parser=parser, var="--threads")
        ),
        dest="threads",
        default=1,
        help=help_msg["threads"],
        required=False,
    )
    optional.add_argument(
        "--debug",
        dest="debug",
        action="store_true",
        help=help_msg["debug"],
        required=False,
    )
    optional.add_argument(
        "--version",
        action="version",
        version="%(prog)s " + __version__,
        help=help_msg["version"],
    )
    optional.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="show this help message and exit",
    )
    # reference_lengths

    args = parser.parse_args(None if sys.argv[1:] else ["-h"])
    return args


@contextmanager
def suppress_stdout():
    """A context manager that redirects stdout and stderr to devnull"""
    with open(devnull, "w") as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)


def applyParallel(dfGrouped, func, threads, parms):
    p = Pool(threads)
    func = partial(func, parms=parms)
    ret_list = tqdm.tqdm(
        p.map(func, [group for name, group in dfGrouped]),
        total=len([group for name, group in dfGrouped]),
    )
    p.close()
    p.join()
    return pd.concat(ret_list)


def fast_flatten(input_list):
    return list(chain.from_iterable(input_list))


def concat_df(frames):
    COLUMN_NAMES = frames[0].columns
    df_dict = dict.fromkeys(COLUMN_NAMES, [])
    for col in COLUMN_NAMES:
        extracted = (frame[col] for frame in frames)
        # Flatten and save to df_dict
        df_dict[col] = fast_flatten(extracted)
    df = pd.DataFrame.from_dict(df_dict)[COLUMN_NAMES]
    return df


def initializer(init_data):
    global parms
    parms = init_data


def do_parallel(parms, lst, func, threads):
    if is_debug():
        dfs = list(map(partial(func, parms=parms), lst))
    else:
        p = Pool(threads, initializer=initializer, initargs=(parms,))
        c_size = calc_chunksize(threads, len(lst))
        dfs = list(
            tqdm.tqdm(
                p.imap_unordered(
                    partial(func, parms=parms),
                    lst,
                    chunksize=c_size,
                ),
                total=len(lst),
                leave=False,
                ncols=80,
                desc=f"Components processed",
            )
        )
        p.close()
        p.join()
    return concat_df(dfs)


# from https://stackoverflow.com/questions/53751050/python-multiprocessing-understanding-logic-behind-chunksize/54032744#54032744
def calc_chunksize(n_workers, len_iterable, factor=4):
    """Calculate chunksize argument for Pool-methods.

    Resembles source-code within `multiprocessing.pool.Pool._map_async`.
    """
    chunksize, extra = divmod(len_iterable, n_workers * factor)
    if extra:
        chunksize += 1
    return chunksize


def get_mgnify_data(sample):
    instrument_model = None
    sequencing_method = None
    investigation_type = None
    for k in sample.sample_metadata:
        if k["key"] == "instrument model":
            instrument_model = k["value"]
        elif k["key"] == "sequencing method":
            sequencing_method = k["value"]
        elif k["key"] == "investigation type":
            investigation_type = k["value"]
    d = {
        "accession": sample.accession,
        "sample_accession": sample.biosample,
        "sample_name": sample.sample_name,
        "longitude": sample.longitude,
        "latitude": sample.latitude,
        "geo_loc_name": sample.geo_loc_name,
        "studies": ",".join([study.accession for study in sample.studies]),
        "biome": sample.biome.id,
        "sample_desc": sample.sample_desc,
        "environment_biome": sample.environment_biome,
        "environment_feature": sample.environment_feature,
        "environment_material": sample.environment_material,
        "instrument_model": instrument_model,
        "sequencing_method": sequencing_method,
        "investigation_type": investigation_type,
    }

    df_mg = pd.DataFrame({k: [v] for k, v in d.items()})
    return df_mg


def get_data(url):
    # Get the data from the API
    retry_strategy = Retry(
        total=5,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504, 403],
        allowed_methods=["HEAD", "GET", "OPTIONS"],
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    http = requests.Session()
    http.mount("https://", adapter)
    http.mount("http://", adapter)
    r = http.get(url, timeout=10)
    # Parse the data
    if r.status_code != 200:
        r = io.StringIO("")
    else:
        r = io.StringIO(r.content.decode("utf-8"))
    return r


def get_ena_data(accession):
    ena_url = f"{API_BASE_ENA}filereport?accession={accession}&result=read_run&fields=sample_accession,study_accession,experiment_accession,run_accession,fastq_ftp&format=tsv&download=true&limit=0"
    data = get_data(ena_url)
    data.seek(0)
    try:
        df = pd.read_csv(data, sep="\t")
    except pd.errors.EmptyDataError:
        df = None
    return df
    # return pd.read_csv(ena_url, sep="\t")
