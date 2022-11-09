from os import makedirs
from os.path import basename
from os.path import join
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import as_completed

from get_biomes.utils import get_arguments, fast_flatten
import logging
import pandas as pd
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import requests
import tqdm
import os
from functools import partial
from multiprocessing import Pool
from pySmartDL import SmartDL
import math


def get_data(url, path):
    url = f"http://{url}"
    filename = basename(url)
    # construct the output path
    outpath = join(path, filename)
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
    r = http.get(url, timeout=10, stream=True)
    # Parse the data
    if r.status_code == 200:
        with open(outpath, "wb") as fastq:
            for chunk in r.iter_content(chunk_size=20480):
                if chunk:
                    fastq.write(chunk)
        return (url, outpath)
    else:
        return (url, None)


log = logging.getLogger("my_logger")


def download(args):

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(levelname)s ::: %(asctime)s ::: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    args = get_arguments()

    # Check if outfile exists and remove it
    # get full path to outfile
    logging.getLogger("my_logger").setLevel(
        logging.DEBUG if args.debug else logging.INFO
    )
    logging.getLogger("urllib3").setLevel(logging.ERROR)
    logging.getLogger("jsonapi_client").setLevel(logging.ERROR)

    biomes = pd.read_csv(args.input, sep="\t", header=0)
    urls = fast_flatten(biomes["fastq_ftp"].str.split(";", expand=True).values.tolist())

    # if output directory does not exist, create it if not delete
    if not os.path.exists(args.outdir):
        makedirs(args.outdir)
    download_report = os.path.join(args.outdir, "download_report.tsv")
    # Check if download_report file exists
    if os.path.exists(download_report):
        report = pd.read_csv(download_report, sep="\t", header=0)

    if args.clean_output and os.path.exists(args.outdir):
        log.info("Output directory already exists, deleting it...")
        for file in os.listdir(args.outdir):
            os.remove(os.path.join(args.outdir, file))

    files = []
    log.info("Downloading files...")
    for url in tqdm.tqdm(
        urls,
        total=len(urls),
        desc="Files downloaded",
        ncols=80,
        leave=False,
    ):
        obj = SmartDL(
            f"http://{url}", args.outdir, threads=args.threads, progress_bar=False
        )
        obj.start(blocking=False)
        with tqdm.tqdm(
            total=obj.get_final_filesize(human=False),
            leave=False,
            desc=basename(url),
            ncols=80,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
        ) as pbar:
            prev = 0
            while not obj.isFinished():
                if obj.get_dl_size(human=False) - prev > 1:
                    pbar.update(obj.get_dl_size(human=False) - prev)
                    prev = obj.get_dl_size(human=False)
        files.append((basename(url), obj.isSuccessful(), obj.get_errors(), url))
    print(files)
    log.info("Done!")
