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
import wget
from multiprocessing.pool import ThreadPool


def download_url(url, path):
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

    file_name = os.path.join(path, basename(url))

    if r.status_code == requests.codes.ok:
        with open(file_name, "wb") as f:
            for data in r:
                f.write(data)
        return basename(url), True, r.status_code, url
    else:
        return basename(url), False, r.status_code, url


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

    # Check if download_report file exists
    download_report = os.path.join(args.outdir, "download_report.tsv")
    if os.path.exists(download_report):
        report = pd.read_csv(download_report, sep="\t", header=0)
        report_success = report[report["success"] == True]
        report_failed = report[report["success"] == False]
        log.info(
            f"Found {report_success.shape[0]} successful downloads and {report_failed.shape[0]} failed downloads"
        )
        urls = list(set(urls) - set(report_success["url"]))
        if len(urls) == 0:
            logging.info("All files have been downloaded successfully")
            return

    if args.clean_output and os.path.exists(args.outdir):
        log.info("Output directory already exists, deleting it...")
        for file in os.listdir(args.outdir):
            os.remove(os.path.join(args.outdir, file))

    files = []
    log.info("Downloading files...")
    # Run 5 multiple threads. Each call will take the next element in urls list
    results = ThreadPool(args.threads).imap_unordered(
        partial(download_url, path=args.outdir), urls
    )

    for r in tqdm.tqdm(
        results, total=len(urls), desc="Files downloaded", ncols=80, leave=False
    ):
        files.append(r)

    if os.path.exists(download_report):
        report = pd.DataFrame(files, columns=["filename", "success", "error", "url"])
        if report_success.shape[0] > 0:
            report = pd.concat([report, report_success], axis=0)
    else:
        report = pd.DataFrame(files, columns=["filename", "success", "error", "url"])

    report.to_csv(download_report, sep="\t", index=False)
    log.info("Done!")
