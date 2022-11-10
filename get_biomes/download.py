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
import urllib

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
        try:
            obj.start(blocking=False)
        except Exception:
            files.append((basename(url), obj.isSuccessful(), obj.get_errors(), url))
            continue

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
    # Update report success with the new report
    if os.path.exists(download_report):
        report = pd.DataFrame(files, columns=["filename", "success", "error", "url"])
        if report_success.shape[0] > 0:
            report = pd.concat([report, report_success], axis=0)
    else:
        report = pd.DataFrame(files, columns=["filename", "success", "error", "url"])

    report.to_csv(download_report, sep="\t", index=False)
    log.info("Done!")
