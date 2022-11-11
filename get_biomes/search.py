"""
This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not,
see <https://www.gnu.org/licenses/>.
"""


import logging
import pandas as pd
from get_biomes.utils import (
    get_arguments,
    concat_df,
    get_mgnify_data,
    get_ena_data,
    create_output_files,
    search_str,
)
from get_biomes.defaults import API_BASE
from jsonapi_client import Filter, Session
from urllib.parse import urlencode
import urllib.parse
import tqdm
from multiprocessing import Pool
import os
from pathlib import Path
from functools import partial

log = logging.getLogger("my_logger")


def search(args):

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

    biomes = pd.read_csv(args.biomes, sep="\t", header=None)
    biomes = biomes[0].tolist()
    output_files = create_output_files(args.prefix, args.biomes, biomes)

    if args.clean_output:
        logging.info("Removing old output files...")
        for file in output_files:
            if os.path.exists(output_files[file]):
                os.remove(output_files[file])
    # API call
    with Session(API_BASE) as session:

        # configure the filters
        if args.mg_filter is not None:
            params = args.mg_filter
            api_filter = Filter(urlencode(params))

        for biome in biomes:
            log.info(f"Processing biome: {biome}...")
            if Path(output_files[biome]).is_file():
                log.info(f"File {output_files[biome]} already exists, skipping")
                continue
            log.info("::: Gathering sample info from MGnify...")
            query = urllib.parse.quote(f"biomes/{biome}/samples")
            if args.mg_filter is not None:
                samples = [
                    get_mgnify_data(x)
                    for x in tqdm.tqdm(
                        session.iterate(query, api_filter),
                        desc="Samples processed",
                        leave=False,
                    )
                ]
            else:
                samples = [
                    get_mgnify_data(x)
                    for x in tqdm.tqdm(
                        session.iterate(query),
                        desc="Samples processed",
                        leave=False,
                    )
                ]
            if len(samples) > 0:
                samples = concat_df(samples)

                log.info("::: Getting ENA links...")

                p = Pool(args.threads)
                dfs = list(
                    tqdm.tqdm(
                        p.imap_unordered(
                            partial(get_ena_data, filter_conditions=args.ena_filter),
                            samples.accession,
                        ),
                        total=len(samples.accession),
                        leave=False,
                        ncols=80,
                        desc="Samples processed",
                    )
                )
                p.close()
                p.join()
                # print(dfs)
                dfs = [x for x in dfs if x is not None]
                if len(dfs) > 0:
                    dfs = concat_df(dfs)
                    dfs = samples.merge(dfs, on="sample_accession")
                else:
                    columns = [
                        "study_accession",
                        "experiment_accession",
                        "run_accession",
                        "read_count",
                        "instrument_model",
                        "instrument_platform",
                        "library_layout",
                        "library_strategy",
                        "library_selection",
                        "fastq_ftp",
                    ]
                    dfs[columns] = None

                dfs["query_biome"] = biome

                if args.exclude_terms is not None:
                    if len(args.exclude_terms) > 1:
                        terms = "|".join(args.exclude_terms)
                    else:
                        terms = args.exclude_terms[0]

                    dfs = search_str(terms, dfs, invert=True)

                logging.info("::: Writing output file...")
                dfs.to_csv(
                    output_files[biome],
                    sep="\t",
                    index=False,
                )
            else:
                logging.info("No samples found for this biome, skipping")

    if args.combine:
        logging.info("Combining output files...")
        if os.path.exists(output_files["combined"]):
            os.remove(output_files["combined"])
        dfs = []
        for file in output_files.values():
            if file == output_files["combined"]:
                continue
            if os.path.exists(file):
                dfs.append(pd.read_csv(file, sep="\t"))
        dfs = concat_df(dfs)
        # A sample might belong to two query biomes (lineage), let's keep the one more detailed
        dfs["query_biome_n"] = dfs["query_biome"].str.split(":").str.len()
        dfs = dfs.loc[
            dfs.reset_index().groupby(["run_accession"])["query_biome_n"].idxmax()
        ]
        dfs = dfs.drop(columns=["query_biome_n"])
        dfs.drop_duplicates().to_csv(output_files["combined"], sep="\t", index=False)

    logging.info("ALL DONE.")
