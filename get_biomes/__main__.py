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
    API_BASE,
    get_mgnify_data,
    get_ena_data,
    create_output_files,
)
from jsonapi_client import Filter, Session
from urllib.parse import urlencode
import urllib.parse
import tqdm
from multiprocessing import Pool
import os
import hashlib
from pathlib import Path
from functools import partial

log = logging.getLogger("my_logger")


def main():

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
        params = args.mg_filter

        api_filter = Filter(urlencode(params))

        for biome in biomes:
            log.info(f"Processing biome: {biome}")
            if Path(output_files[biome]).is_file():
                log.info(f"File {output_files[biome]} already exists, skipping")
                continue
            log.info("::: Gathering sample info from MGnify...")
            query = urllib.parse.quote(f"biomes/{biome}/samples")
            samples = [
                get_mgnify_data(x)
                for x in tqdm.tqdm(
                    session.iterate(query, api_filter),
                    desc="Samples processed",
                    leave=False,
                )
            ]
            samples = concat_df(samples)

            log.info("::: Getting ENA links...")

            p = Pool(args.threads)
            dfs = list(
                tqdm.tqdm(
                    p.imap_unordered(
                        partial(get_ena_data, min_read_count=args.min_read_count),
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
            dfs = concat_df(dfs)
            dfs["query_biome"] = biome
            dfs = dfs.merge(samples, on="sample_accession")
            logging.info("::: Writing output file...")
            dfs.to_csv(
                output_files[biome],
                sep="\t",
                index=False,
            )

    if args.combine:
        logging.info("Combining output files...")
        dfs = []
        for file in output_files.values():
            if file == output_files["combined"]:
                continue
            dfs.append(pd.read_csv(file, sep="\t"))
        dfs = concat_df(dfs)
        dfs.to_csv(output_files["combined"], sep="\t", index=False)

    logging.info("ALL DONE.")


if __name__ == "__main__":
    main()
