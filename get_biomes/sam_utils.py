import pysam
import numpy as np
import re, os, sys
import pandas as pd
from multiprocessing import Pool
import functools
from scipy import stats
import pysam
import tqdm
import logging
import warnings
from get_biomes.utils import (
    is_debug,
    calc_chunksize,
    fast_flatten,
    initializer,
    concat_df,
    __from_nx_to_igraph,
)
import pyranges as pr
import datatable as dt
import networkx as nx
from grakel import GraphKernel, VertexHistogram
from functools import reduce
from scipy import stats

log = logging.getLogger("my_logger")

sys.setrecursionlimit(10**6)

# Function to calculate evenness of coverage
def coverage_evenness(coverage):
    """
    Calculate the evenness of coverage
    """
    # get coverage evenness
    # covEvenness = (
    #     (breadth * 100) / (100.0 * expBreadth * expBreadth) if expBreadth > 0 else 0
    # )
    # C = mean(X)
    # D2 = X[X<=C]
    # N = len(X)
    # n = len(D2)
    # E = 1 - (n - sum(D2) / C) / N
    C = float(round(np.mean(coverage)))
    D2 = [x for x in coverage if x <= C]
    if len(D2) == 0:  # pragma: no cover
        covEvenness = 1.0
    else:
        covEvenness = 1.0 - (len(D2) - sum(D2) / C) / len(coverage)

    return covEvenness


# function to calculate GC content
def calc_gc_content(seq):
    """Calculate GC content of a sequence

    Args:
        seq (str): DNA sequence

    Returns:
        int: Number of GC content
    """
    gc = seq.count("G") + seq.count("C")
    return gc


def create_pyranges(reference, starts, ends, strands):
    """[summary]

    Args:
        reference ([type]): [description]
        starts ([type]): [description]
        ends ([type]): [description]
        strands ([type]): [description]
    """
    chromosomes = [reference] * len(starts)
    chromosomes = pd.Series(chromosomes).astype("category")
    starts = pd.Series(starts)
    ends = pd.Series(ends)
    strands = pd.Series(strands).astype("category")

    return pr.PyRanges(
        pd.DataFrame(
            {"Chromosome": chromosomes, "Start": starts, "End": ends, "Strand": strands}
        )
    )


def get_bam_stats(params, ref_lengths=None):
    """
    Worker function per chromosome
    loop over a bam file and create tuple with lists containing metrics:
    two definitions of the edit distances to the reference genome scaled by aligned read length
    """
    bam, reference = params
    samfile = pysam.AlignmentFile(bam, "rb")
    edit_distances = []
    # edit_distances_md = []
    ani_nm = []
    # ani_md = []
    read_length = []
    read_aligned_length = []
    read_mapq = []
    read_aln_score = []
    read_names = []
    read_gc_content = []
    n_alns = 0
    if ref_lengths is None:
        reference_length = int(samfile.get_reference_length(reference))
        bam_reference_length = reference_length
    else:
        reference_length = int(ref_lengths.loc[reference, "length"])
        bam_reference_length = int(samfile.get_reference_length(reference))

    log.debug(f"Processing reference {reference}")
    log.debug(f"Reference length: {reference_length:,}")
    log.debug(f"BAM reference length: {bam_reference_length:,}")
    starts = []
    ends = []
    strands = []
    for aln in samfile.fetch(reference=reference, multiple_iterators=False):
        n_alns += 1
        if aln.has_tag("AS"):
            read_aln_score.append(aln.get_tag("AS"))
        else:
            read_aln_score.append(np.nan)

        if aln.has_tag("AS"):
            edit_distances.append(aln.get_tag("NM"))
            ani_nm.append((1 - ((aln.get_tag("NM") / aln.infer_query_length()))) * 100)
        else:
            edit_distances.append(np.nan)
            ani_nm.append(np.nan)

        read_gc_content.append(calc_gc_content(aln.query_sequence))
        read_length.append(aln.infer_read_length())
        read_aligned_length.append(aln.query_alignment_length)
        read_mapq.append(aln.mapping_quality)
        read_names.append(aln.query_name)
        # check if strand is reverse
        if aln.is_reverse:
            strand = "-"
        else:
            strand = "+"
        starts.append(aln.reference_start)
        ends.append(aln.reference_end)
        strands.append(strand)

    if n_alns > 1:
        # get bases covered by reads pileup
        cov_pos = [
            pileupcolumn.n
            for pileupcolumn in samfile.pileup(
                reference,
                start=None,
                stop=None,
                region=None,
                stepper="nofilter",
                min_mapping_quality=0,
            )
        ]
        # convert datafrane to pyranges
        ranges = create_pyranges(reference, starts, ends, strands)
        ranges = ranges.merge(strand=False).lengths().to_list()
        max_covered_bases = np.max(ranges)
        mean_covered_bases = np.mean(ranges)
        bases_covered = int(len(cov_pos))
        # get SD from covered bases
        cov_sd = np.std(cov_pos)
        # get average coverage
        mean_coverage = sum(cov_pos) / reference_length
        mean_coverage_covered = sum(cov_pos) / bases_covered

        breadth = bases_covered / reference_length
        exp_breadth = 1 - np.exp(-mean_coverage)
        breadth_exp_ratio = breadth / exp_breadth

        cov_evenness = coverage_evenness(cov_pos)
        gc_content = (np.sum(read_gc_content) / np.sum(read_length)) * 100
        c_v = cov_sd / mean_coverage
        read_mapq = [np.nan if x == 255 else x for x in read_mapq]

        log.debug(f"Number of reads: {len(set(read_names)):,}")
        log.debug(f"Number of alignments: {n_alns:,}")
        log.debug(f"Bases covered: {bases_covered:,}")
        log.debug(f"Mean coverage: {mean_coverage:.2f}")
        log.debug(f"Mean coverage covered: {mean_coverage_covered:.2f}")
        log.debug(f"Max covered bases: {max_covered_bases:,}")
        log.debug(f"Mean covered bases: {mean_covered_bases:.2f}")
        log.debug(f"SD: {cov_sd:.2f}")
        log.debug(f"Breadth: {breadth:.2f}")
        log.debug(f"Exp. breadth: {exp_breadth:.2f}")
        log.debug(f"Breadth/exp. ratio: {breadth_exp_ratio:.2f}")
        log.debug(f"Cov. evenness: {cov_evenness:.2f}")
        log.debug(f"C_v: {c_v:.2f}")
        log.debug(f"Mean mapq: {np.mean(read_mapq):.2f}")
        log.debug(f"GC content: {gc_content:.2f}")
        data = BamAlignment(
            reference=reference,
            n_alns=n_alns,
            reference_length=reference_length,
            bam_reference_length=bam_reference_length,
            mean_coverage=mean_coverage,
            mean_coverage_covered=mean_coverage_covered,
            bases_covered=bases_covered,
            max_covered_bases=max_covered_bases,
            mean_covered_bases=mean_covered_bases,
            cov_evenness=cov_evenness,
            breadth=breadth,
            exp_breadth=exp_breadth,
            breadth_exp_ratio=breadth_exp_ratio,
            c_v=c_v,
            edit_distances=edit_distances,
            # edit_distances_md=edit_distances_md,
            ani_nm=ani_nm,
            # ani_md=ani_md,
            read_length=read_length,
            read_gc_content=read_gc_content,
            read_aligned_length=read_aligned_length,
            mapping_quality=read_mapq,
            read_names=set(read_names),
            read_aln_score=read_aln_score,
        )
    else:
        data = BamAlignment(
            reference=reference,
            n_alns=n_alns,
            reference_length=reference_length,
            bam_reference_length=bam_reference_length,
            mean_coverage=np.nan,
            mean_coverage_covered=np.nan,
            bases_covered=np.nan,
            max_covered_bases=np.nan,
            mean_covered_bases=np.nan,
            cov_evenness=np.nan,
            breadth=np.nan,
            exp_breadth=np.nan,
            breadth_exp_ratio=np.nan,
            c_v=np.nan,
            edit_distances=np.nan,
            # edit_distances_md=np.nan,
            ani_nm=np.nan,
            # ani_md=np.nan,
            read_length=np.nan,
            read_gc_content=np.nan,
            read_aligned_length=np.nan,
            mapping_quality=np.nan,
            read_names=read_names,
            read_aln_score=read_aln_score,
        )
    return data


class BamAlignment:
    """
    Class to store alignment information
    """

    def __init__(
        self,
        reference,
        n_alns,
        read_length,
        read_gc_content,
        read_aligned_length,
        mapping_quality,
        edit_distances,
        # edit_distances_md,
        ani_nm,
        # ani_md,
        bases_covered,
        max_covered_bases,
        mean_covered_bases,
        mean_coverage,
        mean_coverage_covered,
        reference_length,
        bam_reference_length,
        breadth,
        exp_breadth,
        breadth_exp_ratio,
        c_v,
        cov_evenness,
        read_names,
        read_aln_score,
    ):
        self.reference = reference
        self.n_alns = n_alns
        self.read_length = read_length
        self.read_gc_content = read_gc_content
        self.read_aligned_length = read_aligned_length
        self.read_aln_score = read_aln_score
        self.mapping_quality = mapping_quality
        self.edit_distances = edit_distances
        # self.edit_distances_md = edit_distances_md
        self.ani_nm = ani_nm
        # self.ani_md = ani_md
        self.bases_covered = bases_covered
        self.max_covered_bases = max_covered_bases
        self.mean_covered_bases = mean_covered_bases
        self.mean_coverage = mean_coverage
        self.mean_coverage_covered = mean_coverage_covered
        self.reference_length = reference_length
        self.bam_reference_length = bam_reference_length
        self.breadth = breadth
        self.exp_breadth = exp_breadth
        self.breadth_exp_ratio = breadth_exp_ratio
        self.c_v = c_v
        self.cov_evenness = cov_evenness
        self.read_names = read_names
        # function to convert class to dict

    def as_dict(self):
        return {
            "reference": self.reference,
            "n_reads": self.read_names,
            "n_alns": self.n_alns,
            "read_length": self.read_length,
            "read_gc_content": self.read_gc_content,
            "read_aligned_length": self.read_aligned_length,
            "read_aln_score": self.read_aln_score,
            "mapping_quality": self.mapping_quality,
            "edit_distances": self.edit_distances,
            # "edit_distances_md": np.mean(self.edit_distances_md),
            "ani_nm": self.ani_nm,
            # "ani_md": np.mean(self.ani_md),
            "bases_covered": self.bases_covered,
            "max_covered_bases": self.max_covered_bases,
            "mean_covered_bases": self.mean_covered_bases,
            "mean_coverage": self.mean_coverage,
            "mean_coverage_covered": self.mean_coverage_covered,
            "reference_length": self.reference_length,
            "breadth": self.breadth,
            "exp_breadth": self.exp_breadth,
            "breadth_exp_ratio": self.breadth_exp_ratio,
            "c_v": self.c_v,
            "cov_evenness": self.cov_evenness,
        }

    def to_summary(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            read_length_mean = np.mean(self.read_length)
            read_length_median = np.median(self.read_length)
            read_length_mode = stats.mode(self.read_length)[0][0]
            read_aligned_length = np.mean(self.read_aligned_length)
            read_aln_score = np.mean(self.read_aln_score)
            mapping_quality = np.mean(self.mapping_quality)
            edit_distances = np.mean(self.edit_distances)
            read_ani_mean = np.mean(self.ani_nm)
            gc_content = (np.sum(self.read_gc_content) / np.sum(self.read_length)) * 100

        return {
            "reference": self.reference,
            "n_reads": len(self.read_names),
            "n_alns": self.n_alns,
            "read_length_mean": read_length_mean,
            "read_length_median": read_length_median,
            "read_length_mode": read_length_mode,
            "gc_content": gc_content,
            "read_aligned_length": read_aligned_length,
            "read_aln_score": read_aln_score,
            "mapping_quality": mapping_quality,
            "edit_distances": edit_distances,
            # "edit_distances_md": np.mean(self.edit_distances_md),
            "read_ani_mean": read_ani_mean,
            # "ani_md": np.mean(self.ani_md),
            "bases_covered": self.bases_covered,
            "max_covered_bases": self.max_covered_bases,
            "mean_covered_bases": self.mean_covered_bases,
            "coverage_mean": self.mean_coverage,
            "coverage_covered_mean": self.mean_coverage_covered,
            "reference_length": self.reference_length,
            "bam_reference_length": self.bam_reference_length,
            "breadth": self.breadth,
            "exp_breadth": self.exp_breadth,
            "breadth_exp_ratio": self.breadth_exp_ratio,
            "c_v": self.c_v,
            "cov_evenness": self.cov_evenness,
        }


# Inspired from https://gigabaseorgigabyte.wordpress.com/2017/04/14/getting-the-edit-distance-from-a-bam-alignment-a-journey/
def process_bam(bam, read_taxonomy, threads=1):
    """
    Processing function: calls pool of worker functions
    to extract from a bam file two definitions of the edit distances to the reference genome scaled by read length
    Returned in a pandas DataFrame
    """
    dt.options.progress.clear_on_success = True
    dt.options.nthreads = threads
    # logging.info(f"Loading BAM file")
    # save = pysam.set_verbosity(0)
    # samfile = pysam.AlignmentFile(bam, "rb")
    # chr_lengths = []
    # for chrom in samfile.references:
    #     chr_lengths.append(samfile.get_reference_length(chrom))
    # max_chr_length = np.max(chr_lengths)
    # pysam.set_verbosity(save)
    # ref_lengths = None

    # logging.info(f"Loading read taxonomic labels")
    # read_taxonomies = dt.fread(read_taxonomy, nthreads=threads)
    # read_taxonomies.names = ["name", "taxonomy"]
    # nrows = read_taxonomies.shape[0]
    # logging.info(f"Loaded {read_taxonomies.shape[0]:,} read taxonomic labels")
    # logging.info(f"Doing some sanity checks on read taxonomy file")
    # read_taxonomies[:, dt.update(duplicated=(dt.count() > 1)), dt.by(dt.f.name)]
    # if read_taxonomies[dt.f.duplicated > 0, :].shape[0] > 0:
    #     print(read_taxonomies[dt.f.duplicated > 0, :])
    #     logging.error(
    #         f"{read_taxonomies[dt.f.duplicated > 0,:].shape[0]:,} read taxonomies are duplicated"
    #     )
    #     exit(1)
    # del read_taxonomies[:, ["duplicated"]]

    # read_taxonomies = read_taxonomies[
    #     (dt.f.taxonomy != "root")
    #     & (dt.f.taxonomy != "Root")
    #     & (dt.f.taxonomy != "ROOT"),
    #     :,
    # ]
    # if read_taxonomies.shape[0] != nrows:
    #     logging.info(
    #         f"::: {nrows - read_taxonomies.shape[0]:,} removed as they are classified as root"
    #     )
    # if not samfile.has_index():
    #     logging.info(f"BAM index not found. Indexing...")
    #     if max_chr_length > 536870912:
    #         logging.info(f"A reference is longer than 2^29, indexing with csi")
    #         pysam.index(bam, "-c")
    #     else:
    #         pysam.index(bam)

    #     logging.info(f"Reloading BAM file")
    #     samfile = pysam.AlignmentFile(
    #         bam, "rb"
    #     )  # Need to reload the samfile after creating index

    # total_refs = samfile.nreferences
    # logging.info(f"Found {total_refs:,} reference sequences")
    # logging.info(f"Found {samfile.mapped:,} alignments")
    # logging.info(f"Extracting aligned reads from BAM file...")
    # reads = []
    # with pysam.AlignmentFile(bam) as bam_file:
    #     for entry in tqdm.tqdm(
    #         bam_file,
    #         total=samfile.mapped,
    #         leave=False,
    #         ncols=80,
    #         desc=f"Alignments processed",
    #     ):
    #         reads.append(
    #             dt.Frame(
    #                 [
    #                     (
    #                         entry.reference_name,
    #                         bam_file.get_reference_length(entry.reference_name),
    #                         entry.query_name,
    #                         entry.query_length,
    #                         entry.reference_start,
    #                         entry.reference_end,
    #                     )
    #                 ],
    #                 names=[
    #                     "reference",
    #                     "reference_length",
    #                     "name",
    #                     "read_length",
    #                     "start",
    #                     "end",
    #                 ],
    #             )
    #         )
    # # references = samfile.references
    # # params = zip([bam] * len(references), references)
    # # try:
    # #     if is_debug():
    # #         reads = list(map(get_alns, params))
    # #     else:

    # #         p = Pool(threads)
    # #         c_size = calc_chunksize(threads, len(references))
    # #         reads = list(
    # #             tqdm.tqdm(
    # #                 p.imap_unordered(get_alns, params, chunksize=c_size),
    # #                 total=len(references),
    # #                 leave=False,
    # #                 ncols=80,
    # #                 desc=f"References processed",
    # #             )
    # #         )

    # #         p.close()
    # #         p.join()

    # # except KeyboardInterrupt:
    # #     logging.info(f"User canceled the operation. Terminating jobs.")
    # #     p.terminate()
    # #     p.join()
    # #     sys.exit()
    # # # reads = dt.Frame(dict(name=np.unique(reads)))
    # # # reads = reads[:, dt.f[:].extend({"classified": "withTax"})]
    # reads = dt.rbind(reads)
    # logging.info(f"Found {reads.shape[0]:,} reads")

    # # do we lose references?
    # logging.info(f"Removing reads without taxonomic labels")
    # read_taxonomies.key = "name"
    # reads = reads[
    #     :,
    #     :,
    #     dt.join(
    #         read_taxonomies,
    #     ),
    # ][dt.f.taxonomy != None, :].to_pandas()
    # # write table to tsv file
    # logging.info(f"Writing reads to file")
    # reads.to_csv("reads.tsv.gz", sep="\t", compression="gzip", index=False)
    # exit()

    reads = dt.fread("reads.tsv.gz", sep="\t", nthreads=threads).to_pandas()
    logging.info(f"Kept {reads.shape[0]:,} aligned reads with taxonomic labels")
    dfGrouped = reads.groupby("reference")
    logging.info(f"Found {dfGrouped.ngroups:,} references with taxonomies")
    logging.info(f"Creating taxonomic graphs")
    for name, group in dfGrouped:
        # print(
        #     group.groupby(["taxonomy"])
        #     .size()
        #     .reset_index(name="counts")
        #     .sort_values("counts", ascending=False)
        # )
        # print(group.sort_values("start", ascending=True))

        # print(list(sliding_window_iter(group["start"], window=250, overlap=50)))
        pos = np.arange(0, np.unique(group["reference_length"]))
        from collections import defaultdict

        gen_data = defaultdict(dict)
        for a, q, b, c, d in tqdm.tqdm(
            zip(
                group["reference"],
                group["name"],
                group["start"],
                group["end"],
                group["reference_length"],
            ),
            total=group.shape[0],
            leave=False,
            ncols=100,
            desc=f"Alignments processed",
        ):
            b = b - 1

            if a in gen_data:
                gen_data[a]["cov"][b:c] += 1
                gen_data[a]["n_alns"] += 1
            else:
                gen_data[a]["cov"] = np.zeros(d, dtype=int)
                gen_data[a]["cov"][b:c] += 1
                gen_data[a]["n_alns"] = 0
                gen_data[a]["n_alns"] += 1
        mean_cov = np.mean(gen_data[name]["cov"])
        max_chunk_sum = mean_cov
        lst = gen_data[name]["cov"].tolist()
        result = [[lst.pop(0)]]

        for item in lst:
            l = result[-1]
            l.append(item)
            if np.mean(l) > max_chunk_sum:
                result.append([item])
            else:
                result[-1].append(item)
        r_len = []
        for r in result:
            r_len.append(len(r))

        window = np.max(r_len)
        if window > np.unique(group["reference_length"]):
            window = np.unique(group["reference_length"])[0]

        print(window)
        print(mean_cov)
        # window = 500
        overlap = round_half_up(np.mean(group["read_length"]))
        pos = ovl_window(pos, window=window, overlap=overlap)
        p_min = pos[pos.shape[0] - 1, 0]
        p_max = pos[pos.shape[0] - 1, pos.shape[1] - 1]
        pos = [pos[i, :].tolist() for i in range(pos.shape[0])]
        if p_max < np.unique(group["reference_length"]):
            pos_row = np.arange(
                p_min + overlap, np.unique(group["reference_length"])
            ).tolist()
            pos.append(pos_row)
        # print(pos1)
        # n = 250
        # m = 50
        # pos = [pos[i : i + n] for i in range(0, len(pos), n - m)]
        graphs = []
        n_paths = []
        for p in pos:
            gs = []
            dfs = []

            # extract all reads with positions in the interval
            for i in p:
                idx = pd.IntervalIndex.from_arrays(group["start"], group["end"])
                dfs.append(group[idx.contains(i)].copy())
            group_sub = concat_df(dfs).drop_duplicates()

            group_sub["start_w"] = min(p)
            group_sub["end_w"] = max(p)

            tax = group_sub["taxonomy"].to_list()
            root_row = pd.DataFrame({"source": ["root"], "target": ["root"]})

            if len(tax) > 0:
                for t in tax:
                    ts = t.split(";")
                    if ts is None:
                        ts = t
                    if "root" in ts:
                        res = list(zip(ts, ts[1:]))
                        res = pd.DataFrame(res, columns=["source", "target"])
                        res = pd.concat([root_row, res])
                        gs.append(res)
                    else:
                        ts.insert(0, "root")
                        res = list(zip(ts, ts[1:]))
                        res = pd.DataFrame(res, columns=["source", "target"])
                        res = pd.concat([root_row, res])
                        gs.append(res)

                if len(gs) > 1:
                    gs = (
                        concat_df(gs)
                        .groupby(["source", "target"])
                        .size()
                        .reset_index(name="weight")
                    )
                    # remove all taxa with only one read
                    # TODO: make something smarter to remove noise
                    gs_c = gs.copy()
                    gs_c = gs_c[gs_c["weight"] > 1]
                    if len(gs_c.index) > 0:
                        gs = gs_c
                    else:
                        gs = gs
                else:
                    gs = gs[0]
                    gs["weight"] = 1

                G = nx.from_pandas_edgelist(
                    gs,
                    edge_attr=["weight"],
                    create_using=nx.DiGraph(),
                )
                G.remove_edges_from(nx.selfloop_edges(G))
                paths = []

                for node in G:
                    if G.out_degree(node) == 0:  # it's a leaf
                        p = nx.shortest_path(G, "root", node)
                        w = path_cost(G, p)
                        paths.append({"p": p, "w": w, "w": w})
                g_w = []
                for path in paths:
                    g_w.append(path["w"])
                print(g_w)

                zsc = stats.zscore(g_w / np.sum(g_w))

                median_y = np.median(g_w)
                median_absolute_deviation_y = np.median(
                    [np.abs(y - median_y) for y in g_w]
                )
                modified_z_scores = [
                    0.6745 * (y - median_y) / median_absolute_deviation_y for y in g_w
                ]
                for path in paths:
                    path["z"] = modified_z_scores[paths.index(path)]

                # g_w = np.sum(g_w)
                from scipy import stats

                longest_induced_path = [
                    path["p"] for path in paths if abs(path["z"]) < 3.0
                ]
                print("#######")
                print(paths)
                print("longest")
                longest_induced_path = np.unique(fast_flatten(longest_induced_path))
                print(longest_induced_path)
                # longest_induced_path = max(paths, key=lambda x: x["w"])
                G = G.subgraph(longest_induced_path)

                G = nx.relabel_nodes(G, {k: v for k, v in enumerate(tax)})
                labels = {v: v for k, v in enumerate(G.nodes)}
                nx.set_node_attributes(G, labels, "label")
                g = __from_nx_to_igraph(G, directed=True)
                graphs.append(G)
                # n_paths.append(";".join(longest_induced_path]))
            else:
                graphs.append(None)
                n_paths.append(None)
        # for g in graphs:
        #     print(g)
        from grakel.utils import graph_from_networkx

        G = graph_from_networkx(graphs, node_labels_tag="label")
        from grakel.kernels import ShortestPath

        gk = VertexHistogram(normalize=True)
        K = gk.fit_transform(G)
        # convert distance matrix to pandas dataframe of pairs
        K = pd.DataFrame(K)
        K = matrix_to_xy(K, reset_index=True)

        K["from"] = K["from"].astype(str)
        K["to"] = K["to"].astype(str)
        K = K[K["from"] != K["to"]]
        # import natsort

        # print(
        #     K.groupby("row", as_index=False)
        #     .mean()
        #     .sort_values("row", ascending=True, key=natsort.natsort_keygen())
        # )
        # exit()
        n = len(pos)

        # # create from 0 to n-1
        coords = pd.DataFrame(
            {
                "from": np.insert(np.arange(n), 0, 0),
                "to": np.append(np.arange(n), n - 1),
            }
        )

        coords["start"] = [pos[i][0] for i in coords["to"]]
        coords["end"] = [pos[i][-1] for i in coords["to"]]
        coords["from"] = coords["from"].astype(str)
        coords["to"] = coords["to"].astype(str)
        # long_form = K.unstack()
        # long_form.index.rename(["from", "to"], inplace=True)
        # long_form = long_form.to_frame("similarity").reset_index()
        K["reference"] = name
        K["avg_sim"] = np.mean(K["similarity"])
        # is a breakpoint if similarity is smaller than avg_sim
        K["is_breakpoint"] = K["similarity"] < K["avg_sim"]
        # g_taxs = pd.DataFrame(n_paths, columns=["gtax"])
        # g_taxs["to"] = np.arange(len(g_taxs))
        # long_form = pd.merge(long_form, g_taxs, on="to")
        # inner join coords and long_form
        K = pd.merge(coords, K, on=["from", "to"])
        print(K)
    exit()

    try:
        logging.info(f"Getting stats for each reference...")

        if is_debug():
            data = list(
                map(functools.partial(get_bam_stats, ref_lengths=ref_lengths), params)
            )
        else:

            p = Pool(
                threads, initializer=initializer, initargs=([params, ref_lengths],)
            )
            c_size = calc_chunksize(
                n_workers=threads,
                len_iterable=len(references),
            )
            data = list(
                tqdm.tqdm(
                    p.imap_unordered(
                        functools.partial(get_bam_stats, ref_lengths=ref_lengths),
                        params,
                        chunksize=c_size,
                    ),
                    total=len(references),
                    leave=False,
                    ncols=80,
                    desc=f"References processed",
                )
            )

            p.close()
            p.join()

    except KeyboardInterrupt:
        logging.info(f"User canceled the operation. Terminating jobs.")
        p.terminate()
        p.join()
        sys.exit()
    return data


def sliding_window_iter(series, window, overlap):
    """series is a column of a dataframe"""
    for start_row in range(len(series) - window + overlap):
        yield series[start_row : start_row + window]


def get_alns(params):
    bam, reference = params
    samfile = pysam.AlignmentFile(bam, "rb")
    alns = []
    for entry in samfile.fetch(reference=reference, multiple_iterators=True):
        alns.append(
            dt.Frame(
                [
                    (
                        entry.reference_name,
                        reference,
                        entry.query_name,
                        entry.query_length,
                        entry.reference_start,
                        entry.reference_end,
                    )
                ],
                names=[
                    "reference",
                    "reference_length",
                    "name",
                    "read_length",
                    "start",
                    "end",
                ],
            )
        )
        reads = dt.rbind(alns)
    return reads


# from https://stackoverflow.com/a/45730836/15704171
def ovl_window(a, window=4, overlap=2, copy=False):
    sh = (a.size - window + 1, window)
    st = a.strides * 2
    view = np.lib.stride_tricks.as_strided(a, strides=st, shape=sh)[0::overlap]
    if copy:
        return view.copy()
    else:
        return view


def matrix_to_xy(df, columns=None, reset_index=False):
    bool_index = np.triu(np.ones(df.shape)).astype(bool)
    xy = (
        df.where(bool_index).stack().reset_index()
        if reset_index
        else df.where(bool_index).stack()
    )
    if reset_index:
        xy.columns = columns or ["from", "to", "similarity"]
    return xy


# From https://stackoverflow.com/a/69655565/15704171
def round_half_up(x):
    round_lambda = lambda z: (int(z > 0) - int(z < 0)) * int(abs(z) + 0.5)
    if isinstance(x, (np.ndarray, np.generic)):
        return np.vectorize(round_lambda)(x)
    else:
        return round_lambda(x)


# From https://stackoverflow.com/a/56665813/15704171
def path_cost(G, path):
    """
    Compute the cost of a path in a graph.
    There's a penalty the deeper we get in the DAG to account for the
    uncertainty of the tax annotations.

        Parameters
        ----------
        G : networkx.Graph
            A graph.
        path : list
            A list of nodes in the graph.

        Returns
        -------
        cost : float
            The cost of the path.
    """
    return sum(
        [G[path[i]][path[i + 1]]["weight"] / (i + 1) for i in range(len(path) - 1)]
    )
