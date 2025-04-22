import os
import pickle as pkl
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data
import pandas as pd
import argparse
from pathlib import Path
import csv
from collections import defaultdict
from pybedtools import BedTool
import pybedtools

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--cellA", type=str, default="HSC")
    parser.add_argument("--cellB", type=str, default="CMP")
    parser.add_argument("--gff", type=str, default="")
    parser.add_argument("--MAX_DISTANCE", type=int, default=500)
    parser.add_argument("--use_genes", type=int, default=1, choices=[0, 1])
    # parser.add_argument("--genes", type=str, default="mus_musculus_gene.csv")

    args = parser.parse_args()
    cellA = args.cellA
    cellB = args.cellB

    # Maximum Distance from a given genic region to filter
    MAX_DISTANCE = args.MAX_DISTANCE

    # if this is off we will run differential accessibility on regions and not do peak calling on genes
    use_genes = args.use_genes

    counts_file = Path(f"DESEQ_INPUT/counts_deseq.{cellA}_vs_{cellB}.tsv")
    if not counts_file.is_file():
        raise FileNotFoundError(f"Counts file {counts_file} does not exist!")

    metadata = []
    counts = defaultdict(int)
    names = []
    with open("ATAC_bam_files.txt", mode="r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            cell_line = row["cell_type"]
            if cell_line == cellA or cell_line == cellB:
                counts[cell_line] += 1
                names.append(row["bam_file"])
                metadata.append({
                        "sample": row["bam_file"],
                        "condition": "A" if cell_line == cellA else "B",
                        "group": cell_line
                        })
    metadata = pd.DataFrame(metadata).set_index("sample")
    counts_df = pd.read_table(counts_file)
    counts_df_unravel_coordinates = counts_df["coords"].str.split("_", expand=True)
    counts_df_unravel_coordinates.columns = ["seqID", "start", "end"]
    counts_df_unravel_coordinates.loc[:, "seqID"] = counts_df_unravel_coordinates["seqID"].apply(lambda x: x.split("r")[1])
    counts_df_unravel_coordinates_bed = BedTool.from_dataframe(counts_df_unravel_coordinates).sort()

    # map regions to closest genes
    # we do not really need gene names because they provide them in the GFF file
    # genes_names_df = pd.read_csv(args.genes)
    if use_genes:
        genes_bed = BedTool("Mus_musculus.GRCm38.101.gff3")\
                    .filter(lambda x: x[2] == "gene") # .saveas("only_genes.gff")

        closest_genes = pd.read_table(
                                counts_df_unravel_coordinates_bed.closest(genes_bed, k=1, d=True).fn,
                                header=None,
                                names=list(counts_df_unravel_coordinates.columns) \
                                    + ["chromosome", 
                                       "source", 
                                       "compartment", 
                                       "gene_start", 
                                       "gene_end", 
                                       "score", 
                                       "strand", 
                                       "phase", 
                                       "attributes", 
                                       "distance"]
                        )\
                        .query(f"abs(distance) < {MAX_DISTANCE}")

        closest_genes.loc[:, "gene_name"] = closest_genes["attributes"].apply(lambda x: x.split("Name=")[1].split(";")[0])
        closest_genes.loc[:, "coords"] = "chr" + closest_genes["seqID"].astype(str) + "_" \
                                            + closest_genes["start"].astype(str) + "_" \
                                            + closest_genes["end"].astype(str)
        counts_df_unravel_coordinates_merged = closest_genes.merge(
                                                                counts_df,
                                                                left_on="coords",
                                                                right_on="coords",
                                                                how="inner")\
                                            [["coords", "gene_name"] + names]\
                                            .groupby("gene_name")\
                                            .agg({replica: "sum" for replica in names})
        counts_df_unravel_coordinates_merged = counts_df_unravel_coordinates_merged.T
    else:
        counts_df_unravel_coordinates_merged = counts_df.set_index("coords").T

    inference = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(
        counts=counts_df_unravel_coordinates_merged,
        metadata=metadata,
        design="~condition", 
        # column ("B" vs "A")
        refit_cooks=True,
        inference=inference,
    )
    dds.deseq2()

    # Get stats
    stat_res = DeseqStats(
                    dds,  
                    contrast=np.array([0, 1]),
                    alpha=0.05,
                    cooks_filter=True,
                    independent_filter=True,)
    stat_res.summary()

    # save results
    OUTPUT_PATH = Path("DESEQ_OUTPUT_CHROMATIN").resolve()
    OUTPUT_PATH.mkdir(exist_ok=True)
    with open(os.path.join(OUTPUT_PATH, f"{cellA}_vs_{cellB}_stat_results_detailed_pipe.pkl"), "wb") as f:
        pkl.dump(dds, f)

    stat_res.results_df.to_csv(f"{OUTPUT_PATH}/{cellA}_vs_{cellB}_stat_results.deseq_out.txt",
                               sep="\t",
                               mode="w",
                               index=True,
                               header=True)
    # breakpoint()
