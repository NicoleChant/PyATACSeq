import os
import csv
from collections import defaultdict
from pathlib import Path
from termcolor import colored
import sys
import pandas as pd

def load_cell_types() -> dict:
    cell_types = defaultdict(list)
    with open("ATAC_bam_files.txt", mode="r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            cell_types[row['cell_type'], 'bam'].append(row['bam_file'])
            cell_types[row['cell_type'], 'bigBed'].append(row['big_bed_file'])
    return cell_types

if __name__ == "__main__":

    comparisons = []
    with open("cell_lines.txt") as f:
        for line in f:
            line = line.strip()
            cellA, cellB = line.split(",")
            comparisons.append((cellA, cellB))

    print(f"TOTAL COMPARISONS: {len(comparisons)}.")

    cell_types = load_cell_types()

    results_dir_global = Path("ATAC_RESULTS").resolve()
    results_dir_global.mkdir(exist_ok=True)
    total_p = 140

    for cellA, cellB in comparisons:

        # destination directory
        results_dir = results_dir_global.joinpath(f"{cellA}_vs_{cellB}")
        print(colored(f"OUTSOURCING RESULTS TO ----> {results_dir}...", "red"))
        results_dir.mkdir(exist_ok=True)
        # expected output
        counts_result = Path(f"{results_dir}/counts.{cellA}_vs_{cellB}.tsv")
        # if expected output exists skip please
        # if counts_result.is_file():
        #    print(colored(f"Skipping CELL {cellA} vs. CELL {cellB} because output already exists!", "red"))
        #    continue
        print(f"Processing cell {cellA} vs. cell {cellB}...")
        
        # merge bed files to peaks
        # STEP 1: Merge peaks for cell type A using both replicas
        print(f"Step 1: Merging peaks from {cellA}...")
        replica_A, replica_B = cell_types[cellA, 'bigBed']
        replica_A = f"BIGBED_FILES_ATACSeq/{replica_A}.bed"
        replica_B = f"BIGBED_FILES_ATACSeq/{replica_B}.bed"
        os.system(f"cat {replica_A} {replica_B} | sort -k1,1 -k2,2n | bedtools merge -i - > {results_dir}/{cellA}.merged_peaks.bed")
        print("=" * total_p)
        print(colored(f"Step 1: Merging peaks from {cellA}. Status: COMPLETED.\n", "green"))

        # STEP 2: Merged peaks for cell type B using both replicas
        print(f"Step 2: Merging peaks from {cellB}...")
        replica_A, replica_B = cell_types[cellB, 'bigBed']
        replica_A = f"BIGBED_FILES_ATACSeq/{replica_A}.bed"
        replica_B = f"BIGBED_FILES_ATACSeq/{replica_B}.bed"
        os.system(f"cat {replica_A} {replica_B} | sort -k1,1 -k2,2n | bedtools merge -i - > {results_dir}/{cellB}.merged_peaks.bed")
        print("=" * total_p)
        print(colored(f"Step 2: Merging peaks from {cellB}. Status: COMPLETED.\n", "green"))


        # STEP 3: Create union set from merged peaks of both cell types
        print(f"Step 3: Unionizing merged peaks from cellline {cellA} and cellline {cellB}...")
        union_peaks = f"{results_dir}/merged_peaks_{cellA}_vs_{cellB}.union_peaks.bed"
        os.system(f"cat {results_dir}/{cellA}.merged_peaks.bed {results_dir}/{cellB}.merged_peaks.bed | sort -k1,1 -k2,2n | bedtools merge -i - > {union_peaks}")
        print("=" * total_p)
        print(colored(f"Step 3: Unionizing merged peaks from cell line {cellA} and cell line {cellB}. Status: COMPLETED.\n", "green"))

        # STEP 4: Process bam files from both replicas using the unionized peaks from both cell lines
        # and create a count MATRIX which will serve as input for the DESEQ analysis.
        print("Step 4: Extracting counts on the unionized peaks from both cellines using the bam files...")
        bam_files = ' '.join(map(lambda x: f'BAM_FILES_ATACSeq/{x}.bam', cell_types[cellA, 'bam'] + cell_types[cellB, 'bam']))

        # create index if they don't exist
        for bam in bam_files.split(' '):
            bam_index = Path(str(bam) + ".bai")
            if not bam_index.is_file():
                os.system(f"samtools index {bam}")
                print(colored(f"CREATED INDEX FILE FOR: {bam}!", "blue"))

        print(f"Using bam files:\n{bam_files}\n")
        if not counts_result.is_file():
            os.system(f"bedtools multicov -bams {bam_files} -bed {union_peaks} > {counts_result}")
        else:
            print(colored(f"Counts file '{counts_result}' from bam files already exist!", "red"))

        print("=" * total_p)
        print(colored("Step 4: Extracting counts on the unionized peaks from both cell lines using the bam files. Status: COMPLETED.\n", "green"))

        # feature counts ; formatting (prepare for DESEQ2)
        deseq_input = Path("DESEQ_INPUT").resolve()
        deseq_input.mkdir(exist_ok=True)
        deseq_input = deseq_input.joinpath(f"counts_deseq.{cellA}_vs_{cellB}.tsv")
  
        # reformating for deseq
        bam = cell_types[cellA, 'bam'] + cell_types[cellB, 'bam']
        results_df = pd.read_table(counts_result, 
                                   header=None, 
                                   names=["seqID", "start", "end"] + bam
                                   )
        results_df.loc[:, "coords"] = results_df["seqID"].astype(str) + "_" \
                                      + results_df["start"].astype(str) + "_" \
                                      + results_df["end"].astype(str)

        results_df[["coords"] + bam].to_csv(deseq_input, 
                                                  sep="\t", 
                                                  mode="w", 
                                                  index=False)
        # << END
        print("=" * total_p)
        print(colored(f"Processing CELL {cellA} vs. CELL {cellB}. Status: COMPLETED.\n", "green"))
        print(colored(f"RESULTS CAN BE FOUND AT ----> {results_dir}.", "red"))

        input()
