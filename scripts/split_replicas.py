from collections import defaultdict
import os

if __name__ == "__main__":

    # fetch cell lines to study
    cell_lines_to_study = set()
    with open("cell_lines.txt") as f:
        for line in f:
            line = line.strip()
            cellA, cellB = line.split(",")
            cell_lines_to_study.add(cellA)
            cell_lines_to_study.add(cellB)


    # parse encode IDS
    bam_files = defaultdict(list)
    big_bed_files = defaultdict(list)
    flag = False
    with open("encode_ids.txt", mode="r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            colA, colB = line.split(" ")
            if colB.startswith("(") and colB == "(ATAC-seq)":
                flag = True
                cell_type = colA
                continue
            elif colB.startswith("("):
                flag = False

            if flag:
                bam_files[cell_type].append(colA)
                big_bed_files[cell_type].append(colB)

    # download link sample
    # sample: https://www.encodeproject.org/files/<id>/@@download/<id>.<filetag>

    with open("ATAC_bam_files.txt", mode="w", encoding="utf-8") as f:
        f.write("cell_type,bam_file,bam_file_download_link,big_bed_file,big_bed_download_link\n")
        for cell_type in bam_files:
            for bam, big_bed in zip(bam_files[cell_type], big_bed_files[cell_type]):
                bam_file_download_link = f"https://www.encodeproject.org/files/{bam}/@@download/{bam}.bam"
                big_bed_download_link = f"https://www.encodeproject.org/files/{big_bed}/@@download/{big_bed}.bigBed"
                f.write(f"{cell_type},{bam},{bam_file_download_link},{big_bed},{big_bed_download_link}\n")

                if cell_type in cell_lines_to_study:
                    # BAM
                    print(f"Downloading BAM file '{bam}' for cell line {cell_type}...")
                    os.system(f"wget -P BAM_FILES_ATACSeq {bam_file_download_link}")

                    # BIG BED
                    print(f"Downloading BIGBED file '{big_bed}' for cell type {cell_type}...")
                    os.system(f"wget -P BIGBED_FILES_ATACSeq {big_bed_download_link}")