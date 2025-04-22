# ATACSeq in Python

## Introduction

We will use PyDeseq2 to run an ATAC-Seq analysis.

You will find three scripts:

- `split_replicas.py`
- `merged_peaks.py`
- `process_deseq.py`

The script `split_replicas.py` will requires two files:

- `encode_ids.txt` and
- `cell_lines.txt`.

The file `cell_lines.txt` contains all the comparisons you want to do. For instance, HSC vs. Erythroblast.
The file `encode_ids.txt` contains all the information from the class. For ATACSeq we will use only the replicas that correspond to ATACSeq.

The first column contains the BAM files and the second the bigBed files.

The script `split_replicas.py` will download both the bam and the bigBed files from the chosen cell lines and will store the metadata in a new file.


The script `merged_peaks.py` will merge the bigBed files (you first have to convert them to bed), and then, for each cell line it will merge the peaks from each replica.
Finally, we will use both merged peaks from both cell lines, and merge them again in a new bed file.

Subsequently, we will use `bedtools` (you need to install it), to find the coverage from all the four bam files (two for each cell line), in the new file containing the 
merged peaks from both cell lines (both replicas).

We will have now created the deseq input file.

Finally, the script `process_deseq.py` will use the deseq input file created in the previous step to run deseq analysis.

There are two modes you can run using `--use_genes` boolean flag:

- `gene` mode, and
- `chromatin` mode.

The `gene` mode will map each triplet `(seqID, start, end)` to the closest gene, while filtering regions that are not within 500bp (adjustable) of a given gene. 
Thus, before running pydeseq2 we will first aggregate the counts for all the genes. Thus, we will do peak annotation.

The `chromatin` mode will run directly pydeseq2 on chromatin counts without peak annotation. This will help us to do differential chromatin accessibility analysis.

I don't know if I did everything correctly! Hope this helps!

Nikol




