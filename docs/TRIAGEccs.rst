TRIAGEccs
=========

Description
-----------
TRIAGEccs is one of the core functions of the TRIAGE R package v2. It is designed for the unbiased prioritization of regulatory elements - including genes, promoters, enhancers, variants, and other genomic loci of interest - without requiring prior annotation. The function performs genome-wide, single-base resolution scoring of input regions, assigning a Cellular Constraint Score (CCS) to quantify their regulatory importance across a wide range of cell types. The CCS reflects the likelihood that a given genomic region plays a key functional role in regulating cell identity, differentiation, or disease processes.

By enabling CCS calculation for any genomic locus, TRIAGEccs provides a powerful and flexible framework for the discovery and prioritization of regulatory elements, as well as for the interpretation of functional genomics data in complex biological systems. It is particularly useful for analyzing novel genes, long non-coding RNAs (lncRNAs), insertions/deletions (InDels), single-nucleotide polymorphisms (SNPs), and more.

For more details, see: `Sinniah et al., bioRxiv 2024, <https://www.biorxiv.org/content/10.1101/2024.10.28.620690v1.full>`_.


Input and Output
----------------

Input: A BED file containing genomic intervals of interest. These may include any genomic loci of interest.

Output: A tab-delimited file (BED format) with an appended column containing the CCS value for each input interval. Additionally, the function returns a matrix of calculated CCS values for downstream analysis.


Parameters
----------

- `input`: A character string indicating the input BED file.

..

- `rts_genome`: (Optional) A character string specifying the genome-wide repressive tendency score (RTS) BEDGRAPH file. Default is `"TRIAGE_hg19.bedgraph"`. Users working with the hg38 genome build may alternatively use `"TRIAGE_hg38.bedgraph"`.

..

- `output`: (Optional) A character string specifying the name of the output file. Default is `"CCS.bed"`. The output is a BED file with an additional column containing the CCS values.

..

- `exons_only`: (Optional) Logical value indicating whether CCS should be computed only for exon regions. Default is `FALSE`. When `TRUE`, the function uses blockSizes and blockStarts fields (columns 11 and 12 in BED format) to define exons.


Usage Examples
--------------

.. code-block:: R

    # Example 1: Run TRIAGEccs with default RTS file ("TRIAGE_hg19.bedgraph") and output name
    ccs <- TRIAGEccs("input.bed")

    # Example 2: Specify RTS genome file and output file
    ccs <- TRIAGEccs("input.bed", rts_genome = "TRIAGE_hg38.bedgraph", output = "output_file.txt")

    # Example 3: Calculate CCS only for exon regions
    ccs <- TRIAGEccs("input.bed", output = "output_file.txt", exons_only = TRUE)
