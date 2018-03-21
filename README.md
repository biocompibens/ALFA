# ALFA: Annotation Landscape For Aligned reads

ALFA provides a global overview of features distribution composing NGS dataset(s). Given a set of aligned reads (BAM files) and an annotation file (GTF format), the tool produces plots of the raw and normalized distributions of those reads among genomic categories (stop codon, 5'-UTR, CDS, intergenic, etc.) and biotypes (protein coding genes, miRNA, tRNA, etc.). Whatever the sequencing technique, whatever the organism.

## Table of contents
- [**Outputs**](#outputs)
  - [Plots](#plots)
  - [Count files](#count-files)
- [**Installation**](#installation)
  - [Requirements](#requirements)
- [**Quick start**](#quick-start)
- [**Detailed example**](#detailed-example)
- [**Manual**](#manual)
  - [ALFA usages](#alfa-usages)
  - [Categories depth](#categories-depth)
  - [Priorities](#priorities)
  - [Unknown feature](#unknown-feature)

## Outputs
ALFA produces two outputs: plots and count files.

### Plots
![Output example](https://github.com/biocompibens/ALFA/blob/master/Images/outputs_demo.svg)

**Legend**: Two images display the nucleotides distributions among the different features, one for the categories (A) and another one for the biotypes (B). Each image is composed of two plots: the top ones (1) represents the raw nucleotides fraction of the reads that mapped to each feature in the samples whereas the bottom ones (2) represents the feature enrichment/depletion relative to the genome ”surface” occupied by this feature.

### Count files
For each input sample, a count table with the following information is produced:
 - category/biotype pair
 - count of nucleotides mapping to this feature pair in the sample
 - count of nucleotides belonging to this feature pair in the genome

> Using -c/--count option, this file can be used to plot again the results and avoid the whole program execution.

![Count file example](https://github.com/biocompibens/ALFA/blob/master/Images/counts.head.png)

## Installation
### GitHub
You can either download the code by clicking the ZIP link on this webpage or clone the project using:

    git clone https://github.com/biocompibens/ALFA

### Pip
You can install ALFA through the command:

    pip install ALFA

### Galaxy
You can run ALFA on Galaxy, the instance administrator can find the progam on the mail tool shed.

### Requirements
RAM: 1 GB should be sufficient for a BAM file aligned on the human genome

Dependencies:
  - Bedtools suite (v2.20.0 and above)
  - python libs: argparse, pysam, copy, subprocess, matplotlib, re, progressbar, multiprocessing

## Quick start
There is a [toy dataset](https://github.com/biocompibens/ALFA/tree/master/Toy_dataset) in the respository in order to test the tool installation. To do so, one can go in the directory and run the following command:

    python ../ALFA.py -a toy.gtf -g toy_reference -i toy.bam toy --chr_len chr_len.txt
![Quick start terminal](https://github.com/biocompibens/ALFA/blob/master/Images/quick_start_terminal.png)
![Quick start plots](https://github.com/biocompibens/ALFA/blob/master/Images/quick_start_plots.png)

## Detailed example
Here is an illustrated detailed example produced by ALFA from fake input files.

The BAM file contains 10 reads fully mapped to unique genomic categories referring to a GTF file describing a genome made of only of one gene (without introns).
The figure shows an illustration of the input BAM file reads distribution on the genome. These reads are converted to proportions on the top plot produced by ALFA. As an example, 60% of the reads (n=6) are mapped to a CDS region. This plot is then normalized according to the "surface" that each category occupies within a given genome described in the input GTF file. Finally, the bottom plot produced by ALFA shows the enrichment/depletion of the different categories. For instance, the CDS regions are enriched by a factor of 1.2 since 60% of the nucleotides from the reads map to this feature although the genome is only composed of 50% of CDS regions.
![Detailed example](https://github.com/biocompibens/ALFA/blob/master/Images/detailed.png)

## Manual
### ALFA usages
The basic ALFA workflow consists in 2 steps performed at once or separately:

* Generating genome index files (stranded and unstranded)

> The user supplies an annotations file (in GTF format) to generate indexes that will be used in the 2nd step. The genome index files are saved and need to be generated only once for each annotation file.

* Intersecting mapped reads with the genome index files

> The user provides the previously generated genome indexes as well as mapped reads file(s) (BAM format) that are intersected to determine the proportion of reads (at nucleotide resolution) within each of the genomic features. As an output, ALFA produces a count file for each input as well as plots based on these files.

#### Generating index files
Usage:

    ALFA.py -a GTF_FILE [-g GENOME_INDEX] [--chr_len CHR_LENGTHS_FILE]

Arguments:
* _**-a/--annotation**_ specifies the path to the genomic annotation file (GTF format) to generate indexes.
* _**-g/--genome_index**_ defines index files basename. In absence of this option, the annotation file basename will be used.
* _**-p/--processors**_ specifies the number of processors used
* _**--chr_len specifies**_ the path to the tab-delimited text file defining chromosome names and lengths. In absence of this option, lengths will be estimated using the GTF file. Chr_len file example:

> “Chr12    100000”

Important: the GTF file has to be sorted by position. Otherwise, you can use the following command line to sort it

    sort -k1,1 -k4,4n -k5,5nr GTF_FILE > SORTED_GTF_FILE

#### Processing reads files
Usage:

    ALFA.py -g GENOME_INDEX --bam BAM_FILE1 LABEL1 [BAM_FILE2 LABEL2 …]
                         [-s STRAND] [-n]
                         [-d {1,2,3,4}] [--pdf output.pdf]

Arguments:
* _**-g/--genome_index**_ specifies path and basename of existing index files
* _**--bam**_ specifies BAM files paths and associated labels (labels are used within output filenames and plots legends)
* _**-s/--strandness**_ specifies the strandness of the library. Authorized values are: ‘unstranded’ (default), ’forward’/’fr-firststrand’ and ‘reverse’/’fr-secondstrand’
* _**-d/--categories_depth**_ specifies the depth for the categories (see [Categories depth](#categories-depth))
* _**-p/--processors**_ specifies the number of processors used (will use at most twice the BAM file number)
* _**-t/--threshold**_ set the y-axis limits
* _**--pdf or --svg or --png**_ specifies the path to save the plots in the chosen format
* _**-n/--no_display**_ do not create and show the plots

Important: BAM files have to be sorted by position. Otherwise, you can use the 'sort' module of samtools suite

    samtools sort BAM_FILE -T aln.sorted -O bam -o SORTED_BAM_FILE

#### Advanced possibilities
* *Indexing + processing*

> create the genome indexes and process your BAM files at once using both _**-a/--annotation**_ and _**--bam**_ options.

* *Processing bedgraph files*

> provide the data in the coverage bedgraph format to skip the bam to bedgraph and coverageBed steps using _**--bedgraph**_ flag.

* *Running the tool from counts*

> specify the count files previously generated to avoid running the script again on already processed datasets using the _**-c/--counts**_ option (instead of _**--bam**_).

### Categories depth
ALFA can assign categories to nucleotides according to different hierarchical levels considered in the GTF file using the _**-d/--categories_depth**_ option.
Here are the features considered in the 4 different levels:

1. Gene / intergenic
2. Exon / intron / intergenic
3. 5’-UTR / CDS / 3’-UTR / intron / intergenic (default)
4. 5’-UTR / start_codon / CDS / stop_codon / 3’-UTR / intron / intergenic

*Warning: using a non-homogeneous GTF file in term of deep level annotations may lead to inconsistent results due to the fact that, for instance, reads mapping to genes without UTR annotation will increase the CDS category count whereas on the other genes, the UTR categories may be increased.*

### Priorities
By default, as GTF files are built on a hierarchical way, some assumptions are made on categories priorities.
> start_codon/stop_codon > five_prime_utr/three_prime_utr/CDS > exon > transcript > gene

This means, for example, that a nucleotide found in a *gene* as well as in a *transcript* will be counted within the *transcript* category.
In case of a nucleotide found in two categories of equal priority, the count is split between them.
Overlapping biotype priorities are first solved with the associated category, in case of an equality, the previous principle is applied.

### Unknown feature
If ALFA meets a category that is not reference in its code, it won’t take it into account as its priority is tricky to assess. However, an unknown biotype is added on the fly and will be processed.
