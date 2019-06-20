# ALFA: Annotation Landscape For Aligned reads

ALFA provides a global overview of features distribution composing NGS dataset(s). Given a set of aligned reads (BAM files) and an annotation file (GTF format with biotypes), the tool produces plots of the raw and normalized distributions of those reads among genomic categories (stop codon, 5'-UTR, CDS, intergenic, etc.) and biotypes (protein coding genes, miRNA, tRNA, etc.). Whatever the sequencing technique, whatever the organism.

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
![Output example](https://github.com/biocompibens/ALFA/blob/master/Images/outputs_demo.improved.png)

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

    pip install alfa

### Conda
You can install ALFA through the command (might need to add bioconda channel first with "conda config --add channels bioconda"):

    conda install -c biocomp alfa

### Galaxy
You can run ALFA on Galaxy, the instance administrator can find the progam on the mail tool shed.

### Requirements
RAM: 1 GB should be sufficient for a BAM file aligned on the human genome

Dependencies:
  - Bedtools suite (v2.20.0 and above)
  - python libs: argparse, pysam, copy, subprocess, matplotlib, re, progressbar2, multiprocessing

## Quick start
There is a [quick start](https://github.com/biocompibens/ALFA/tree/master/Quick_start) in the respository in order to test the tool installation. To do so, one can go in the directory and run the following command:

    python alfa.py -a quick_start.gtf --chr_len quick_start.chr_len.txt -g quick_start --bam quick_start.bam quick_start
![Quick start terminal](https://github.com/biocompibens/ALFA/blob/master/Images/quick_start.terminal.png)
![Quick start plots](https://github.com/biocompibens/ALFA/blob/master/Images/quick_start.improved.png)

## Detailed example
Here is an illustrated detailed example produced by ALFA from fake input files.

The BAM file contains 10 reads fully mapped to unique genomic categories referring to a GTF file describing a genome made of only one gene (without introns).
The figure shows an illustration of the input BAM file reads distribution on the genome (A). These reads are converted to proportions on the top plot produced by ALFA (B). As an example, 60% of the reads (n=6) are mapped to a CDS region. This plot is then normalized according to the "surface" that each category occupies within a given genome described in the input GTF file (C). Finally, the bottom plot produced by ALFA shows the enrichment/depletion of the different categories (D). For instance, the CDS regions are enriched by a factor of 1.2 since 60% of the nucleotides from the reads map to this feature although the genome is only composed of 50% of CDS regions.
![Detailed example](https://github.com/biocompibens/ALFA/blob/master/Images/detailed_example.png)

## Manual
### ALFA usages
ALFA can be installed from pip or conda, in this case, an executable is provided meaning that the user can call the program with "alfa" directly. Otherwise, if ALFA is installed from a GitHub clone, the user has to type "python alfa.py".

The basic ALFA workflow consists in 2 steps performed at once or separately:

* Generating ALFA genome index files (stranded and unstranded)

> The user supplies an annotations file (GTF format) to generate ALFA indexes that will be used in the 2nd step. The ALFA genome index files are saved and need to be generated only once for each annotation file.

* Intersecting mapped reads with the ALFA genome index files

> The user provides the previously generated ALFA genome indexes as well as mapped reads file(s) (BAM format) that are intersected to determine the proportion of reads (at nucleotide resolution) within each of the genomic features. As an output, ALFA produces a count file for each input as well as plots based on these files.

#### Generating ALFA index files
Usage:

    alfa -a GTF_FILE [-g GENOME_INDEX] [--chr_len CHR_LENGTHS_FILE] [-p NB_CPUS]

Arguments:
* _**-a/--annotation**_ specifies the path to the genomic annotation file (GTF format).
* _**-g/--genome_index**_ defines ALFA index files basename. In absence of this option, the annotation file basename will be used.
* _**--chr_len specifies**_ the path to the tab-delimited text file defining chromosome names and lengths. In absence of this option, lengths will be estimated using the GTF file. Chr_len file example:

> “Chr12    100000”

Important: the GTF file has to be sorted by position. Otherwise, you can use the following command line to sort it:

    sort -k1,1 -k4,4n -k5,5nr GTF_FILE > SORTED_GTF_FILE

#### Processing reads files
Usage:

    alfa -g GENOME_INDEX --bam BAM_FILE1 LABEL1 [BAM_FILE2 LABEL2 …]
                         [-s STRAND] [-d {1,2,3,4}] [-t YMIN YMAX]
                         [-n] [--{pdf, svg, png} output.{pdf, svg, png}]
                         [--keep_ambiguous] [-p NB_CPUS]

Arguments:
* _**-g/--genome_index**_ specifies path and basename of existing ALFA index files
* _**--bam**_ specifies BAM files paths and associated labels (labels are used within output filenames and plots legends)
* _**-s/--strandness**_ specifies the strandness of the library. Authorized values are: ‘unstranded’ (default), ’forward’/’fr-firststrand’ and ‘reverse’/’fr-secondstrand’
* _**-d/--categories_depth**_ specifies the depth for the categories (see [Categories depth](#categories-depth))
* _**-t/--threshold**_ set the y-axis limits
* _**--pdf or --svg or --png**_ specifies the path to save the plots in the chosen format
* _**-n/--no_display**_ do not create and show the plots
* _**--keep_ambiguous**_ specifies that the nucleotides mapping to more than one feature (category or biotype) is equally split between the different features instead of being discarded

Important: BAM files have to be sorted by position. Otherwise, you can use the 'sort' module of samtools suite:

    samtools sort BAM_FILE -T aln.sorted -O bam -o SORTED_BAM_FILE

#### Advanced possibilities
* *Indexing + processing*

> create the ALFA genome indexes and process your BAM files at once using both _**-a/--annotation**_ and _**--bam**_ options. Example:
```
alfa -a quick_start.gtf -g quick_start --bam quick_start.bam quick_start
```

* *Processing bedgraph files*

> provide the data in the coverage bedgraph format to skip the bam to bedgraph and coverageBed steps using _**--bedgraph**_ flag. Example (unstranded and stranded):
```
alfa -g quick_start --bedgraph quick_start.bedgraph quick_start
alfa -g quick_start --bedgraph quick_start.plus.bedgraph quick_start.plus.bedgraph quick_start -s forward
```

* *Running the tool from counts*

> specify the count files previously generated to avoid running the script again on already processed datasets using the _**-c/--counts**_ option (instead of _**--bam**_). Example:
```
alfa -c quick_start.feature_counts.tsv
```

#### Miscellaneous arguments
* _**-p/--processors**_ specifies the number of processors used
* _**-o**_ specifies an output directory where the files are created (default is the current directory, it will be created if it doesn't exist)

### Categories depth
ALFA can assign categories to nucleotides according to different hierarchical levels considered in the GTF file using the _**-d/--categories_depth**_ option.
Here are the features considered in the 4 different levels:

1. Gene / intergenic
2. Exon / intron / intergenic
3. (Default) 5’-UTR / CDS / 3’-UTR / intron / intergenic
4. 5’-UTR / start_codon / CDS / stop_codon / 3’-UTR / intron / intergenic

*Warning: using a non-homogeneous GTF file in term of deep level annotations may lead to inconsistent results due to the fact that, for instance, reads mapping to genes without UTR annotation will increase the CDS category count whereas on the other genes, the UTR categories may be increased.*

### Priorities
By default, as GTF files are built on a hierarchical way, some assumptions are made on categories priorities:
> start_codon/stop_codon > five_prime_utr/three_prime_utr/CDS > exon > transcript > gene

This means, for example, that a nucleotide associated to an *exon* as well as to a *CDS* will be accounted for the *CDS* category.

### Ambiguity
In case of a nucleotide mapping to overlapping categories, by default, it is discarded. Although with the option _**--keep_ambiguous**_, the count will be split equally between the different categories. Overlapping biotype priorities are first solved with the associated category, in case of an equality, the previous principle is applied.

### Unknown feature
If ALFA meets a category that is not referenced in its code, it will be ignored since its priority is tricky to assess. However, an unknown biotype is added on the fly and will be processed.

### Acknowledgments
We thank Felipe Delestro Matos from [IBENS Bioinformatics platform](https://www.ibens.ens.fr/?rubrique55) for the figures and plots design.
