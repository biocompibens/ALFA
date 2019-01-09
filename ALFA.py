#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "noel & bahin"
""" ALFA provides a global overview of features distribution composing NGS dataset(s). """

import os
import sys
import re
import numpy as np
import collections
import copy
import argparse
import pysam
import pybedtools
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from multiprocessing import Pool
import progressbar

# To correctly embed the texts when saving plots in svg format
matplotlib.rcParams["svg.fonttype"] = "none"


##########################################################################
#                         FUNCTIONS                                      #
##########################################################################

def init_dict(d, key, init):
    if key not in d:
        d[key] = init


def tryint(s):
    """ Function called by "alphanum_key" function to sort the chromosome names. """
    try:
        return int(s)
    except ValueError:
        return s


def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [tryint(c) for c in re.split("([0-9]+)", s)]


def existing_file(filename):
    """ Checks if filename already exists and exit if so. """
    if os.path.isfile(filename):
        sys.exit("Error: The file '" + filename + "' is about to be produced but already exists in the directory. \n### End of program")


def get_chromosome_names_in_GTF(options):
    """ Function to get the list of chromosome names present in the provided GTF file. """
    chr_list = []
    with open(options.annotation, "r") as GTF_file:
        for line in GTF_file:
            if not line.startswith("#"):
                chr = line.split("\t")[0]
                if chr not in chr_list:
                    chr_list.append(chr)
    return sorted(chr_list)


def get_chromosome_names_in_index(genome_index, index_chrom_list):
    """ Returns the chromosome names in a genome index file. """
    with open(genome_index, "r") as index_fh:
        for line in index_fh:
            if not line.startswith("#") and (line.split("\t")[0] not in index_chrom_list):
                index_chrom_list.append(line.split("\t")[0])
    return index_chrom_list


def GTF_splitter(GTF_file, chunk_basename, size=10000):
    """ Function to split a GTF file into chunks of one chromosome or several chromosomes/scaffolds up to N (default=10k) lines. """
    if os.path.isfile(chunk_basename + "1.gtf"):
        sys.exit("Error: There is already a file called '" + chunk_basename + "1.gtf' in the directory. Running the command would crush this file. Aborting")
    prev_chr = ""  # Chr/scaffold previously processed
    prev_cpt = 0  # Currently building chunk file line counter
    cpt = 0  # Processed chr/scaffold line counter
    cpt_chunk = 1  # Chunk counter
    current_file = open("current.gtf", "w")  # New piece to add to the building chunk file (one chromosome/scaffold)
    # Processing the input GTF file
    with open(GTF_file, "r") as input_file:
        for line in input_file:
            # Burning header lines
            if line.startswith("#"):
                continue
            # Getting the chromosome/scaffold
            chr = line.split("\t")[0]  # Processed chr/scaffold
            # Reaching another chromosome/scaffold
            if (chr != prev_chr) and (prev_chr != ""):
                if cpt > size:
                    # Packing up the processed chr/scaffold
                    current_file.close()
                    os.rename("current.gtf", chunk_basename + str(cpt_chunk) + ".gtf")
                    current_file = open("current.gtf", "w")
                    # Updating counters
                    cpt_chunk += 1
                else:
                    if cpt + prev_cpt > size:
                        # Packing up the currently building chunk file without the last chr/scaffold
                        os.rename("old.gtf", chunk_basename + str(cpt_chunk) + ".gtf")
                        # Updating counters
                        cpt_chunk += 1
                        prev_cpt = 0
                    # Moving the new piece to the currently building chunk file
                    current_file.close()
                    with open("current.gtf", "r") as input_file, open("old.gtf", "a") as output_file:
                        for line in input_file:
                            output_file.write(line)
                    current_file = open("current.gtf", "w")
                    # Updating counters
                    prev_cpt += cpt
                # Updating the processed chr/scaffold line counter and the previous chromosome
                cpt = 0
                prev_chr = chr
                current_file.write(line)
            else:  # First content line or another line for the processed chr/scaffold
                current_file.write(line)
                prev_chr = chr
                cpt += 1
        # Processing the last chromosome/scaffold
        current_file.close()
        if prev_cpt == 0:  # There was only one chromosome/scaffold in the annotation file
            # Packing up the processed chr/scaffold
            os.rename("current.gtf", chunk_basename + str(cpt_chunk) + ".gtf")
        else:
            if prev_cpt + cpt > size:
                # Packing up the processed chr/scaffold
                os.rename("current.gtf", chunk_basename + str(cpt_chunk) + ".gtf")
                cpt_chunk += 1
            else:
                # Moving the new piece to the currently building chunk file
                with open("current.gtf", "r") as input_file, open("old.gtf", "a") as output_file:
                    for line in input_file:
                        output_file.write(line)
                os.remove("current.gtf")
            # Packing up the currently building chunk file without the last chr/scaffold
            os.rename("old.gtf", chunk_basename + str(cpt_chunk) + ".gtf")


def get_chromosome_lengths(options):
    """
    Parse the file containing the chromosome lengths.
    If no length file is provided, browse the annotation file (GTF) to estimate the chromosome sizes.
    """
    lengths = {}
    gtf_chrom_names = set()
    # If the user provided a chromosome length file
    if options.chr_len:
        # Getting the chromosome lengths from the chromosome lengths file
        with open(options.chr_len, "r") as chr_len_fh:
            for line in chr_len_fh:
                try:
                    lengths[line.split("\t")[0]] = int(line.rstrip().split("\t")[1])
                except IndexError:
                    sys.exit("Error: The chromosome lengths file is not correctly formed. It is supposed to be tabulated file with two fields per line.")
        # Getting the chromosome lengths from the GTF file
        with open(options.annotation, "r") as gtf_fh:
            for line in gtf_fh:
                if not line.startswith("#"):
                    gtf_chrom_names.add(line.split("\t")[0])
        # Checking if the chromosomes from the chromosome lengths file are present in the GTF file
        for chrom in lengths:
            if chrom not in gtf_chrom_names:
                print("Warning: chromosome '" + chrom + "' of the chromosome lengths file does not match any chromosome name in the GTF file provided and was ignored.", file=sys.stderr)
        # Checking if the chromosomes from the GTF file are present in the lengths file
        for chrom in gtf_chrom_names:
            if chrom not in lengths:
                print("Warning: at least one chromosome ('" + chrom + "') was found in the GTF file and does not match any chromosome provided in the lengths file.", file=sys.stderr)
                print("\t=> All the chromosome lengths will be approximated using annotations in the GTF file.", file=sys.stderr)  ## MB: better than the first end of line??
                break
        else:
            return lengths
    # If no chromosome lengths file was provided or if at least one chromosome was missing in the file, the end of the last annotation of the chromosome in the GTF file is considered as the chromosome length
    with open(options.annotation, "r") as gtf_fh:
        for line in gtf_fh:
            if not line.startswith("#"):
                chrom = line.split("\t")[0]
                end = int(line.split("\t")[4])
                init_dict(lengths, chrom, 0)
                lengths[chrom] = max(lengths[chrom], end)
    print("The chromosome lengths have been approximated using the GTF file annotations (the stop position of the last annotation of each chromosome is considered as its length).")
    return lengths


def write_index_line(feat, chrom, start, stop, sign, fh):
    """ Write a new line in an index file. """
    # Formatting the features info
    feat_by_biotype = []
    for biot, cat in feat.items():
        feat_by_biotype.append(":".join((str(biot), ",".join(sorted(cat)))))
    # Writing the features info in the index file
    fh.write("\t".join((chrom, start, stop, sign)) + "\t" + "\t".join(feat_by_biotype) + "\n")


def write_index(feat_values, chrom, start, stop, stranded_genome_index, unstranded_genome_index):
    """ Writing the features info in the proper index files. """
    # Writing info to the stranded indexes
    if feat_values[0] != {}:
        write_index_line(feat_values[0], chrom, start, stop, "+", stranded_genome_index)
    else:
        stranded_genome_index.write("\t".join((chrom, start, stop, "+", "antisense\n")))
    if feat_values[1] != {}:
        write_index_line(feat_values[1], chrom, start, stop, "-", stranded_genome_index)
    else:
        stranded_genome_index.write("\t".join((chrom, start, stop, "-", "antisense\n")))
    # Writing info to the unstranded index
    unstranded_feat = dict(feat_values[0], **feat_values[1])
    for name in set(feat_values[0]) & set(feat_values[1]):
        unstranded_feat[name] += feat_values[0][name]
    write_index_line(unstranded_feat, chrom, start, stop, ".", unstranded_genome_index)


def register_interval(features_dict, chrom, stranded_index_fh, unstranded_index_fh):
    """ Write the interval features info into the genome index files. """
    # Writing the interval in the index file
    with open(unstranded_index_fh, "a") as unstranded_index_fh, open(stranded_index_fh, "a") as stranded_index_fh:
        # Initializing the first interval start and features
        sorted_pos = sorted(features_dict["+"].keys())
        interval_start = sorted_pos[0]
        features_plus = features_dict["+"][interval_start]
        features_minus = features_dict["-"][interval_start]
        # Browsing the interval boundaries
        for interval_stop in sorted_pos[1:]:
            # Writing the current interval features to the indexes
            write_index([features_plus, features_minus], chrom, str(interval_start), str(interval_stop), stranded_index_fh, unstranded_index_fh)
            # Initializing the new interval start and features
            interval_start = interval_stop
            # Store current features
            prev_features_plus = features_plus
            prev_features_minus = features_minus
            # Update features
            features_plus = features_dict["+"][interval_start]
            features_minus = features_dict["-"][interval_start]
            # If feature == transcript and prev interval's feature is exon => add intron feature
            for biotype, categ in features_plus.items():
                if categ == ["gene", "transcript"]:
                    if "exon" in prev_features_plus[biotype] or "intron" in prev_features_plus[biotype]:
                        categ.append("intron")
                else:
                    continue
            for biotype, categ in features_minus.items():
                if categ == ["gene", "transcript"]:
                    if "exon" in prev_features_minus[biotype]:
                        categ.append("intron")
                else:
                    continue


def generate_genome_index_1chr(arguments):
    # Getting the arguments
    (annotation, chunk_basename, reverse_strand) = arguments
    # Setting the annotation file basename
    annotation_basename = re.sub(".gtf$", "", annotation)
    # Processing the annotation file
    with open(annotation, "r") as gtf_fh:
        max_value = -1  # Maximum value of the currently processed interval
        intervals_dict = {}
        prev_chrom = ""
        for line in gtf_fh:
            # Processing lines except comment ones
            if not line.startswith("#"):
                # Getting the line info
                line_split = line.rstrip().split("\t")
                chrom = line_split[0]
                cat = line_split[2]
                start = int(line_split[3]) - 1
                stop = int(line_split[4])
                strand = line_split[6]
                antisense = reverse_strand[strand]
                biotype = line_split[8].split("gene_biotype")[1].split(";")[0].strip('" ')
                # Registering stored features info in the genome index file(s) if the new line concerns a new chromosome or the new line concerns an annotation not overlapping previously recorded ones
                if (start > max_value) or (chrom != prev_chrom):
                    # Write the previous features
                    if intervals_dict:
                        register_interval(intervals_dict, prev_chrom, annotation_basename + ".stranded.ALFA_index", annotation_basename + ".unstranded.ALFA_index")
                    if chrom != prev_chrom:
                        with open(chunk_basename + "txt", "a") as input_file:
                            input_file.write(chrom + "\n")
                    prev_chrom = chrom
                    # (Re)Initializing the intervals info dict
                    intervals_dict = {strand: {start: {biotype: [cat]}, stop: {}}, antisense: {start: {}, stop: {}}}
                    max_value = stop

                # Update the dictionary which represents intervals for every distinct annotation
                else:
                    # Storing the intervals on the strand of the current feature
                    stranded_intervals = intervals_dict[strand]
                    added_info = False  # Variable to know if the features info were already added
                    # Browsing the existing boundaries
                    for boundary in sorted(stranded_intervals):
                        # While the GTF line start is after the browsed boundary: keep the browsed boundary features info in case the GTF line start is before the next boundary
                        if boundary < start:
                            stored_feat_strand = copy.deepcopy(dict(stranded_intervals[boundary]))
                            stored_feat_antisense = copy.deepcopy(dict(intervals_dict[antisense][boundary]))
                            continue

                        # The GTF line start is already an existing boundary: store the existing features info (to manage after the GTF line stop) and update it with the GTF line features info
                        elif boundary == start:
                            stored_feat_strand = copy.deepcopy(dict(stranded_intervals[boundary]))
                            stored_feat_antisense = copy.deepcopy(dict(intervals_dict[antisense][boundary]))
                            # Adding the GTF line features info to the interval
                            try:
                                if cat not in stranded_intervals[boundary][biotype]:
                                    stranded_intervals[boundary][biotype].append(cat)
                            except KeyError:  # If the GTF line features info regards a new biotype
                                stranded_intervals[boundary][biotype] = [cat]
                            added_info = True  # The features info were added
                            continue

                        # The browsed boundary is after the GTF line start: add the GTF line features info
                        elif boundary > start:
                            # Create a new boundary for the GTF line start if necessary (if it is between 2 existing boundaries, it was not created before)
                            if not added_info:
                                stranded_intervals[start] = copy.deepcopy(stored_feat_strand)
                                try:
                                    if cat not in stranded_intervals[start][biotype]:
                                        stranded_intervals[start][biotype].append(cat)
                                except KeyError:
                                    stranded_intervals[start][biotype] = [cat]
                                intervals_dict[antisense][start] = copy.deepcopy(stored_feat_antisense)
                                added_info = True  # The features info were added
                            # While the browsed boundary is before the GTF line stop: store the existing features info (to manage after the GTF line stop) and update it with the GTF line features info
                            if boundary < stop:
                                stored_feat_strand = copy.deepcopy(dict(stranded_intervals[boundary]))
                                stored_feat_antisense = copy.deepcopy(dict(intervals_dict[antisense][boundary]))
                                try:
                                    if cat not in stranded_intervals[boundary][biotype]:
                                        stranded_intervals[boundary][biotype].append(cat)
                                except KeyError:
                                    stranded_intervals[boundary][biotype] = [cat]
                            # The GTF line stop is already exists, nothing more to do, the GTF line features info are integrated
                            elif boundary == stop:
                                break
                            # The browsed boundary is after the GTF line stop: create a new boundary for the GTF line stop (with the stored features info)
                            else:  # boundary > stop
                                stranded_intervals[stop] = copy.deepcopy(stored_feat_strand)
                                intervals_dict[antisense][stop] = copy.deepcopy(stored_feat_antisense)
                                break  # The GTF line features info are integrated
                    # If the GTF line stop is after the last boundary, extend the dictionary
                    if stop > max_value:
                        max_value = stop
                        stranded_intervals[stop] = {}
                        intervals_dict[antisense][stop] = {}

        # Store the categories of the last chromosome
        register_interval(intervals_dict, chrom, annotation_basename + ".stranded.ALFA_index", annotation_basename + ".unstranded.ALFA_index")
        if chrom != prev_chrom:
            with open(chunk_basename + "txt", "a") as input_file:
                input_file.write(chrom + "\n")
    return None


def generate_genome_index(unstranded_genome_index, stranded_genome_index, chrom_sizes, chunk_basename, options, reverse_strand):
    """ Create an index of the genome annotations and save it in a file. """
    # Write the chromosome lengths as comment lines before the genome index
    with open(unstranded_genome_index, "w") as unstranded_index_fh, open(stranded_genome_index, "w") as stranded_index_fh:
        for key, value in list(chrom_sizes.items()):
            unstranded_index_fh.write("#%s\t%s\n" % (key, value))
            stranded_index_fh.write("#%s\t%s\n" % (key, value))
    # Chunk file list creation
    chunks = np.array([f for f in os.listdir(".") if f.startswith(chunk_basename) and f.endswith(".gtf")])
    # Sorting the chunks by file size
    file_sizes = np.array([os.stat(f).st_size for f in chunks])
    chunks = chunks[file_sizes.argsort()]
    # Progress bar to track the genome indexes creation
    pbar = progressbar.ProgressBar(widgets=["Indexing the genome ", progressbar.Percentage(), " ", progressbar.Bar(), progressbar.Timer()], maxval=len(chunks)).start()
    pool = Pool(options.nb_processors)
    list(pbar(pool.imap_unordered(generate_genome_index_1chr, zip(chunks, [chunk_basename] * len(chunks), [reverse_strand] * len(chunks)))))
    """
    # Non-parallel version for debugging
    for f in chunks:
        generate_genome_index_1chr(f)
    """


def merge_index_chunks(unstranded_genome_index, stranded_genome_index, chunk_basename):
    """ Merges the genome index chunks into a single file. """
    for fh, strandness in zip([unstranded_genome_index, stranded_genome_index], ["unstranded", "stranded"]):
        files = [f for f in os.listdir(".") if f.startswith(chunk_basename) and f.endswith("." + strandness + ".ALFA_index")]
        with open(fh, "a") as output_file:
            for file in sorted(files):
                with open(file, "r") as input_file:
                    for line in input_file:
                        output_file.write(line)


def chunks_cleaner(chunk_basename):
    """ Cleans the chunks created to index the genome. """
    for f in os.listdir("."):
        if f.startswith(chunk_basename):
            os.remove(f)


def count_genome_features(cpt, features, start, stop, biotype_prios, options, prios, unknown_cat, coverage=1):
    """ Reads genome index and registers feature counts. """
    # If no biotype priority: category with the highest priority for each found biotype has the same weight (1/n_biotypes)
    if not biotype_prios:
        nb_biot = len(features)
        if options.keep_ambiguous and nb_biot != 1:
            # Increment "ambiguous" counter if more than 1 biotype
            try:
                cpt[("ambiguous", "ambiguous")] += (int(stop) - int(start)) * coverage
            except:
                cpt[("ambiguous", "ambiguous")] = (int(stop) - int(start)) * coverage
        else:
            # For each categ(s)/biotype pairs
            for feat in features:
                cur_prio = 0
                # Separate categorie(s) and biotype
                try:
                    biot, cats = feat.split(":")
                # Error if the feature is "antisense": update the "antisense/antisense" counts
                except ValueError:
                    try:
                        cpt[("opposite_strand", "opposite_strand")] += (int(stop) - int(start)) * coverage / float(nb_biot)
                    except KeyError:
                        cpt[("opposite_strand", "opposite_strand")] = (int(stop) - int(start)) * coverage / float(nb_biot)
                    return None
                # Browse the categories and get only the one(s) with highest priority
                for cat in cats.split(","):
                    try:
                        prio = prios[cat]
                    except KeyError:
                        if cat not in unknown_cat:
                            print("Warning: Unknown categorie '%s' found and ignored.\n" % cat, file=sys.stderr)
                        unknown_cat.add(cat)
                        continue
                    # Check if the category has a highest priority than the current one
                    if prio > cur_prio:
                        cur_prio = prio
                        cur_cat = {cat}
                    if prio == cur_prio:
                        cur_cat.add(cat)

                nb_cat = len(cur_cat)
                if options.keep_ambiguous and nb_cat != 1:
                    # Increment "ambiguous" counter if more than 1 category
                    try:
                        cpt[("ambiguous", "ambiguous")] += (int(stop) - int(start)) * coverage / (float(nb_biot))
                    except KeyError:
                        cpt[("ambiguous", "ambiguous")] = (int(stop) - int(start)) * coverage / (float(nb_biot))
                else:
                    # Increase each counts by the coverage divided by the number of categories and biotypes
                    for cat in cur_cat:
                        try:
                            cpt[(cat, biot)] += (int(stop) - int(start)) * coverage / (float(nb_biot * nb_cat))
                        except KeyError:
                            cpt[(cat, biot)] = (int(stop) - int(start)) * coverage / (float(nb_biot * nb_cat))
    else:
        # TODO: Add an option to provide biotype priorities and handle it
        pass


def read_index(genome_index, lengths, cpt_genome, options, prios, index_chrom_list, biotype_prios, unknown_cat):
    """ Parse index files info (chromosomes list, lengths and genome features). """
    with open(genome_index, "r") as index_fh:
        for line in index_fh:
            if line.startswith("#"):
                lengths[line.split("\t")[0][1:]] = int(line.split("\t")[1])
            else:
                chrom = line.split("\t")[0]
                if chrom not in index_chrom_list:
                    index_chrom_list.append(chrom)
                count_genome_features(cpt_genome, line.rstrip().split("\t")[4:], line.split("\t")[1], line.split("\t")[2], biotype_prios, options, prios, unknown_cat)


def run_genomecov(arguments):
    """ Run genomecov (from Bedtools through pybedtools lib) for a set of parameters to produce a BedGraph file. """
    (strand, bam_file, sample_label, name, bedgraph_extension) = arguments
    # Load the BAM file
    input_file = pybedtools.BedTool(bam_file)
    # Run genomecov
    if strand == "":
        input_file.genome_coverage(bg=True, split=True).saveas(sample_label + name + bedgraph_extension)
    else:
        input_file.genome_coverage(bg=True, split=True, strand=strand).saveas(sample_label + name + bedgraph_extension)
    return None


def generate_bedgraph_files_parallel(sample_labels, bam_files, options, bedgraph_extension):
    """ Creates, through multi-processors, BedGraph files from BAM ones. """
    # Sorting the BAM file on size to process the biggest first
    files = list(zip(sample_labels, bam_files, [os.stat(i).st_size for i in bam_files]))
    files.sort(key=lambda p: p[2], reverse=True)
    # Defining parameters sets to provide to the genomecov instances to run
    parameter_sets = []
    for l, b, s in files:
        # If the dataset is stranded, one BedGraph file for each strand is created
        if options.strandness in ["forward", "fr-firststrand"]:
            parameter_sets.append(["+", b, l, ".plus", bedgraph_extension])
            parameter_sets.append(["-", b, l, ".minus", bedgraph_extension])
        elif options.strandness in ["reverse", "fr-secondstrand", bedgraph_extension]:
            parameter_sets.append(["-", b, l, ".plus", bedgraph_extension])
            parameter_sets.append(["+", b, l, ".minus", bedgraph_extension])
        else:
            parameter_sets.append(["", b, l, "", bedgraph_extension])
    # Setting the progressbar
    pbar = progressbar.ProgressBar(widgets=["Generating the BedGraph files ", progressbar.Percentage(), progressbar.Bar(), progressbar.SimpleProgress(), "|", progressbar.Timer()], maxval=len(parameter_sets)).start()
    # Setting the processors number
    pool = Pool(options.nb_processors)
    # Running the instances
    list(pbar(pool.imap_unordered(run_genomecov, parameter_sets)))
    pybedtools.cleanup()  # If pybedtools can't finish properly but the program is not stopped, the /tmp will be cleaned anyway from pybedtools temp files
    return None


def read_gtf(gtf_index_file, sign):
    global gtf_line, gtf_chrom, gtf_start, gtf_stop, gtf_features, endGTF
    strand = ""
    while strand != sign:
        gtf_line = gtf_index_file.readline()
        # If the GTF file is finished
        if not gtf_line:
            endGTF = True
            return endGTF
        splitline = gtf_line.rstrip().split("\t")
        try:
            strand = splitline[3]
        # strand information can not be found in the file header
        except IndexError:
            pass
    gtf_chrom = splitline[0]
    gtf_features = splitline[4:]
    gtf_start, gtf_stop = list(map(int, splitline[1:3]))
    return endGTF


def intersect_bedgraphs_and_index_to_count_categories_1_file(arguments):
    global gtf_line, gtf_chrom, gtf_start, gtf_stop, gtf_cat, endGTF
    (sample_labels, bedgraph_files, biotype_prios, bedgraph_extension, genome_index, options, prios, index_chrom_list, unknown_cat, strand, sign) = arguments  # For Python 3
    unknown_chrom = []
    cpt = {}  # Counter for the nucleotides in the BAM input file(s)
    prev_chrom = ""
    endGTF = False  # Reaching the next chr or the end of the GTF index
    intergenic_adds = 0.0
    # Checking if the BedGraph file is empty
    if os.stat(bedgraph_files + strand + bedgraph_extension).st_size == 0:
        return sample_labels, sign, {}, []
    with open(bedgraph_files + strand + bedgraph_extension, "r") as bedgraph_fh:
        # Running through the BedGraph file
        for bam_line in bedgraph_fh:
            # Getting the BAM line info
            bam_chrom = bam_line.split("\t")[0]
            bam_start, bam_stop, bam_cpt = list(map(float, bam_line.split("\t")[1:4]))
            # Skip the line if the chromosome is not in the index
            if bam_chrom not in index_chrom_list:
                if bam_chrom not in unknown_chrom:
                    unknown_chrom.append(bam_chrom)
                continue
            # If this is a new chromosome (or the first one)
            if bam_chrom != prev_chrom:
                intergenic_adds = 0.0
                # Closing the GTF file if it was open (exception caught only for the first chr)
                try:
                    gtf_index_file.close()
                except UnboundLocalError:
                    pass
                # (Re)opening the GTF index and looking for the first line of the matching chr
                gtf_index_file = open(genome_index, "r")
                endGTF = False
                read_gtf(gtf_index_file, sign)
                while bam_chrom != gtf_chrom:
                    read_gtf(gtf_index_file, sign)
                    if endGTF:
                        break
                prev_chrom = bam_chrom

            # Looking for the first matching annotation in the GTF index
            while (not endGTF) and (gtf_chrom == bam_chrom) and (bam_start >= gtf_stop):
                read_gtf(gtf_index_file, sign)
                if gtf_chrom != bam_chrom:
                    endGTF = True
            # Processing BAM lines before the first GTF annotation if there are
            if bam_start < gtf_start:
                # Increase the "intergenic" category counter with all or part of the BAM interval
                try:
                    intergenic_adds += min(bam_stop, gtf_start) - bam_start
                    cpt[("intergenic", "intergenic")] += (min(bam_stop, gtf_start) - bam_start) * bam_cpt
                except KeyError:
                    cpt[("intergenic", "intergenic")] = (min(bam_stop, gtf_start) - bam_start) * bam_cpt
                # Go to next line if the BAM interval was totally upstream the first GTF annotation, carry on with the remaining part otherwise
                if endGTF or (bam_stop <= gtf_start):
                    continue
                else:
                    bam_start = gtf_start

            # We can start the crossover
            while not endGTF:
                # Update category counter
                count_genome_features(cpt, gtf_features, bam_start, min(bam_stop, gtf_stop), biotype_prios, options, prios, unknown_cat, coverage=bam_cpt)
                # Read the next GTF file line if the BAM line is not entirely covered
                if bam_stop > gtf_stop:
                    # Update the BAM start pointer
                    bam_start = gtf_stop
                    endGTF = read_gtf(gtf_index_file, sign)
                    # If we read a new chromosome in the GTF file then it is considered finished
                    if bam_chrom != gtf_chrom:
                        endGTF = True
                    if endGTF:
                        break
                else:
                    # Next if the BAM line is entirely covered
                    bam_start = bam_stop
                    break

            # Processing the end of the BAM line if necessary
            if endGTF and (bam_stop > bam_start):
                try:
                    cpt[("intergenic", "intergenic")] += (bam_stop - bam_start) * bam_cpt
                except KeyError:
                    cpt[("intergenic", "intergenic")] = (bam_stop - bam_start) * bam_cpt
        # In stranded mode, if one of the BedGraph files doesn't have any of the chromosomes from the reference file, the error is not detected during the preprocessing (a chromosome is found within the other BedGraph file)
        try:
            gtf_index_file.close()
        except UnboundLocalError:
            # Then the file is not opened and can't be closed
            pass
    return sample_labels, sign, cpt, unknown_chrom


def intersect_bedgraphs_and_index_to_count_categories(sample_labels, bedgraph_files, options, bedgraph_extension, genome_index, prios, index_chrom_list, unknown_cat, biotype_prios=None): ## MB: To review
    # Initializing variables
    unknown_chrom = []
    cpt = {}  # Counter for the nucleotides in the BAM input file(s)
    if bedgraph_files == []:
        bedgraph_files = sample_labels
    if options.strandness == "unstranded":
        strands = [("", ".")]
    else:
        strands = [(".plus", "+"), (".minus", "-")]

    # Initializing the progress bar
    pbar = progressbar.ProgressBar(widgets=["Intersecting BAM and genome ", progressbar.Percentage(), " ", progressbar.Bar(), progressbar.SimpleProgress(), "|", progressbar.Timer()], maxval=len(sample_labels) * len(strands)).start()
    pool = Pool(options.nb_processors)
    inputs = [sample + strand for sample in zip(sample_labels, bedgraph_files, [biotype_prios] * len(sample_labels), [bedgraph_extension] * len(sample_labels), [genome_index] * len(sample_labels), [options] * len(sample_labels), [prios] * len(sample_labels), [index_chrom_list] * len(sample_labels), [unknown_cat] * len(sample_labels)) for strand in strands]

    # Running the intersection in parallel
    results = list(pbar(pool.imap_unordered(intersect_bedgraphs_and_index_to_count_categories_1_file, inputs)))
    """
    # Non-parallel version for debugging
    results = []
    for i in inputs:
        results.append(intersect_bedgraphs_and_index_to_count_categories_1_file(i))
    """

    # Reformatting the results
    for res in results:
        init_dict(cpt, res[0], {})
        if res[1] != {}:
            init_dict(cpt[res[0]], res[1], res[2])
        unknown_chrom.append(res[3])
    # Merging strands counts for the same samples
    final_cpt = {}
    for sample in cpt:
        final_cpt[sample] = {}
        for strand in strands:
            for feat in cpt[sample][strand[1]]:
                try:
                    final_cpt[sample][feat] += cpt[sample][strand[1]][feat]
                except KeyError:
                    final_cpt[sample][feat] = cpt[sample][strand[1]][feat]

    print("Unknown chromosomes: " + str(set([i for u in unknown_chrom for i in u])) + ".")
    return final_cpt


def write_counts_in_files(cpt, genome_counts):
    """ Writes the biotype/category counts in an output file. """
    for sample_label, counters in list(cpt.items()):
        sample_label = "_".join(re.findall(r"[\w\-']+", sample_label))
        with open(sample_label + ".ALFA_feature_counts.tsv", "w") as output_fh:
            output_fh.write("#Category,biotype\tCounts_in_BAM/BedGraph\tSize_in_genome\n")
            for features_pair, counts in list(counters.items()):
                output_fh.write("%s\t%s\t%s\n" % (",".join(features_pair), counts, genome_counts[features_pair]))


def read_counts(sample_labels, counts_files):
    """ Reads the counts from an input file. """
    cpt = {}
    cpt_genome = {}
    for sample_label, filename in zip(sample_labels, counts_files):
        cpt[sample_label] = {}
        with open(filename, "r") as counts_fh:
            for line in counts_fh:
                if not line.startswith("#"):
                    feature = tuple(line.split("\t")[0].split(","))
                    cpt[sample_label][feature] = float(line.split("\t")[1])
                    cpt_genome[feature] = float(line.rstrip().split("\t")[2])
    return cpt, cpt_genome


def group_counts_by_categ(cpt, cpt_genome, final, selected_biotype):
    final_cat_cpt = {}
    final_genome_cpt = {}
    filtered_cat_cpt = {}
    for f in cpt:
        final_cat_cpt[f] = {}
        filtered_cat_cpt[f] = {}
        for final_cat in final:
            tot = 0
            tot_filter = 0
            tot_genome = 0
            for cat in final[final_cat]:
                for key, value in list(cpt[f].items()):
                    if key[0] == cat:
                        tot += value
                        tot_genome += cpt_genome[key]
                        if key[1] == selected_biotype:
                            tot_filter += value
            # output_file.write('\t'.join((final_cat, str(tot))) + '\n')
            # print '\t'.join((final_cat, str(tot)))
            final_cat_cpt[f][final_cat] = tot
            if tot_genome == 0:
                final_genome_cpt[final_cat] = 1e-100
            else:
                final_genome_cpt[final_cat] = tot_genome
            filtered_cat_cpt[f][final_cat] = tot_filter
    # if "antisense" in final_genome_cpt: final_genome_cpt["antisense"] = 0
    return final_cat_cpt, final_genome_cpt, filtered_cat_cpt


def group_counts_by_biotype(cpt, cpt_genome, biotypes):
    final_cpt = {}
    final_genome_cpt = {}
    for f in cpt:
        final_cpt[f] = {}
        for biot in biotypes:
            tot = 0
            tot_genome = 0
            try:
                for final_biot in biotypes[biot]:
                    for key, value in list(cpt[f].items()):
                        if key[1] == final_biot:
                            tot += value
                            # if key[1] != 'antisense':
                            tot_genome += cpt_genome[key]
            except:
                for key, value in list(cpt[f].items()):
                    if key[1] == biot:
                        tot += value
                        tot_genome += cpt_genome[key]
            if tot != 0:
                final_cpt[f][biot] = tot
                final_genome_cpt[biot] = tot_genome
    return final_cpt, final_genome_cpt


def display_percentage_of_ambiguous(cpt, count_files_option=False):
    if count_files_option:
        print("INFO: Ambiguous counts were discarded in at least one sample\n" \
+              "      (see --ambiguous option for more information)")
    else:
        print("INFO: Reads matching ambiguous annotation have been discarded.\n" \
              + "      To change this option, please see \"--ambiguous\" help.")
    print("INFO: Percentage of ambiguous counts:")

    # Loop for each sample
    for sample, counts in cpt.items():
        # Compute and display the percentage of ambiguous counts
        try:
            ambiguous = counts[('ambiguous', 'ambiguous')]
            total = sum([count for feat, count in counts.items() if 'intergenic' not in feat])
            print("    {!s:25.25} {:5.2f}% of ambiguous".format(sample, float(ambiguous / total) * 100))
        # If ambiguous is not in the count file
        except KeyError:
            if count_files_option:
                print(
                    "    {!s:25.25} {:3.2f}% of ambiguous (this sample may have been processed with --ambiguous option)".format(sample, 0))
            else:
                print("    {!s:25.25} {:3.2f}% of ambiguous".format(sample, 0))


def recategorize_the_counts(cpt, cpt_genome, final):
    final_cat_cpt = {}
    final_genome_cpt = {}
    for f in cpt:
        # print "\nFinal categories for",f,"sample"
        final_cat_cpt[f] = {}
        for final_cat in final:
            tot = 0
            tot_genome = 0
            for cat in final[final_cat]:
                tot += cpt[f][cat]
                tot_genome += cpt_genome[cat]
            # output_file.write('\t'.join((final_cat, str(tot))) + '\n')
            # print '\t'.join((final_cat, str(tot)))
            final_cat_cpt[f][final_cat] = tot
            final_genome_cpt[final_cat] = tot_genome
    return final_cat_cpt, final_genome_cpt


def one_sample_plot(ordered_categs, percentages, enrichment, n_cat, index, index_enrichment, bar_width, counts_type,
                    title, sample_labels):
    ### Initialization
    fig = plt.figure(figsize=(13, 9))
    ax1 = plt.subplot2grid((2, 4), (0, 0), colspan=2)
    ax2 = plt.subplot2grid((2, 4), (1, 0), colspan=2)
    cmap = plt.get_cmap("Spectral")
    cols = [cmap(x) for x in range(0, 256, 256 / n_cat)]
    if title:
        ax1.set_title(title + "in: %s" % sample_labels[0])
    else:
        ax1.set_title(counts_type + " distribution in mapped reads in: %s" % sample_labels[0])
    ax2.set_title("Normalized counts of " + counts_type)

    ### Barplots
    # First barplot: percentage of reads in each categorie
    ax1.bar(index, percentages, bar_width,
            color=cols)
    # Second barplot: enrichment relative to the genome for each categ
    # (the reads count in a categ is divided by the categ size in the genome)
    ax2.bar(index_enrichment, enrichment, bar_width,
            color=cols, )
    ### Piecharts
    pielabels = [ordered_categs[i] if percentages[i] > 0.025 else "" for i in range(n_cat)]
    sum_enrichment = np.sum(enrichment)
    pielabels_enrichment = [ordered_categs[i] if enrichment[i] / sum_enrichment > 0.025 else "" for i in range(n_cat)]
    # Categories piechart
    ax3 = plt.subplot2grid((2, 4), (0, 2))
    pie_wedge_collection, texts = ax3.pie(percentages, labels=pielabels, shadow=True, colors=cols)
    # Enrichment piechart
    ax4 = plt.subplot2grid((2, 4), (1, 2))
    pie_wedge_collection, texts = ax4.pie(enrichment, labels=pielabels_enrichment, shadow=True, colors=cols)
    # Add legends (append percentages to labels)
    labels = [" ".join((ordered_categs[i], "({:.1%})".format(percentages[i]))) for i in range(len(ordered_categs))]
    ax3.legend(pie_wedge_collection, labels, loc="center", fancybox=True, shadow=True, prop={"size": "medium"},
               bbox_to_anchor=(1.7, 0.5))
    labels = [" ".join((ordered_categs[i], "({:.1%})".format(enrichment[i] / sum_enrichment))) for i in
              range(len(ordered_categs))]  # if ordered_categs[i] != "antisense"]
    ax4.legend(pie_wedge_collection, labels, loc="center", fancybox=True, shadow=True, prop={"size": "medium"},
               bbox_to_anchor=(1.7, 0.5))
    # Set aspect ratio to be equal so that pie is drawn as a circle
    ax3.set_aspect("equal")
    ax4.set_aspect("equal")
    return fig, ax1, ax2


def make_plot(sample_labels, ordered_categs, categ_counts, genome_counts, counts_type, options, title=None, categ_groups=[]):

    #Test matplotlib version. If __version__ >= 2, use a shift value to correct the positions of bars and xticks
    if int(matplotlib.__version__[0]) == 2:
        shift_mpl = 0.5
    else:
        shift_mpl = 0
    # From ordered_categs, keep only the features (categs or biotypes) that we can find in at least one sample.
    existing_categs = set()
    for sample in list(categ_counts.values()):
        existing_categs |= set(sample.keys())
    ordered_categs = list(filter(existing_categs.__contains__, ordered_categs))
    xlabels = [cat if len(cat.split("_")) == 1 else "\n".join(cat.split("_")) if cat.split("_")[0] != 'undescribed' else "\n".join(["und.",cat.split("_")[1]]) for cat in ordered_categs]
    n_cat = len(ordered_categs)
    #n_exp = len(sample_names)
    nb_samples = len(categ_counts)
    ##Initialization of the matrix of counts (nrow=nb_experiements, ncol=nb_categorie)
    #counts = np.matrix(np.zeros(shape=(n_exp, n_cat)))
    counts = np.matrix(np.zeros(shape=(nb_samples, n_cat)))
    '''
    for exp in xrange(len(sample_names)):
        for cat in xrange(len(ordered_categs)):
            try:
                counts[exp, cat] = categ_counts[sample_names[exp]][ordered_categs[cat]]
            except:
                pass
    '''
    for sample_label in sample_labels:
        for cat in range(len(ordered_categs)):
            try:
                counts[sample_labels.index(sample_label), cat] = categ_counts[sample_label][ordered_categs[cat]]
            except:
                pass

    ##Normalize the categorie sizes by the total size to get percentages
    sizes = []
    sizes_sum = 0
    for cat in ordered_categs:
        sizes.append(genome_counts[cat])
        sizes_sum += genome_counts[cat]
    if "opposite_strand" in ordered_categs:
        antisense_pos = ordered_categs.index("opposite_strand")
        sizes[antisense_pos] = 1e-100
    for cpt in range(len(sizes)):
        sizes[cpt] /= float(sizes_sum)

    ## Create array which contains the percentage of reads in each categ for every sample
    percentages = np.array(counts / np.sum(counts, axis=1))
    ## Create the enrichment array (counts divided by the categorie sizes in the genome)
    enrichment = np.array(percentages / sizes)
    if "antisense_pos" in locals():
        '''
        for i in xrange(len(sample_names)):
            enrichment[i][antisense_pos] = 0
        '''
        for n in range(nb_samples):
            enrichment[n][antisense_pos] = 0

    # enrichment=np.log(np.array(percentages/sizes))
    #for exp in xrange(n_exp):
    for n in range(nb_samples):
        for i in range(n_cat):
            val = enrichment[n][i]
            if val > 1:
                enrichment[n][i] = val - 1
            elif val == 1 or val == 0:
                enrichment[n][i] = 0
            else:
                enrichment[n][i] = -1 / val + 1

    #### Finally, produce the plot

    ##Get the colors from the colormap
    ncolor = 16
    cmap = ["#e47878", "#68b4e5", "#a3ea9b", "#ea9cf3", "#e5c957", "#a3ecd1", "#e97ca0", "#66d985", "#8e7ae5",
            "#b3e04b", "#b884e4", "#e4e758", "#738ee3", "#e76688", "#70dddd", "#e49261"]
    '''
    if n_exp > ncolor:
        cmap = plt.get_cmap("Set3", n_exp)
        cmap = [cmap(i) for i in xrange(n_exp)]
    '''
    if nb_samples > ncolor:
        cmap = plt.get_cmap("tab20", nb_samples)
        cmap = [cmap(i) for i in range(nb_samples)]

    ## Parameters for the plot
    opacity = 1
    # Create a vector which contains the position of each bar
    index = np.arange(n_cat)
    # Size of the bars (depends on the categs number)
    #bar_width = 0.9 / n_exp
    bar_width = 0.9 / nb_samples

    ##Initialise the subplot
    # if there is only one sample, also plot piecharts
    # if n_exp == 1 and counts_type.lower() == 'categories':
    # fig, ax1, ax2 = one_sample_plot(ordered_categs, percentages[0], enrichment[0], n_cat, index, bar_width, counts_type, title)
    ## If more than one sample
    # else:
    if counts_type.lower() != "categories":
        #fig, (ax1, ax2) = plt.subplots(2, figsize=(5 + (n_cat + 2 * n_exp) / 3, 10))
        fig, (ax1, ax2) = plt.subplots(2, figsize=(5 + (n_cat + 2 * nb_samples) / 3, 10))
    else:
        #fig, (ax1, ax2) = plt.subplots(2, figsize=(5 + (n_cat + 2 * n_exp) / 3, 10))
        fig, (ax1, ax2) = plt.subplots(2, figsize=(5 + (n_cat + 2 * nb_samples) / 3, 10))
    # Store the bars objects for percentages plot
    rects = []
    # Store the bars objects for enrichment plot
    rects_enrichment = []
    # For each sample/experiment
    #for i in range(n_exp):
    for sample_label in sample_labels:
        # First barplot: percentage of reads in each categorie
        n = sample_labels.index(sample_label)
        #ax1.bar(index + i * bar_width, percentages[i], bar_width,
        rects.append(ax1.bar(index + n * bar_width + shift_mpl/nb_samples, percentages[n], bar_width,
                alpha=opacity,
                #color=cmap[i],
                color=cmap[n],
                #label=sample_names[i], edgecolor="#FFFFFF", lw=0)
                label=sample_label, edgecolor="#FFFFFF", lw=0))
        # Second barplot: enrichment relative to the genome for each categ
        # (the reads count in a categ is divided by the categ size in the genome)
        rects_enrichment.append(ax2.bar(index + n * bar_width + shift_mpl/nb_samples, enrichment[n], bar_width,
                             alpha=opacity,
                             #color=cmap[i],
                             color=cmap[n],
                             #label=sample_names[i], edgecolor=cmap[i], lw=0))
                             label=sample_label, edgecolor=cmap[n], lw=0))

    ## Graphical options for the plot
    # Adding of the legend
    #if n_exp < 10:
    if nb_samples < 10:
        ax1.legend(loc="best", frameon=False)
        legend_ncol = 1
    #elif n_exp < 19:
    elif nb_samples < 19:
        legend_ncol = 2
    else:
        legend_ncol = 3
    ax1.legend(loc="best", frameon=False, ncol=legend_ncol)
    ax2.legend(loc="best", frameon=False, ncol=legend_ncol)
    # ax2.legend(loc='upper center',bbox_to_anchor=(0.5,-0.1), fancybox=True, shadow=True)
    # Main titles
    if title:
        ax1.set_title(title)
    else:
        ax1.set_title(counts_type + " counts")
    ax2.set_title(counts_type + " normalized counts")

    # Adding enrichment baseline
    # ax2.axhline(y=0,color='black',linestyle='dashed',linewidth='1.5')
    # Axes limits
    ax1.set_xlim(-0.1, len(ordered_categs) + 0.1)
    if len(sizes) == 1: ax1.set_xlim(-2, 3)
    ax2.set_xlim(ax1.get_xlim())
    # Set axis limits (max/min values + 5% margin)
    ax2_ymin = []
    ax2_ymax = []
    for sample_values in enrichment:
        ax2_ymin.append(min(sample_values))
        ax2_ymax.append(max(sample_values))
    ax2_ymax = max(ax2_ymax)
    ax2_ymin = min(ax2_ymin)
    margin_top, margin_bottom = (abs(0.05 * ax2_ymax), abs(0.05 * ax2_ymin))
    ax1.set_ylim(0, ax1.get_ylim()[1] * 1.05)
    if options.threshold:
        threshold_bottom = -abs(float(options.threshold[0])) + 1
        threshold_top = abs(float(options.threshold[1]) - 1)

        #for i in xrange(n_exp):
        for n in range(nb_samples):
            for y in range(n_cat):
                #val = enrichment[i][y]
                val = enrichment[n][y]
                if not np.isnan(val) and not (threshold_bottom < val < threshold_top):
                    #rect = rects_enrichment[i][y]
                    rect = rects_enrichment[n][y]
                    rect_height = rect.get_height()
                    if rect.get_y() < 0:
                        diff = rect_height + threshold_bottom
                        rect.set_y(threshold_bottom)
                        ax2_ymin = threshold_bottom
                        margin_bottom = 0
                    else:
                        diff = rect_height - threshold_top
                        ax2_ymax = threshold_top
                        margin_top = 0
                    rect.set_height(rect.get_height() - diff)
    if margin_top != 0 and margin_bottom != 0:
        margin_top, margin_bottom = [max(margin_top, margin_bottom) for i in range(2)]
    ax2.set_ylim(ax2_ymin - margin_bottom, ax2_ymax + margin_top)
    # Y axis title
    ax1.set_ylabel("Proportion of reads (%)")
    ax2.set_ylabel("Enrichment relative to genome")

    # Add the categories on the x-axis
    #ax1.set_xticks(index + bar_width * n_exp / 2)
    ax1.set_xticks(index + bar_width * nb_samples / 2)
    #ax2.set_xticks(index + bar_width * n_exp / 2)
    ax2.set_xticks(index + bar_width * nb_samples / 2)
    if counts_type.lower() != "categories":
        ax1.set_xticklabels(ordered_categs, rotation="30", ha="right")
        ax2.set_xticklabels(ordered_categs, rotation="30", ha="right")
    else:
        ax1.set_xticklabels(xlabels)
        ax2.set_xticklabels(xlabels)

    # Display fractions values in percentages
    ax1.set_yticklabels([str(int(i * 100)) for i in ax1.get_yticks()])
    # Correct y-axis ticks labels for enrichment subplot
    # ax2.set_yticklabels([str(i+1)+"$^{+1}$" if i>0 else 1 if i==0 else str(-(i-1))+"$^{-1}$" for i in ax2.get_yticks()])
    yticks = list(ax2.get_yticks())
    yticks = [yticks[i] - 1 if yticks[i] > 9 else yticks[i] + 1 if yticks[i] < -9 else yticks[i] for i in range(len(yticks))]
    ax2.set_yticks(yticks)
    ax2.set_yticklabels([str(int(i + 1)) if i > 0 and i % 1 == 0 else str(
        i + 1) if i > 0 else 1 if i == 0 else str(
        int(-(i - 1))) + "$^{-1}$" if i < 0 and i % 1 == 0 else str(-(i - 1)) + "$^{-1}$" for i in ax2.get_yticks()])
    # ax2.set_yticklabels([i+1 if i>0 else 1 if i==0 else "$\\frac{1}{%s}$" %-(i-1) for i in ax2.get_yticks()])
    # Change appearance of 'antisense' bars on enrichment plot since we cannot calculate an enrichment for this artificial category
    if "antisense_pos" in locals():  # ax2.text(antisense_pos+bar_width/2,ax2.get_ylim()[1]/10,'NA')
        #for i in xrange(n_exp):
        for n in range(nb_samples):
            #rect = rects_enrichment[i][antisense_pos]
            rect = rects_enrichment[n][antisense_pos]
            rect.set_y(ax2.get_ylim()[0])
            rect.set_height(ax2.get_ylim()[1] - ax2.get_ylim()[0])
            rect.set_hatch("/")
            rect.set_fill(False)
            rect.set_linewidth(0)
            # rect.set_color('lightgrey')
            # rect.set_edgecolor('#EDEDED')
            rect.set_color("#EDEDED")
        #ax2.text(index[antisense_pos] + bar_width * n_exp / 2 - 0.1, (ax2_ymax + ax2_ymin) / 2, "NA")
        ax2.text(index[antisense_pos] + bar_width * nb_samples / 2 - 0.1, (ax2_ymax + ax2_ymin) / 2, "NA")

    # Add text for features absent in sample & correct for bars too small to be seen
    for n in range(nb_samples):
        for y in range(n_cat):
            # if no counts in sample for this feature, display "Abs. in sample"
            if percentages[n][y] == 0:
                txt = ax1.text(y + bar_width * (n + 0.5), 0.02, "Abs.", rotation="vertical", color=cmap[n],
                               horizontalalignment="center", verticalalignment="bottom")
                txt.set_path_effects([PathEffects.Stroke(linewidth=0.5), PathEffects.Normal()])
            else:
                # if percentage value is lower than 1% of the plot height, modify it to fit this minimum value
                if rects[n][y].get_height() < 5e-3 * (ax1.get_ylim()[1] - ax1.get_ylim()[0]):
                    rects[n][y].set_height(5e-3 * (ax1.get_ylim()[1] - ax1.get_ylim()[0]))
                # if enrichment value equal to 0, increase the line width to see the bar on the plot
                if enrichment[n][y] == 0:
                    #rects_enrichment[i][y].set_linewidth(1)
                    rects_enrichment[n][y].set_linewidth(1)
                # if enrichment value is too small to be seen, increase the bar height to 1% of the plot height
                elif abs(rects_enrichment[n][y].get_height()) < 1e-2 * (ax2_ymax - ax2_ymin):
                    # Correction for negative value in Matplotlib v1
                    if rects_enrichment[n][y].get_y() < 0:
                        rects_enrichment[n][y].set_height(1e-2 * (ax2_ymax - ax2_ymin))
                        rects_enrichment[n][y].set_y(-1e-2 * (ax2_ymax - ax2_ymin))
                    # Correction for negative value in Matplotlib v2
                    elif rects_enrichment[n][y].get_height() < 0:
                        rects_enrichment[n][y].set_height(-1e-2 * (ax2_ymax - ax2_ymin))
                    # Correction for positive value
                    else:
                        rects_enrichment[n][y].set_height(1e-2 * (ax2_ymax - ax2_ymin))
    # Remove top/right/bottom axes
    for ax in [ax1, ax2]:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        #ax.tick_params(axis="x", which="both", bottom="on", top="off", labelbottom="on")
        ax.tick_params(axis="x", which="both", bottom=True, top=False, labelbottom=True)
        #ax.tick_params(axis="y", which="both", left="on", right="off", labelleft="on")
        ax.tick_params(axis="y", which="both", left=True, right=False, labelleft=True)

    ### Add second axis with categ groups
    annotate_group(categ_groups, label=None, ax=ax1)
    annotate_group(categ_groups, label=None, ax=ax2)

    ### Adjust figure margins to
    if counts_type.lower() == "categories":
        plt.tight_layout(h_pad=5.0)
        fig.subplots_adjust(bottom=0.1)
    else:
        plt.tight_layout()

    ## Displaying or saving the plot
    if not options.pdf and not options.svg and not options.png:
        plt.show()
    else:  # If any of the 3 plot output format is set
        for output_basename, output_format in [(options.pdf, "pdf"), (options.svg, "svg"), (options.png, "png")]:
            if output_basename:
                # Checking if the file extension have been specified and removing it if so
                if output_basename.endswith("." + output_format):
                    output_basename = output_basename[:-4]
                # Saving the plot
                plt.savefig(".".join((output_basename.rstrip("." + output_format), counts_type, output_format)))
    plt.close()


def annotate_group(groups, ax=None, label=None, labeloffset=30):
    """Annotates the categories with their parent group and add x-axis label"""

    def annotate(ax, name, left, right, y, pad):
        """Draw the group annotation"""
        arrow = ax.annotate(name, xy=(left, y), xycoords="data",
                            xytext=(right, y - pad), textcoords="data",
                            annotation_clip=False, verticalalignment="top",
                            horizontalalignment="center", linespacing=2.0,
                            arrowprops={'arrowstyle': "-", 'shrinkA': 0, 'shrinkB': 0,
                                        'connectionstyle': "angle,angleB=90,angleA=0,rad=5"}
                            )
        return arrow

    if ax is None:
        ax = plt.gca()
    level = 0
    for level in range(len(groups)):
        grp = groups[level]
        for name, coord in list(grp.items()):
            ymin = ax.get_ylim()[0] - np.ptp(ax.get_ylim()) * 0.12 - np.ptp(ax.get_ylim()) * 0.05 * (level)
            ypad = 0.01 * np.ptp(ax.get_ylim())
            xcenter = np.mean(coord)
            annotate(ax, name, coord[0], xcenter, ymin, ypad)
            annotate(ax, name, coord[1], xcenter, ymin, ypad)

    if label is not None:
        # Define xlabel and position it according to the number of group levels
        ax.annotate(label,
                    xy=(0.5, 0), xycoords="axes fraction",
                    xytext=(0, -labeloffset - (level + 1) * 15), textcoords="offset points",
                    verticalalignment="top", horizontalalignment="center")

    return


def filter_categs_on_biotype(selected_biotype, cpt):
    filtered_cpt = {}
    for sample in cpt:
        filtered_cpt[sample] = {}
        for feature, count in list(cpt[sample].items()):
            if feature[1] == selected_biotype:
                filtered_cpt[sample][feature[0]] = count
    return filtered_cpt

def usage_message():
    return """
    README on GitHub: https://github.com/biocompibens/ALFA/blob/master/README.md
    
    Generate ALFA genome indexes:
         python ALFA.py -a GTF [-g GENOME_INDEX_BASENAME]
                                        [--chr_len CHR_LENGTHS_FILE]
                                        [-p NB_PROCESSORS]
    Process BAM file(s):
         python ALFA.py -g GENOME_INDEX_BASENAME --bam BAM1 LABEL1 [BAM2 LABEL2 ...]
                                        [--bedgraph] [-s STRAND]
                                        [-d {1,2,3,4}] [-t YMIN YMAX]
                                        [--keep_ambiguous]
                                        [-n] [--pdf output.pdf] [--svg output.svg] [--png output.png]
                                        [-p NB_PROCESSORS]
    Index genome + process BAM files(s):
         python ALFA.py -a GTF [-g GENOME_INDEX_BASENAME] [--chr_len CHR_LENGTHS_FILE]
                                        --bam BAM1 LABEL1 [BAM2 LABEL2 ...]
                                        [--bedgraph][-s STRAND]
                                        [-d {1,2,3,4}] [-t YMIN YMAX]
                                        [--keep_ambiguous]
                                        [-n] [--pdf output.pdf] [--svg output.svg] [--png output.png]
                                        [-p NB_PROCESSORS]

    Process previously created ALFA count file(s):
         python ALFA.py -c COUNTS1 [COUNTS2 ...]
                                        [-s STRAND]
                                        [-d {1,2,3,4}] [-t YMIN YMAX]
                                        [-n] [--pdf output.pdf] [--svg output.svg] [--png output.png]

        """


##########################################################################
#                           MAIN                                         #
##########################################################################


def main():

    #### Parse command line arguments and store them in the variable options
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, usage=usage_message())
    parser.add_argument("--version", action="version", version="version 0.25",
                        help="Show ALFA version number and exit.\n\n-----------\n\n")
    # Options regarding the index
    parser.add_argument("-g", "--genome_index", metavar="GENOME_INDEX_BASENAME",
                        help="Genome index files path and basename for existing index, or path and basename for new index creation.\n\n")
    parser.add_argument("-a", "--annotation", metavar="GTF_FILE", help="Genomic annotations file (GTF format).\n\n")
    parser.add_argument("--chr_len", help="Tabulated file containing chromosome names and lengths.\n\n-----------\n\n")

    # Options regarding the intersection step
    parser.add_argument("--bam", metavar=("BAM1 LABEL1", ""), nargs="+",
                        help="Input BAM file(s) and label(s). The BAM files must be sorted by position.\n\n")
    parser.add_argument("--bedgraph", metavar=("BEDGRAPH1 LABEL1", ""), nargs="+", help="Use this options if your input(s) is/are BedGraph file(s). If stranded, provide the BedGraph files\nfor each strand for all samples (e.g. '--bedgraph file.plus.bedgraph file.minus.bedgraph LABEL').\n\n")
    parser.add_argument("-c", "--counts", metavar=("COUNTS1", ""), nargs="+",
                        help="Use this options instead of '--bam/--bedgraph' to provide ALFA counts files as input \ninstead of bam/bedgraph files.\n\n")
    parser.add_argument("-s", "--strandness", default="unstranded",
                        choices=["unstranded", "forward", "reverse", "fr-firststrand", "fr-secondstrand"], metavar="",
                        help="Library orientation. Choose within: 'unstranded', "
                             "'forward'/'fr-firststrand' \nor 'reverse'/'fr-secondstrand'. "
                             "(Default: 'unstranded')\n\n-----------\n\n")

    # Options regarding the plot
    #parser.add_argument("--biotype_filter", help=argparse.SUPPRESS)  # "Make an extra plot of categories distribution using only counts of the specified biotype."
    parser.add_argument("-d", "--categories_depth", type=int, default=3, choices=list(range(1, 5)),
                        help="Use this option to set the hierarchical level that will be considered in the GTF file (default=3): \n(1) gene,intergenic; \n(2) intron,exon,intergenic; \n(3) 5'UTR,CDS,3'UTR,intron,intergenic; \n(4) start_codon,5'UTR,CDS,3'UTR,stop_codon,intron,intergenic. \n\n")
    parser.add_argument("--pdf", nargs="?", const="ALFA_plots.pdf",
                        help="Save produced plots in PDF format at the specified path ('categories_plots.pdf' if no argument provided).\n\n")
    parser.add_argument("--png", nargs="?", const="ALFA_plots.png",
                        help="Save produced plots in PNG format with the provided argument as basename \n('categories.png' and 'biotypes.png' if no argument provided).\n\n")
    parser.add_argument("--svg", nargs="?", const="ALFA_plots.svg",
                        help="Save produced plots in SVG format with the provided argument as basename \nor 'categories.svg' and 'biotypes.svg' if no argument provided.\n\n")
    parser.add_argument("-n", "--no_display", action="store_const", const=True, default=False, help="Do not display plots.\n\n") # We have to add "const=None" to avoid a bug in argparse
    parser.add_argument("-t", "--threshold", dest="threshold", nargs=2, metavar=("YMIN", "YMAX"), type=float,
                        help="Set ordinate axis limits for enrichment plots.\n\n")
    parser.add_argument("-p", "--processors", dest="nb_processors", type=int, default=1, help="Set the number of processors used for multi-processing operations.\n\n")
    parser.add_argument("--keep_ambiguous", action="store_const", const=False, default=True, help="Keep reads mapping to different features (discarded by default).\n\n")
    parser.add_argument("--temp_dir", dest="temp_dir", default="/tmp/", help="Temp directory to store pybedtools files ('/tmp/' by default).\n\n")

    if len(sys.argv) == 1:
        parser.print_usage()
        sys.exit(1)

    options = parser.parse_args()
    print("### ALFA ###")

    # Sample labels and file paths
    labels = []
    bams = []
    bedgraphs = []
    bedgraph_extension = ".bedgraph"
    count_files = []
    lengths = {}
    index_chrom_list = []  # Not a set because we need to sort it according to the chromosome names later

    # Booleans for steps to be executed
    generate_index = False
    generate_BedGraph = False
    intersect_indexes_BedGraph = False
    generate_plot = False

    # Checking whether the script will be able to write in the current directory
    if not os.access(".", os.W_OK):
        # The only exception would be if the user already has the count files and wants to display the plots without saving them
        if not options.counts or options.pdf or options.svg or options.png:
            sys.exit("Error: write permission denied in the directory, ALFA will not be able to run correctly.\n### End of program")

    #### Check arguments conformity and define which steps have to be performed
    print("### Checking parameters")
    # Checking whether there is at least one parameter among "-a", "--bam/bedgraph" and "-c"
    if not (options.counts or options.bam or options.bedgraph or options.annotation):
        sys.exit("Error: argument(s) are missing. At least '-a', '--bam', '--bedgraph' or '-c' is required. Please refer to help (-h/--help) and usage cases for more details.\n### End of program")
    # Checking the parameters related to the counts file(s)
    if options.counts:
        # Checking whether the counts file(s) exist
        for filename in options.counts:
            if not os.path.isfile(filename):
                sys.exit("Error: the file '" + filename + "' doesn't exist.\n### End of program")
            # Checking whether the counts file(s) have the correct "ALFA_feature_counts.tsv" extension
            if not filename.endswith(".ALFA_feature_counts.tsv"):
                sys.exit("Error: the counts file '" + filename + "' doesn't have the extension 'ALFA_feature_counts.tsv'.\n### End of program")
            # Registering the sample labels
            label = os.path.basename(filename)
            label = re.sub(".ALFA_feature_counts.tsv", "", label)
            label = "_".join(re.findall(r"[\w\-']+", label))
            labels.append(label)
            # Registering the counts filename
            count_files.append(filename)
    else:
        # If the counts are not provided, then at least '-a' or '-g' arguments are mandatory
        if not options.annotation and not options.genome_index:
            sys.exit("Error: at least '-a' or '-g' argument is missing.\n###End of program")
        # Declare genome_index variables (either from the parameter "-g" of from the annotation filename)
        if options.genome_index:
            genome_index_basename = options.genome_index
        else:
            # Otherwise the GTF filename without extension will be used as the basename
            genome_index_basename = options.annotation.split("/")[-1].split(".gtf")[0]
        # Setting the stranded and unstranded ALFA index files and choosing the one to use
        stranded_genome_index = genome_index_basename + ".stranded.ALFA_index"
        unstranded_genome_index = genome_index_basename + ".unstranded.ALFA_index"
        if options.strandness == "unstranded":
            genome_index = unstranded_genome_index
        else:
            genome_index = stranded_genome_index
        # Checking the parameters related to genome annotation and ALFA index files
        if options.annotation:
            # Checking whether the annotation file exists
            if not os.path.isfile(options.annotation):
                sys.exit("Error: the file '" + options.annotation + "' doesn't exist.\n### End of program")
            # Checking whether the chromosomes lengths file exists
            if options.chr_len and not os.path.isfile(options.chr_len):
                sys.exit("Error: the file '" + options.chr_len + "' doesn't exist.\n### End of program")
            # Checking if the annotation file is correctly formatted (a biotype associated to each line)
            with open (options.annotation, 'r') as input_file:
                for line in input_file:
                    if not line.startswith("#"):
                        try:
                            biot_split = line.split("gene_biotype")[1]
                        except IndexError:
                            sys.exit("Error: at least one feature in the annotation file doesn't have a biotype description. ALFA won't be able to work robustly.\n=>" + line.rstrip() + "\n### End of program")
            # If "--genome_index" parameter is present, setting the future ALFA index basename
            if options.genome_index:
                genome_index_basename = options.genome_index
            else:
                # Otherwise the GTF filename without extension will be used as the basename
                genome_index_basename = options.annotation.split("/")[-1].split(".gtf")[0]
            # Set this step as a task to process
            generate_index = True
            # Check if the ALFA indexes already exist and warning if BAM/BedGrap(s) is/are provided, raise an error otherwise
            if os.path.isfile(genome_index_basename + ".stranded.ALFA_index") or os.path.isfile(
                            genome_index_basename + ".unstranded.ALFA_index"):
                if options.bam or options.bedgraph:
                    print("Warning: an ALFA index file named '%s' already exists and will be used. If you want to create a new index, please delete this file or specify another path." % (genome_index_basename + ".(un)stranded.ALFA_index"), file=sys.stderr)
                    generate_index = False
                else:
                    sys.exit("Error: an ALFA index file named '%s' already exists. If you want to create a new index, please delete this file or specify an other path.\n### End of program" % (
                            genome_index_basename + ".(un)stranded.ALFA_index"))
        elif options.genome_index:  # Checking whether the ALFA index files exist if no annotation file was provided
            if not os.path.isfile(options.genome_index + ".stranded.ALFA_index") or not os.path.isfile(
                            options.genome_index + ".unstranded.ALFA_index"):
                sys.exit("Error: the file '" + options.genome_index + ".stranded.ALFA_index" + "' and/or the file '" + options.genome_index + ".unstranded.ALFA_index" + "' doesn't exist.\n### End of program")
            genome_index_basename = options.genome_index

        # Getting the reference genome chromosome list to check whether there is at least one in common with each BAM/BedGraph file
        if options.bam or options.bedgraph:
            if options.annotation:
                # Checking the chromosomes list from GTF file
                reference_chr_list = get_chromosome_names_in_GTF(options)
            else:
                # Checking chromosome list from genome index
                reference_chr_list = get_chromosome_names_in_index(genome_index, index_chrom_list)

        # Checking parameters related to the BAM file(s)
        if options.bam:
            # Checking the input BAM file(s)
            if len(options.bam) % 2 != 0:
                sys.exit("Error: Make sure to follow the expected format: --bam BAM_file1 Label1 [BAM_file2 Label2 ...].\n### End of program ###")
            for sample_package_nb in range(0, len(options.bam), 2):
                # Check whether the BAM file exists
                if not os.path.isfile(options.bam[sample_package_nb]):
                    sys.exit("Error: the file '" + options.bam[sample_package_nb] + "' doesn't exist.\n### End of program")
                # Check whether the BAM file has a correct extension
                if not options.bam[sample_package_nb].endswith(".bam"):
                    sys.exit("Error: at least one of the BAM file(s) doesn't have a '.bam' extension.\n"
                             "Make sure to follow the expected format: --bam BAM_file1 Label1 [BAM_file2 Label2 ...].\n### End of program ###")
                # Registering the BAM filename
                bams.append(options.bam[sample_package_nb])
                # Registering the label(s) (all that is not a character or a "minus" will be transformer into a "_")
                label = "_".join(re.findall(r"[\w\-']+", options.bam[sample_package_nb + 1]))
                labels.append(label)
                # Checking whether the BedGraph file(s) don't already exist
                if options.strandness == "unstranded":
                    existing_file(label + bedgraph_extension)
                else:
                    existing_file(label + ".plus" + bedgraph_extension)
                    existing_file(label + ".minus" + bedgraph_extension)
                # Checking whether the counts file(s) that will be created already exist
                if os.path.isfile(label + ".ALFA_feature_counts.tsv"):
                    sys.exit("Error: The file '" + label + ".ALFA_feature_counts.tsv' is about to be produced but already exists in the directory. \n### End of program")
                # Listing the BAM chromosome(s) to check whether there is at least one common with the reference genome
                BAM_chr_list = pysam.AlignmentFile(options.bam[sample_package_nb], "r").references
                # Checking if there is at least one common chromosome name between the reference genome and the processed BAM file
                if not any(i in reference_chr_list for i in BAM_chr_list):
                    print(("Reference genome chromosome(s): " + str(reference_chr_list)))
                    print(("BAM file chromosome(s): " + str(list(BAM_chr_list))))
                    sys.exit("Error: no matching chromosome between the BAM file '" + options.bam[sample_package_nb] + "' and the reference genome.\n### End of program")
            # Set these steps as a tasks to process
            generate_BedGraph = True
            intersect_indexes_BedGraph = True

        # Checking parameters related to the BedGraph file(s)
        if options.bedgraph:
            # Determining the number of files (BedGraph + label(s)) expected
            if options.strandness == "unstranded":
                sample_file_nb = 2
            else:
                sample_file_nb = 3
            # Setting and checking the BedGraph extension
            bedgraph_extension = "." + options.bedgraph[0].split(".")[-1]
            # Checking the input BedGraph file(s)
            for sample_package_nb in range(0, len(options.bedgraph), sample_file_nb):
                for sample_file in range(0, sample_file_nb - 1):
                    # Check whether the BedGraph file exists
                    if not os.path.isfile(options.bedgraph[sample_package_nb + sample_file]):
                        sys.exit("Error: the file '" + options.bedgraph[
                            sample_package_nb + sample_file] + "' doesn't exist.\n### End of program")
                    # Check whether the BedGraph file has a correct extension
                    if options.bedgraph[sample_package_nb + sample_file].split(".")[-1] not in ("bedgraph", "bg"):
                        sys.exit("Error: at least one of the BedGraph files doesn't have a '.bedgraph'/'.bg' extension."
                                 "Make sure to follow the expected format: --bedgraph BedGraph_file1 Label1 [BedGraph_file2 Label2 ...].\n"
                                 "Or for stranded samples: --bedgraph BedGraph_file1_plus BedGraph_file2_minus Label1 [...].\n### End of program ###")
                    # Listing the BedGraph chromosome(s) to check whether there is at least one common with the reference genome
                    bedgraph_chr_list = []
                    with open(options.bedgraph[sample_package_nb + sample_file], "r") as bedgraph_file:
                        for line in bedgraph_file:
                            chr = line.split("\t")[0]
                            if chr not in bedgraph_chr_list:
                                bedgraph_chr_list.append(chr)
                    # Checking if there is at least one common chromosome name between the reference genome and the processed BedGraph file
                    if not any(i in reference_chr_list for i in bedgraph_chr_list):
                        print(("Reference genome chromosome(s): " + str(reference_chr_list)))
                        print(("BedGraph file chromosome(s): " + str(list(bedgraph_chr_list))))
                        sys.exit("Error: no matching chromosome between the BedGraph file '" + options.bedgraph[
                            sample_package_nb + sample_file] + "' and the reference genome.\n### End of program")
                    # Register the BedGraph filename(s)
                    bedgraphs.append(re.sub("(.(plus|minus))?" + bedgraph_extension, "", options.bedgraph[sample_package_nb + sample_file]))
                # Registering the label(s) (all that is not a character or a "minus" will be transformer into a "_")
                label = "_".join(re.findall(r"[\w\-']+", options.bedgraph[sample_package_nb + sample_file_nb - 1]))
                labels.append(label)
                # Checking whether the count file(s) that will be created already exist
                existing_file(label + ".ALFA_feature_counts.tsv")
            # Set this step as a task to process
            intersect_indexes_BedGraph = True

    # Checking that there is no duplicated labels
    if len(labels) != len(set(labels)):
        sys.exit("Error: at least one label is duplicated.\n### End of program")

    # Setting whether plots should be generated
    try:
        # Checking if a X server environment variable is set
        x_server = os.environ['DISPLAY']
        # Plots are generated except if the flag "no_display" was used or if there is only the "-a" and eventually the "-g" argument
        if not options.no_display and (options.counts or options.svg or options.pdf or options.png or options.bam or options.bedgraph):
            generate_plot = True
            # Checking the sample number for the colors
            if len(labels) > 20:
                print("Warning: there are more than 20 samples, some colors on the plot will be duplicated.", file=sys.stderr)
    except KeyError:
        print("Warning: your current configuration does not allow graphical interface ('$DISPLAY' variable is not set). Plotting step will not be performed.", file=sys.stderr)

    # Setting the temp directory if specified
    if options.temp_dir:
        pybedtools.set_tempdir(options.temp_dir)


    #### Initialization of some variables

    # Miscellaneous variables
    reverse_strand = {"+": "-", "-": "+"}
    samples = collections.OrderedDict()  # Structure: {<label>: [<filename1>(, <filename2>)]} # Example: {'Toy': ['toy.bam']}

    # Initializing the category priority order, coding biotypes and the final list
    prios = {"start_codon": 4, "stop_codon": 4, "five_prime_utr": 3, "three_prime_utr": 3, "UTR": 3, "CDS": 3,
             "exon": 2, "intron": 2, "transcript": 1.5, "gene": 1, "opposite_strand": 0, "intergenic": -1}

    biotype_prios = None
    # biotype_prios = {"protein_coding":1, "miRNA":2}

    categs_level1 = {"gene": ["five_prime_utr", "three_prime_utr", "UTR", "CDS", "exon", "intron", "start_codon",
                             "stop_codon", "transcript", "gene"],
                    "intergenic": ["intergenic"],
                    "opposite_strand": ["opposite_strand"]}

    categs_level2 = {"exons": ["five_prime_utr", "three_prime_utr", "UTR", "CDS", "exon", "start_codon", "stop_codon"],
                    "introns": ["intron"],
                    "undescribed_genes": ["transcript", "gene"],
                    "intergenic": ["intergenic"],
                    "opposite_strand": ["opposite_strand"]}

    categs_level3 = {"5UTR": ["five_prime_utr", "UTR"],
                    "CDS": ["CDS", "start_codon", "stop_codon"],
                    "3UTR": ["three_prime_utr"],
                    "undescribed_exons": ["exon"],
                    "introns": ["intron"],
                    "undescribed_genes": ["transcript", "gene"],
                    "intergenic": ["intergenic"],
                    "opposite_strand": ["opposite_strand"]}

    categs_level4 = {"5UTR": ["five_prime_utr", "UTR"],
                     "start": ["start_codon"],
                     "stop": ["stop_codon"],
                     "CDS_body": ["CDS"],
                     "3UTR": ["three_prime_utr"],
                     "undescribed_exons": ["exon"],
                     "introns": ["intron"],
                     "undescribed_genes": ["transcript", "gene"],
                     "intergenic": ["intergenic"],
                     "opposite_strand": ["opposite_strand"]}

    # categs_groups = [categs_group4, categs_group3, categs_group2, categs_group1]  # Order and merging for the final plot
    categs_levels = [categs_level1, categs_level2, categs_level3, categs_level4]

    parent_categ_level1 = []
    parent_categ_level2 = [{"gene": [0.5, 2.5]}]
    parent_categ_level3 = [{"exon": [0.5, 3.5]}, {"gene": [0.5, 5.5]}]
    parent_categ_level4 = [{"CDS": [1.5, 3.5]}, {"exon": [0.5, 5.5]}, {"gene": [0.5, 7.5]}]
    parent_categ_groups = [parent_categ_level1, parent_categ_level2, parent_categ_level3, parent_categ_level4]

    cat_list = ["5UTR", "start", "CDS", "CDS_body", "stop", "3UTR", "exons", "undescribed_exons", "introns", "gene", "undescribed_genes", "intergenic", "opposite_strand", "ambiguous"]

    # biotypes list
    biotypes = {"protein_coding", "polymorphic_pseudogene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene", "IG_C_gene",
                "IG_D_gene", "IG_J_gene", "IG_V_gene", "3prime_overlapping_ncrna", "lincRNA", "macro_lncRNA", "miRNA",
                "misc_RNA", "Mt_rRNA", "Mt_tRNA", "processed_transcript", "ribozyme", "rRNA", "scaRNA", "sense_intronic",
                "sense_overlapping", "snoRNA", "snRNA", "sRNA", "TEC", "vaultRNA", "opposite_strand",
                "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene",
                "translated_unprocessed_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene", "unitary_pseudogene",
                "unprocessed_pseudogene", "processed_pseudogene", "IG_C_pseudogene", "IG_J_pseudogene", "IG_V_pseudogene",
                "pseudogene", "ncRNA", "tRNA"}  # Type: set (to access quickly)

    # Grouping of biotypes:
    biotypes_group1 = {"protein_coding": ["protein_coding"],
                       "pseudogenes": ["polymorphic_pseudogene", "transcribed_processed_pseudogene",
                                       "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene",
                                       "translated_unprocessed_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene",
                                       "unitary_pseudogene", "unprocessed_pseudogene", "processed_pseudogene",
                                       "IG_C_pseudogene", "IG_J_pseudogene", "IG_V_pseudogene", "pseudogene"],
                       "TR": ["TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene"],
                       "IG": ["IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene"],
                       "MT_RNA": ["Mt_rRNA", "Mt_tRNA"],
                       "ncRNA": ["lincRNA", "macro_lncRNA", "3prime_overlapping_ncrna", "ncRNA"],
                       "others": ["misc_RNA", "processed_transcript", "ribozyme", "scaRNA", "sense_intronic",
                                  "sense_overlapping", "TEC", "vaultRNA"],
                       "opposite_strand": ["opposite_strand"]}
    for biot in ["miRNA", "snoRNA", "snRNA", "rRNA", "sRNA", "tRNA"]:
        biotypes_group1[biot] = [biot]

    # Initializing the unknown features list
    unknown_cat = set()
    unknown_biot = set()

    # Initializing the genome category counter dict
    cpt_genome = {}


    ## Executing the step(s)

    # Indexes generation
    if generate_index:
        # Running the index generation commands
        print("# Generating the genome index files")
        # Getting the PID of the process as a unique random number
        pid = os.getpgrp()
        chunk_basename = "chunk.ALFA." + str(pid) + "."
        # Splitting the GTF file into chunks
        GTF_splitter(options.annotation, chunk_basename)
        # Getting chromosomes lengths
        lengths = get_chromosome_lengths(options)
        # Generating the index files
        generate_genome_index(unstranded_genome_index, stranded_genome_index, lengths, chunk_basename, options, reverse_strand)
        # Merging the genome index chunks
        merge_index_chunks(unstranded_genome_index, stranded_genome_index, chunk_basename)
        # Displaying the list of indexed chromosomes
        for f in os.listdir("."):
            if f.startswith(chunk_basename) and f.endswith(".txt"):
                with open(f, "r") as input_file:
                    for line in input_file:
                        index_chrom_list.append(line.rstrip())
        index_chrom_list.sort(key=alphanum_key)
        print("Indexed chromosomes: " + ", ".join(index_chrom_list))
        chunks_cleaner(chunk_basename)
    if generate_index or not options.counts:
        # Getting index info
        read_index(genome_index, lengths, cpt_genome, options, prios, index_chrom_list, biotype_prios, unknown_cat)
        # Computing the genome intergenic count: sum of the chr lengths minus sum of the genome annotated intervals
        cpt_genome[("intergenic", "intergenic")] = sum(lengths.values()) - sum([v for x, v in cpt_genome.items() if x != ("opposite_strand", "opposite_strand")])

    # BedGraph files generation
    if generate_BedGraph:
        print("# Generating the BedGraph files")
        pybedtools.set_tempdir(options.temp_dir)
        generate_bedgraph_files_parallel(labels, bams, options, bedgraph_extension)

    # Indexes and BedGraph files intersection
    if intersect_indexes_BedGraph:
        print("# Intersecting index and BedGraph files")
        cpt = intersect_bedgraphs_and_index_to_count_categories(labels, bedgraphs, options, bedgraph_extension, genome_index, prios, index_chrom_list, unknown_cat)  # TODO: Write the counts to an output file
        write_counts_in_files(cpt, cpt_genome)

    ## Plot generation ## MB: all the section still to review
    if generate_plot:
        print("# Generating plots")
        # If input files are the categories counts, the first step is to load them
        if options.counts:
            #cpt, cpt_genome, sample_names = read_counts(options.counts)
            cpt, cpt_genome = read_counts(labels, count_files)
        # Managing the unknown biotypes
        for sample_label, counters in list(cpt.items()):
            for (cat, biot) in counters:
                if biot not in biotypes:
                    biotypes.add(biot)
                    biotypes_group1["others"].append(biot)
        biotypes = sorted(biotypes)
        # Moving antisense cat to the end of the list
        biotypes.remove("opposite_strand")
        biotypes.append("opposite_strand")
        # Do not plot ambiguous on biotypes plot
        try:
            biotypes.remove("ambiguous")
        except ValueError:
            pass
        biotypes_group1 = sorted(biotypes_group1)
        # Filtering biotypes if necessary
        filtered_biotype = None
        """
        if options.biotype_filter:
            for sample_label in cpt:
                for feature in cpt[sample_label]:
                    biotype = feature[1]
                    if options.biotype_filter.lower() == biotype.lower():
                        selected_biotype = biotype
                        break
            if filtered_biotype:
                print "\nWarning: biotype '" + options.biotype_filter + "' not found. Please check the biotype name and that this biotype exists in your sample(s)."
        """

        ## Generate the categories plot
        # Recategorizing within the final categories and plot generation
        final_cats = categs_levels[options.categories_depth - 1]
        parent_categs = parent_categ_groups[options.categories_depth - 1]
        final_cat_cpt, final_genome_cpt, filtered_cat_cpt = group_counts_by_categ(cpt, cpt_genome, final_cats, filtered_biotype)
        # If ambiguous features were discarded, print the percentage for each sample
        if options.keep_ambiguous and not options.counts:
            display_percentage_of_ambiguous(cpt)
        # If only counts are provided, check whether 'ambiguous' feature exists in at least one sample and then display the percentages
        elif options.counts and any([('ambiguous', 'ambiguous') in features for features in list(cpt.values())]):
            display_percentage_of_ambiguous(cpt, options.counts)
        # Remove the "opposite_strand" category if the library type is "unstranded" ## MB: if options.strandness == "unstranded": cat_list.remove("opposite_strand")??
        for dic in cpt.values():
            if ("opposite_strand", "opposite_strand") in dic.keys(): break
        else:
            cat_list.remove("opposite_strand")
        make_plot(labels, cat_list, final_cat_cpt, final_genome_cpt, "Categories", options, categ_groups=parent_categs)
        if filtered_biotype:
            make_plot(labels, cat_list, filtered_cat_cpt, final_genome_cpt, "Categories", options, title="Categories distribution for '" + filtered_biotype + "' biotype", categ_groups= parent_categs)
        ## Generate the biotypes plot
        # Recategorization within the final biotypes and plot generation
        final_cat_cpt, final_genome_cpt = group_counts_by_biotype(cpt, cpt_genome, biotypes)
        make_plot(labels, biotypes, final_cat_cpt, final_genome_cpt, "Biotypes", options)

    print("### End of program ###")


if __name__ == "__main__":
    main()