#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "noel & bahin"
""" ALFA provides a global overview of features distribution composing NGS dataset(s). """

import argparse
import os
import numpy
import copy
import sys
import subprocess
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colors
import matplotlib.patheffects as PathEffects
import re
from matplotlib.backends.backend_pdf import PdfPages
# To correctly embed the texts when saving plots in svg format
import matplotlib
# import progressbar
import collections
import matplotlib as mpl
import numpy as np

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
    return [ tryint(c) for c in re.split("([0-9]+)", s) ]


def required_arg(arg, aliases):
    """ Function to display the help and quit if a required argument is missing. """
    if not arg:
        print >> sys.stderr, "\nError: %s argument is missing.\n" % aliases
        parser.print_usage()
        sys.exit(1)


def existing_file(filename):
    """ Checks if filename already exists and exit if so. """
    if os.path.isfile(filename):
        sys.exit("Error: The file '" + filename + "' is about to be produced but already exists in the directory. \n### End of program")


def get_chromosome_lengths():
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
                print >> sys.stderr, "Warning: chromosome '" + chrom + "' of the chromosome lengths file does not match any chromosome name in the GTF file provided and was ignored."
        # Checking if the chromosomes from the GTF file are present in the lengths file
        for chrom in gtf_chrom_names:
            if chrom not in lengths:
                print >> sys.stderr, "Warning: at least one chromosome ('" + chrom + "') was found in the GTF file and does not match any chromosome provided in the lengths file."
                print >> sys.stderr, "\t=> All the chromosome lengths will be approximated using annotations in the GTF file."
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
    print "The chromosome lengths have been approximated using the GTF file annotations (the stop position of the last annotation of each chromosome is considered as its length)."
    return lengths


def write_feature_on_index(feat, chrom, start, stop, sign, stranded_genome_index):
    """ Write one new line in the stranded index file and, if necessary, the unstranded index file. """
    grouped_by_biotype_features = []
    for biotype, categs in feat.iteritems():
        categ_list = []
        for cat in set(categs):
            categ_list.append(cat)
        grouped_by_biotype_features.append(":".join((str(biotype), ",".join(categ_list))))

    #stranded_genome_index.write('\t'.join((chrom, start, stop, sign, '')) + '\t'.join(grouped_by_biotype_features) + '\n')
    stranded_genome_index.write(
        '\t'.join((chrom, start, stop, sign)) + '\t' + '\t'.join(grouped_by_biotype_features) + '\n')
    if unstranded_genome_index: ## MB: Why? Not always unstranded and stranded??
        #unstranded_genome_index.write('\t'.join((chrom, start, stop, '.', '')) + '\t'.join(grouped_by_biotype_features) + '\n')
        unstranded_genome_index.write(
            '\t'.join((chrom, start, stop, '.')) + '\t' + '\t'.join(grouped_by_biotype_features) + '\n')


def write_index_line(feat, chrom, start, stop, sign, fh):
    """ Write a new line in an index file. """
    # Formatting the features info
    feat_by_biotype = []
    for biot, cat in feat.iteritems():
        #feat_by_biotype.append(":".join((str(biot), ",".join(set(cat)))))
        feat_by_biotype.append(":".join((str(biot), ",".join(set(cat)))))
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


def count_genome_features(cpt, features, start, stop, coverage=1):
    """ Reads genome index and registers feature counts. """
    # If no biotype priority: category with the highest priority for each found biotype has the same weight (1/n_biotypes)
    if not biotype_prios:
        nb_biot = len(features)
        # For each categ(s)/biotype pairs
        for feat in features:
            cur_prio = 0
            # Separate categorie(s) and biotype
            try:
                biot, cats = feat.split(":")
            # Error if the feature is "antisense": update the "antisense/antisense" counts
            except ValueError:
                try:
                    cpt[("antisense", "antisense")] += (int(stop) - int(start)) * coverage / float(nb_biot)
                except KeyError:
                    cpt[("antisense", "antisense")] = (int(stop) - int(start)) * coverage / float(nb_biot)
                return None
            # Browse the categories and get only the one(s) with highest priority
            for cat in cats.split(","):
                try:
                    prio = prios[cat]
                except KeyError:
                    # TODO: Find a way to add unknown categories
                    if cat not in unknown_cat:
                        print >> sys.stderr, "Warning: Unknown categorie '%s' found and ignored.\n" % cat,
                    unknown_cat.add(cat)
                    continue
                # Check if the category has a highest priority than the current one
                if prio > cur_prio:
                    cur_prio = prio
                    cur_cat = {cat}
                if prio == cur_prio:
                    cur_cat.add(cat)
            # Increase each counts by the coverage divided by the number of categories and biotypes
            nb_cat = len(cur_cat)
            for cat in cur_cat:
                try:
                    cpt[(cat, biot)] += (int(stop) - int(start)) * coverage / (float(nb_biot * nb_cat))
                except KeyError:
                    cpt[(cat, biot)] = (int(stop) - int(start)) * coverage / (float(nb_biot * nb_cat))
    else:
        # TODO: Add an option to provide biotype priorities and handle it
        pass

'''
def add_info(cpt, feat_values, start, stop, chrom=None, unstranded_genome_index=None, stranded_genome_index=None,
             biotype_prios=None, coverage=1, categ_prios=None):
    """
    From an annotated genomic interval (start/stop positions and associated feature: one or more category(ies)/biotype pair(s) )
    - If a file is provided: write the  information at the end of the index file being generated;
    - else: browse the features and update the counts of categories/biotypes found in the genome.
    """
    ## Writing in the file if it was provided
    if stranded_genome_index:
        # If only one strand has a feature, this feature will directly be written on the unstranded index
        if feat_values[0] == {}:
            stranded_genome_index.write('\t'.join((chrom, start, stop, '+', 'antisense\n')))
        elif feat_values[1] == {}:
            stranded_genome_index.write('\t'.join((chrom, start, stop, '-', 'antisense\n')))

        write_feature_on_index(feat_values[0], chrom, start, stop, '+', stranded_genome_index)
        write_feature_on_index(feat_values[1], chrom, start, stop, '-', stranded_genome_index)
        unstranded_feat = dict(feat_values[0], **feat_values[1])
        for name in set(feat_values[0]) & set(feat_values[1]):
            unstranded_feat[name] += feat_values[0][name]
        write_feature_on_index(unstranded_feat, chrom, start, stop, '.', unstranded_genome_index)

        if feat_values[0] == {}:
            # An interval with no feature corresponds to a region annotated only on the reverse strand: update 'antisense' counts
            stranded_genome_index.write('\t'.join((chrom, start, stop, '+', 'antisense\n')))
            #write_feature_on_index(feat_values[1], chrom, start, stop, '-', stranded_genome_index, unstranded_genome_index)
            write_feature_on_index(feat_values[1], chrom, start, stop, '-', stranded_genome_index)
            write_feature_on_index(feat_values[1], chrom, start, stop, '.', unstranded_genome_index)
        elif feat_values[1] == {}:
                #write_feature_on_index(feat_values[0], chrom, start, stop, '+', stranded_genome_index, unstranded_genome_index)
                write_feature_on_index(feat_values[0], chrom, start, stop, '+', stranded_genome_index)
                write_feature_on_index(feat_values[0], chrom, start, stop, '.', unstranded_genome_index)
                stranded_genome_index.write('\t'.join((chrom, start, stop, '-', 'antisense\n')))
        # Else, the unstranded index should contain the union of plus and minus features
        else:
            write_feature_on_index(feat_values[0], chrom, start, stop, '+', stranded_genome_index)
            write_feature_on_index(feat_values[1], chrom, start, stop, '-', stranded_genome_index)
            unstranded_feat = dict(feat_values[0], **feat_values[1])
            for name in set(feat_values[0]) & set(feat_values[1]):
                unstranded_feat[name] += feat_values[0][name]
            write_feature_on_index(unstranded_feat, chrom, start, stop, '.', unstranded_genome_index)

    ## Increasing category counter(s)
    else: ##MB Why increase if file not provided??
        # Default behavior if no biotype priorities : category with the highest priority for each found biotype has the same weight (1/n_biotypes)
        if not biotype_prios:
            nb_feat = len(feat_values) ## MB: nb biotypes??
            # For every categ(s)/biotype pairs
            for feat in feat_values: ## MB: for each biotype
                cur_prio = 0
                selected_categ = []
                # Separate categorie(s) and biotype
                try:
                    biot, cats = feat.split(":")
                # Error if feature equal "antisense" : update the 'antisense/antisense' counts
                except ValueError:
                    try:
                        cpt[(feat, feat)] += (int(stop) - int(start)) * coverage / float(nb_feat) ## MB: clearly write 'antisense'??
                    except: ## MB: add a KeyError exception??
                        cpt[(feat, feat)] = (int(stop) - int(start)) * coverage / float(nb_feat)
                    return
                # Browse the categories and get only the one(s) with highest priority
                for cat in cats.split(','):
                    try:
                        prio = prios[cat]
                    except: ## MB: KeyError??
                        # TODO Find a way to add unknown categories
                        if cat not in unknown_feature:
                            print >> sys.stderr, "Warning: Unknown categorie %s found and ignored.\n" % cat,
                            unknown_feature.append(cat)
                        continue
                    if prio > cur_prio:
                        cur_prio = prio
                        selected_categ = [cat]
                    if prio == cur_prio:
                        if cat != selected_categ: ## MB: work with set??
                            try:
                                if cat not in selected_categ:
                                    selected_categ.append(cat)
                            except TypeError: ## What TypeError??
                                selected_categ = [selected_categ, cat]
                # Increase each counts by the coverage divided by the number of categories and biotypes
                nb_cats = len(selected_categ)
                for cat in selected_categ:
                    try:
                        cpt[(cat, biot)] += (int(stop) - int(start)) * coverage / (float(nb_feat * nb_cats))
                    except KeyError:
                        cpt[(cat, biot)] = (int(stop) - int(start)) * coverage / (float(nb_feat * nb_cats))
                    # else :
                    # cpt[(cats,biot)] = (int(stop) - int(start)) / float(nb_feat) * coverage
        # Else, apply biotype selection according to the priority set
        else:
            # TODO Add an option to pass biotype priorities and handle it
            pass
'''

def register_interval(features_dict, chrom, stranded_index_fh, unstranded_index_fh):
    """ Write the interval features info into the genome index files. """
    # Adding the chromosome to the list if not present
    if chrom not in index_chrom_list:
        index_chrom_list.append(chrom)
    # Writing the chromosome in the index file
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
            for biotype, categ in features_plus.iteritems():
                if set(categ) == {"gene", "transcript"}:
                    if "exon" in prev_features_plus[biotype]:
                        categ.append("intron")
                else:
                    continue
            for biotype, categ in features_minus.iteritems():
                if set(categ) == {"gene", "transcript"}:
                    if "exon" in prev_features_minus[biotype]:
                        categ.append("intron")
                else:
                    continue



def generate_genome_index(annotation, unstranded_genome_index, stranded_genome_index, chrom_sizes):
    """ Create an index of the genome annotations and save it in a file. """
    # Initializations
    intervals_dict = {}
    max_value = -1
    prev_chrom = ""
    i = 0  # Line counter
    # Write the chromosome lengths as comment lines before the genome index
    with open(unstranded_genome_index, "w") as unstranded_index_fh, open(stranded_genome_index, "w") as stranded_index_fh:
        for key, value in chrom_sizes.items():
            unstranded_index_fh.write("#%s\t%s\n" % (key, value))
            stranded_index_fh.write("#%s\t%s\n" % (key, value))
    # Progress bar to track the genome indexes creation
    nb_lines = sum(1 for _ in open(annotation))
    # pbar = progressbar.ProgressBar(widgets=["Indexing the genome ", progressbar.Percentage(), " ", progressbar.Bar(), progressbar.Timer()], maxval=nb_lines).start()
    # Browsing the GTF file and writing into genome index files
    with open(annotation, "r") as gtf_fh:
        for line in gtf_fh:
            i += 1
            # Update the progressbar every 1k lines
            # if i % 1000 == 1:
                # pbar.update(i)
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
                biotype = line_split[8].split("biotype")[1].split(";")[0].strip('" ')
                # Registering stored features info in the genome index file(s) if the new line concerns a new chromosome or the new line concerns an annotation not overlapping previously recorded ones
                if start > max_value or chrom != prev_chrom:
                    # Write the previous features
                    if intervals_dict:
                        register_interval(intervals_dict, prev_chrom, stranded_genome_index, unstranded_genome_index)
                    prev_chrom = chrom
                    # (Re)Initializing the intervals info dict
                    intervals_dict = {strand: {start: {biotype: [cat]}, stop: {}}, antisense: {start: {}, stop: {}}}
                    max_value = stop

                # Update the dictionary which represents intervals for every distinct annotation
                else:
                    # Storing the intervals on the strand of the current feature
                    stranded_intervals = intervals_dict[strand]
                    added_info = False # Variable to know if the features info were already added
                    # Browsing the existing boundaries
                    for boundary in sorted(stranded_intervals):
                        # While the GTF line start is after the browsed boundary: keep the browsed boundary features info in case the GTF line start is before the next boundary
                        if boundary < start:
                            stored_feat_strand, stored_feat_antisense = [dict(stranded_intervals[boundary]), dict(intervals_dict[antisense][boundary])]

                        # The GTF line start is already an existing boundary: store the existing features info (to manage after the GTF line stop) and update it with the GTF line features info
                        elif boundary == start:
                            stored_feat_strand, stored_feat_antisense = [dict(stranded_intervals[boundary]), dict(intervals_dict[antisense][boundary])]
                            # Adding the GTF line features info to the interval
                            try:
                                stranded_intervals[boundary][biotype] = stranded_intervals[boundary][biotype] + [cat]
                            except KeyError: # If the GTF line features info regard an unregistered biotype
                                stranded_intervals[boundary][biotype] = [cat]
                            added_info = True # The features info were added

                        # The browsed boundary is after the GTF line start: add the GTF line features info
                        elif boundary > start:
                            # Create a new boundary for the GTF line start if necessary (if it is between 2 existing boundaries, it was not created before)
                            if not added_info:
                                stranded_intervals[start] = copy.deepcopy(stored_feat_strand)
                                #stranded_intervals[start][biotype] = [cat]
                                try:
                                    stranded_intervals[start][biotype].append(cat)
                                except KeyError:
                                    stranded_intervals[start][biotype] = [cat]
                                intervals_dict[antisense][start] = copy.deepcopy(stored_feat_antisense)
                                added_info = True # The features info were added
                            # While the browsed boundary is before the GTF line stop: store the existing features info (to manage after the GTF line stop) and update it with the GTF line features info
                            if boundary < stop:
                                stored_feat_strand, stored_feat_antisense = [dict(stranded_intervals[boundary]), dict(intervals_dict[antisense][boundary])]
                                try:
                                    stranded_intervals[boundary][biotype] = stranded_intervals[boundary][biotype] + [cat]
                                except KeyError:
                                    stranded_intervals[boundary][biotype] = [cat]
                            # The GTF line stop is already exists, nothing more to do, the GTF line features info are integrated
                            elif boundary == stop:
                                break
                            # The browsed boundary is after the GTF line stop: create a new boundary for the GTF line stop (with the stored features info)
                            else: # boundary > stop
                                stranded_intervals[stop] = copy.deepcopy(stored_feat_strand)
                                intervals_dict[antisense][stop] = copy.deepcopy(stored_feat_antisense)
                                break # The GTF line features info are integrated
                    # If the GTF line stop is after the last boundary, extend the dictionary
                    if stop > max_value:
                        max_value = stop
                        stranded_intervals[stop] = {}
                        intervals_dict[antisense][stop] = {}

        # Store the categories of the last chromosome
        register_interval(intervals_dict, chrom, stranded_genome_index, unstranded_genome_index)
    # pbar.finish()


def generate_bedgraph_files(sample_labels, bam_files):
    """ Creates BedGraph files from BAM ones and return filenames and labels. """
    #sample_files = []
    #sample_labels = []
    # Progress bar to track the BedGraph file creation
    #pbar = progressbar.ProgressBar(widgets=["Generating the BedGraph files ", progressbar.Percentage(), progressbar.Bar(), progressbar.Timer()], max_value=len(bam_files)+1).start()
    # pbar = progressbar.ProgressBar(widgets=["Generating the BedGraph files ", progressbar.Percentage(), progressbar.Bar(), progressbar.Timer()], maxval=len(sample_labels)+1).start()
    n = 1
    # pbar.update(n)
    #for n in range(0, len(bam_files), 2):
    for sample_label, bam_file in zip(sample_labels, bam_files):
        # Get the label for this sample
        #sample_labels.append(bam_files[n + 1])
        # Modify it to contain only alphanumeric characters (avoid files generation with dangerous names)
        #modified_label = "_".join(re.findall(r"[\w']+", bam_files[n + 1]))
        if options.strandness in ["forward", "fr-firststrand"]:
            #subprocess.call("bedtools genomecov -bg -split -strand + -ibam " + bam_files[n] + " > " + modified_label + ".plus.bedgraph", shell=True)
            subprocess.call("bedtools genomecov -bg -split -strand + -ibam " + bam_file + " > " + sample_label + ".plus.bedgraph", shell=True)
            # pbar.update(n + 0.5)
            #subprocess.call("bedtools genomecov -bg -split -strand - -ibam " + bam_files[n] + " > " + modified_label + ".minus.bedgraph", shell=True)
            subprocess.call("bedtools genomecov -bg -split -strand - -ibam " + bam_file + " > " + sample_label + ".minus.bedgraph", shell=True)
            # pbar.update(n + 0.5)
        elif options.strandness in ["reverse", "fr-secondstrand"]:
            #subprocess.call("bedtools genomecov -bg -split -strand - -ibam " + bam_files[n] + " > " + modified_label + ".plus.bedgraph", shell=True)
            subprocess.call("bedtools genomecov -bg -split -strand - -ibam " + bam_file + " > " + sample_label + ".plus.bedgraph", shell=True)
            # pbar.update(n + 0.5)
            #subprocess.call("bedtools genomecov -bg -split -strand + -ibam " + bam_files[n] + " > " + modified_label + ".minus.bedgraph", shell=True)
            subprocess.call("bedtools genomecov -bg -split -strand + -ibam " + bam_file + " > " + sample_label + ".minus.bedgraph", shell=True)
            # pbar.update(n + 0.5)
        else:
            #subprocess.call("bedtools genomecov -bg -split -ibam " + bam_files[n] + " > " + modified_label + ".bedgraph", shell=True)
            subprocess.call("bedtools genomecov -bg -split -ibam " + bam_file + " > " + sample_label + ".bedgraph", shell=True)
            # pbar.update(n + 1)
        #sample_files.append(modified_label)
    # pbar.finish()
    #return sample_files, sample_labels
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
        # strand information can not be found in the file file header
        except IndexError:
            pass
    gtf_chrom = splitline[0]
    gtf_features = splitline[4:]
    gtf_start, gtf_stop = map(int, splitline[1:3])
    return endGTF


#def read_counts(counts_files):
def read_counts(sample_labels, counts_files):
    """ Reads the counts from an input file. """
    cpt = {}
    cpt_genome = {}
    #for fcounts in counts_files:
    for sample_label, filename in zip(sample_labels, counts_files):
        #label = os.path.splitext(os.path.basename(fcounts))[0]
        #labels.append(label)
        #cpt[label] = {}
        cpt[sample_label] = {}
        #with open(fcounts, "r") as counts_fh:
        with open(filename, "r") as counts_fh:
            for line in counts_fh:
                if not line.startswith("#"):
                    feature = tuple(line.split("\t")[0].split(","))
                    #cpt[label][feature] = float(line.split("\t")[1])
                    cpt[sample_label][feature] = float(line.split("\t")[1])
                    cpt_genome[feature] = float(line.rstrip().split("\t")[2])
    #return cpt, cpt_genome, labels
    return cpt, cpt_genome


def get_chromosome_names_in_index(genome_index):
    """ Returns the chromosome names in a genome index file. """
    with open(genome_index, "r") as index_fh:
        for line in index_fh:
            if not line.startswith("#") and (line.split("\t")[0] not in index_chrom_list):
                index_chrom_list.append(line.split("\t")[0])
    return index_chrom_list


def read_index():
    """ Parse index files info (chromosomes list, lengths and genome features). """
    with open(genome_index, "r") as index_fh:
        for line in index_fh:
            if line.startswith("#"):
                lengths[line.split("\t")[0][1:]] = int(line.split("\t")[1])
            else:
                chrom = line.split("\t")[0]
                if chrom not in index_chrom_list:
                    index_chrom_list.append(chrom)
                count_genome_features(cpt_genome, line.rstrip().split("\t")[4:], line.split("\t")[1], line.split("\t")[2])


#def intersect_bedgraphs_and_index_to_counts_categories(sample_files, sample_labels, prios, genome_index, biotype_prios=None): ## MB: To review
def intersect_bedgraphs_and_index_to_counts_categories(sample_labels, bedgraph_files, biotype_prios=None): ## MB: To review
    global gtf_line, gtf_chrom, gtf_start, gtf_stop, gtf_cat, endGTF
    unknown_chrom = []
    cpt = {}  # Counter for the nucleotides in the BAM input file(s)
    #for n in range(len(sample_files)):
    if bedgraph_files == []:
        bedgraph_files = sample_labels
    for sample_label, bedgraph_basename in zip(sample_labels, bedgraph_files):
        #sample_file = sample_files[n]
        #sample_name = sample_labels[n]
        # Initializing the category counter dict for this sample and the number of lines to process for the progress bar
        #init_dict(cpt, sample_name, {})
        init_dict(cpt, sample_label, {})
        bg_extension = ".bedgraph"
        if options.strandness == "unstranded":
            strands = [("", ".")]
            #nb_lines = sum(1 for _ in open(sample_file + ".bedgraph"))
            try :
                nb_lines = sum(1 for _ in open(bedgraph_basename + bg_extension))
            except IOError:
                bg_extension = "bg"
                nb_lines = sum(1 for _ in open(bedgraph_basename + bg_extension))
        else:
            strands = [(".plus", "+"), (".minus", "-")]
            #nb_lines = sum(1 for _ in open(sample_file + ".plus.bedgraph")) + sum(1 for _ in open(sample_file + ".minus.bedgraph"))
            try:
                nb_lines = sum(1 for _ in open(bedgraph_basename + ".plus" + bg_extension)) + sum(1 for _ in open(bedgraph_basename + ".minus" + bg_extension))
            except IOError:
                nb_lines = sum(1 for _ in open(bedgraph_basename + ".plus" + bg_extension)) + sum(
                    1 for _ in open(bedgraph_basename + ".minus" + bg_extension))

        # Progress bar to track the BedGraph and index intersection
        #pbar = progressbar.ProgressBar(widgets=["Processing " + sample_file + " ", progressbar.Percentage(), progressbar.Bar(), progressbar.Timer()], max_value=nb_lines).start()
        # pbar = progressbar.ProgressBar(widgets=["Processing " + sample_label + " ", progressbar.Percentage(), progressbar.Bar(), progressbar.Timer()], maxval=nb_lines).start()
        i = 0

        # Intersecting the BedGraph and index files
        for strand, sign in strands:
            prev_chrom = ""
            endGTF = False  # Reaching the next chr or the end of the GTF index
            intergenic_adds = 0.0
            #with open(sample_file + strand + ".bedgraph", "r") as bedgraph_fh:
            with open(bedgraph_basename + strand + bg_extension, "r") as bedgraph_fh:
                # Running through the BedGraph file
                for bam_line in bedgraph_fh:
                    i += 1
                    # if i % 10000 == 0:
                        # pbar.update(i)
                    # Getting the BAM line info
                    bam_chrom = bam_line.split("\t")[0]
                    bam_start, bam_stop, bam_cpt = map(float, bam_line.split("\t")[1:4])
                    # Skip the line if the chromosome is not in the index
                    if bam_chrom not in index_chrom_list:
                        if bam_chrom not in unknown_chrom:
                            unknown_chrom.append(bam_chrom)
                            print "\r                          \r Chromosome '" + bam_chrom + "' not found in index." # MB: to adapt with the progress bar
                        continue
                    # If this is a new chromosome (or the first one)
                    if bam_chrom != prev_chrom:
                        intergenic_adds = 0.0
                        # (Re)opening the GTF index and looking for the first line of the matching chr
                        try:
                            gtf_index_file.close()
                        except UnboundLocalError:
                            pass
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
                            #cpt[sample_name][("intergenic", "intergenic")] += (min(bam_stop,
                            cpt[sample_label][("intergenic", "intergenic")] += (min(bam_stop,
                                                                                    gtf_start) - bam_start) * bam_cpt
                        except KeyError:
                            #cpt[sample_name][("intergenic", "intergenic")] = (min(bam_stop,
                            cpt[sample_label][("intergenic", "intergenic")] = (min(bam_stop,
                                                                                   gtf_start) - bam_start) * bam_cpt
                        # Go to next line if the BAM interval was totally upstream the first GTF annotation, carry on with the remaining part otherwise
                        if endGTF or (bam_stop <= gtf_start):
                            continue
                        else:
                            bam_start = gtf_start

                    # We can start the crossover
                    while not endGTF:
                        # Update category counter
                        #add_info(cpt[sample_name], gtf_features, bam_start, min(bam_stop, gtf_stop), coverage=bam_cpt)
                        #count_genome_features(cpt[sample_name], gtf_features, bam_start, min(bam_stop, gtf_stop), coverage=bam_cpt)
                        count_genome_features(cpt[sample_label], gtf_features, bam_start, min(bam_stop, gtf_stop), coverage=bam_cpt)
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
                            #cpt[sample_name][("intergenic", "intergenic")] += (bam_stop - bam_start) * bam_cpt
                            cpt[sample_label][("intergenic", "intergenic")] += (bam_stop - bam_start) * bam_cpt
                        except KeyError:
                            #cpt[sample_name][("intergenic", "intergenic")] = (bam_stop - bam_start) * bam_cpt
                            cpt[sample_label][("intergenic", "intergenic")] = (bam_stop - bam_start) * bam_cpt
                gtf_index_file.close()
        # pbar.finish()
    return cpt


def write_counts_in_files(cpt, genome_counts):
    """ Writes the biotype/category counts in an output file. """
    for sample_label, counters in cpt.items():
        sample_label = "_".join(re.findall(r"[\w\-']+", sample_label))
        with open(sample_label + ".feature_counts.tsv", "w") as output_fh:
            output_fh.write("#Category,biotype\tCounts_in_bam\tSize_in_genome\n")
            for features_pair, counts in counters.items():
                output_fh.write("%s\t%s\t%s\n" % (",".join(features_pair), counts, genome_counts[features_pair]))


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
                for key, value in cpt[f].items():
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
                    for key, value in cpt[f].items():
                        if key[1] == final_biot:
                            tot += value
                            # if key[1] != 'antisense':
                            tot_genome += cpt_genome[key]
            except:
                for key, value in cpt[f].items():
                    if key[1] == biot:
                        tot += value
                        tot_genome += cpt_genome[key]
            if tot != 0:
                final_cpt[f][biot] = tot
                final_genome_cpt[biot] = tot_genome
    return final_cpt, final_genome_cpt


# def get_cmap(N):
# '''Returns a function that maps each index in 0, 1, ... N-1 to a distinct
# RGB color.'''
# color_norm  = colors.Normalize(vmin=0, vmax=N-1)
# scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv')
# def map_index_to_rgb_color(index):
# return scalar_map.to_rgba(index)
# return map_index_to_rgb_color

def one_sample_plot(ordered_categs, percentages, enrichment, n_cat, index, index_enrichment, bar_width, counts_type,
                    title, sample_labels):
    ### Initialization
    fig = plt.figure(figsize=(13, 9))
    ax1 = plt.subplot2grid((2, 4), (0, 0), colspan=2)
    ax2 = plt.subplot2grid((2, 4), (1, 0), colspan=2)
    cmap = plt.get_cmap("Spectral")
    cols = [cmap(x) for x in xrange(0, 256, 256 / n_cat)]
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
    pielabels = [ordered_categs[i] if percentages[i] > 0.025 else "" for i in xrange(n_cat)]
    sum_enrichment = numpy.sum(enrichment)
    pielabels_enrichment = [ordered_categs[i] if enrichment[i] / sum_enrichment > 0.025 else "" for i in xrange(n_cat)]
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


#def make_plot(ordered_categs, sample_names, categ_counts, genome_counts, pdf, counts_type, threshold, title=None,
def make_plot(sample_labels, ordered_categs, categ_counts, genome_counts, pdf, counts_type, threshold, title=None, svg=None, png=None,
              categ_groups=[]):  # MB: to review
    # From ordered_categs, keep only the features (categs or biotypes) that we can find in at least one sample.
    existing_categs = set()
    for sample in categ_counts.values():
        existing_categs |= set(sample.keys())
    ordered_categs = filter(existing_categs.__contains__, ordered_categs)
    xlabels = [cat if len(cat.split("_")) == 1 else "\n".join(cat.split("_")) if cat.split("_")[0] != 'undescribed' else "\n".join(["und.",cat.split("_")[1]]) for cat in ordered_categs]
    n_cat = len(ordered_categs)
    #n_exp = len(sample_names)
    nb_samples = len(categ_counts)
    ##Initialization of the matrix of counts (nrow=nb_experiements, ncol=nb_categorie)
    #counts = numpy.matrix(numpy.zeros(shape=(n_exp, n_cat)))
    counts = numpy.matrix(numpy.zeros(shape=(nb_samples, n_cat)))
    '''
    for exp in xrange(len(sample_names)):
        for cat in xrange(len(ordered_categs)):
            try:
                counts[exp, cat] = categ_counts[sample_names[exp]][ordered_categs[cat]]
            except:
                pass
    '''
    for sample_label in sample_labels:
        for cat in xrange(len(ordered_categs)):
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
    if "antisense" in ordered_categs:
        antisense_pos = ordered_categs.index("antisense")
        sizes[antisense_pos] = 1e-100
    for cpt in xrange(len(sizes)):
        sizes[cpt] /= float(sizes_sum)

    ## Create array which contains the percentage of reads in each categ for every sample
    percentages = numpy.array(counts / numpy.sum(counts, axis=1))
    ## Create the enrichment array (counts divided by the categorie sizes in the genome)
    enrichment = numpy.array(percentages / sizes)
    if "antisense_pos" in locals():
        '''
        for i in xrange(len(sample_names)):
            enrichment[i][antisense_pos] = 0
        '''
        for n in xrange(nb_samples):
            enrichment[n][antisense_pos] = 0

    # enrichment=numpy.log(numpy.array(percentages/sizes))
    #for exp in xrange(n_exp):
    for n in xrange(nb_samples):
        for i in xrange(n_cat):
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
        cmap = plt.get_cmap("Set3", nb_samples)
        cmap = [cmap(i) for i in xrange(nb_samples)]

    ## Parameters for the plot
    opacity = 1
    # Create a vector which contains the position of each bar
    index = numpy.arange(n_cat)
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
    # Store the bars objects for enrichment plot
    rects = []
    # For each sample/experiment
    #for i in range(n_exp):
    for sample_label in sample_labels:
        # First barplot: percentage of reads in each categorie
        n = sample_labels.index(sample_label)
        #ax1.bar(index + i * bar_width, percentages[i], bar_width,
        ax1.bar(index + n * bar_width, percentages[n], bar_width,
                alpha=opacity,
                #color=cmap[i],
                color=cmap[n],
                #label=sample_names[i], edgecolor="#FFFFFF", lw=0)
                label=sample_label, edgecolor="#FFFFFF", lw=0)
        # Second barplot: enrichment relative to the genome for each categ
        # (the reads count in a categ is divided by the categ size in the genome)
        rects.append(ax2.bar(index + n * bar_width, enrichment[n], bar_width,
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
    if threshold:
        threshold_bottom = -abs(float(threshold[0])) + 1
        threshold_top = float(threshold[1]) - 1

        #for i in xrange(n_exp):
        for n in xrange(nb_samples):
            for y in xrange(n_cat):
                #val = enrichment[i][y]
                val = enrichment[n][y]
                if not numpy.isnan(val) and not (threshold_bottom < val < threshold_top):
                    #rect = rects[i][y]
                    rect = rects[n][y]
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
        margin_top, margin_bottom = [max(margin_top, margin_bottom) for i in xrange(2)]
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
    yticks = [yticks[i] - 1 if yticks[i] > 9 else yticks[i] + 1 if yticks[i] < -9 else yticks[i] for i in
              xrange(len(yticks))]
    ax2.set_yticks(yticks)
    ax2.set_yticklabels([str(int(i + 1)) + "$^{+1}$" if i > 0 and i % 1 == 0 else str(
        i + 1) + "$^{+1}$" if i > 0 else 1 if i == 0 else str(
        int(-(i - 1))) + "$^{-1}$" if i < 0 and i % 1 == 0 else str(-(i - 1)) + "$^{-1}$" for i in ax2.get_yticks()])
    # ax2.set_yticklabels([i+1 if i>0 else 1 if i==0 else "$\\frac{1}{%s}$" %-(i-1) for i in ax2.get_yticks()])
    # Change appearance of 'antisense' bars on enrichment plot since we cannot calculate an enrichment for this artificial category
    if "antisense_pos" in locals():  # ax2.text(antisense_pos+bar_width/2,ax2.get_ylim()[1]/10,'NA')
        #for i in xrange(n_exp):
        for n in xrange(nb_samples):
            #rect = rects[i][antisense_pos]
            rect = rects[n][antisense_pos]
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
    # Add text for features absent in sample
    #for i in xrange(n_exp):
    for n in xrange(nb_samples):
        for y in xrange(n_cat):
            #if percentages[i][y] == 0:
            if percentages[n][y] == 0:
                txt = ax1.text(y + bar_width * (n + 0.5), 0.02, "Abs.", rotation="vertical", color=cmap[n],
                               horizontalalignment="center", verticalalignment="bottom")
                txt.set_path_effects([PathEffects.Stroke(linewidth=0.5), PathEffects.Normal()])
            #elif enrichment[i][y] == 0:
            elif enrichment[n][y] == 0:
                #rects[i][y].set_linewidth(1)
                rects[n][y].set_linewidth(1)

    # Remove top/right/bottom axes
    for ax in [ax1, ax2]:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(axis="x", which="both", bottom="on", top="off", labelbottom="on")
        ax.tick_params(axis="y", which="both", left="on", right="off", labelleft="on")


    ### Add second axis with categ groups
    annotate_group(categ_groups, label=None, ax=ax1)
    annotate_group(categ_groups, label=None, ax=ax2)

    ### Adjust figure margins to
    if counts_type.lower() == "categories":
        plt.tight_layout(h_pad=5.0)
        fig.subplots_adjust(bottom=0.1)
    else:
        plt.tight_layout()

    ## Showing the plot
    if pdf:  ## TODO: allow several output formats
        pdf.savefig()
        plt.close()
    elif svg:
        if svg == True:
            plt.savefig(counts_type + ".svg")
        else:
            if os.path.splitext(svg)[1] == ".svg":
                plt.savefig(".".join((os.path.splitext(svg)[0], counts_type, "svg")))
            else:
                plt.savefig(".".join((svg, counts_type, "svg")))
    elif png:
        if png == True:
            plt.savefig(counts_type + ".png")
        else:
            if os.path.splitext(png)[1] == ".png":
                plt.savefig(".".join((os.path.splitext(png)[0], counts_type, "png")))
            else:
                plt.savefig(".".join((png, counts_type, "png")))
    else:
        plt.show()
    ## Save on disk it as a png image
    # fig.savefig('image_output.png', dpi=300, format='png', bbox_extra_artists=(lgd,), bbox_inches='tight')

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
        for name, coord in grp.items():
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
        for feature, count in cpt[sample].items():
            if feature[1] == selected_biotype:
                filtered_cpt[sample][feature[0]] = count
    return filtered_cpt


def unnecessary_param(param, message):
    """ Function to display a warning on the standard error. """
    if param:
        print >> sys.stderr, message


def usage_message():
    return """
    Generate genome indexes:
         python ALFA.py -a GTF_FILE  [-g GENOME_INDEX]
                                         [--chr_len CHR_LENGTHS_FILE]
    Process BAM file(s):
         python ALFA.py -g GENOME_INDEX -i BAM1 LABEL1 [BAM2 LABEL2 ...]
                                         [--bedgraph] [-s STRAND]
                                         [-n] [--pdf output.pdf]
                                         [-d {1,2,3,4}] [-t YMIN YMAX]
    Index genome + process BAM:
         python ALFA.py -a GTF_FILE [-g GENOME_INDEX]
                            -i BAM1 LABEL1 [BAM2 LABEL2 ...]
                            [--chr_len CHR_LENGTHS_FILE]
                            [--bedgraph][-s STRAND]
                            [-n] [--pdf output.pdf]
                            [-d {1,2,3,4}] [-t YMIN YMAX]

    Process previously created ALFA counts file(s):
         python ALFA.py -c COUNTS1 [COUNTS2 ...]
                            [-s STRAND]
                            [-n] [--pdf output.pdf]
                            [-d {1,2,3,4}] [-t YMIN YMAX]

        """

##########################################################################
#                           MAIN                                         #
##########################################################################


if __name__ == "__main__":

    #### Parse command line arguments and store them in the variable options
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, usage=usage_message())
    parser.add_argument("--version", action="version", version="version 1.0",
                        help="Show ALFA version number and exit\n\n-----------\n\n")
    # Options regarding the index
    parser.add_argument("-g", "--genome_index",
                        help="Genome index files path and basename for existing index, or path and basename for new index creation\n\n")
    parser.add_argument("-a", "--annotation", metavar="GTF_FILE", help="Genomic annotations file (GTF format)\n\n")
    parser.add_argument("--chr_len", help="Tabulated file containing chromosome names and lengths\n\n-----------\n\n")

    # Options regarding the intersection step
    #parser.add_argument('-i', '--input', '--bam', dest='input', metavar=('BAM_FILE1 LABEL1', ""), nargs='+',
                        #help='Input BAM file(s) and label(s). The BAM files must be sorted by position.\n\n')
    parser.add_argument("--bam", metavar=("BAM_FILE1 LABEL1", ""), nargs="+",
                        help="Input BAM file(s) and label(s). The BAM files must be sorted by position.\n\n") ## MB: position AND chr??
    # parser.add_argument('--bedgraph', action='store_const',default = False, const = True, help="Use this options if your input file(s) is(are) already in bedgraph format\n\n")
    #parser.add_argument('--bedgraph', dest='bedgraph', action='store_const', default=False, const=True,
                        #help="Use this options if your input file(s) is(are) already in bedgraph format\n\n")
    parser.add_argument("--bedgraph", metavar=("BEDGRAPH_FILE1 LABEL1", ""), nargs="+",
                        help="Use this options if your input is/are BedGraph file(s).\n\n")
    parser.add_argument("-c", "--counts", metavar=("COUNTS_FILE", ""), nargs="+",
                        help="Use this options instead of '-i/--input' to provide ALFA counts files as input \ninstead of bam/bedgraph files.\n\n")
    #parser.add_argument('-s', '--strandness', dest="strandness", nargs=1, action='store', default=['unstranded'], choices=['unstranded', 'forward', 'reverse', 'fr-firststrand', 'fr-secondstrand'], metavar="", help="Library orientation. Choose within: 'unstranded', 'forward'/'fr-firststrand' \nor 'reverse'/'fr-secondstrand'. (Default: 'unstranded')\n\n-----------\n\n")
    parser.add_argument("-s", "--strandness", default="unstranded",
                        choices=["unstranded", "forward", "reverse", "fr-firststrand", "fr-secondstrand"], metavar="",
                        help="Library orientation. Choose within: 'unstranded', "
                             "'forward'/'fr-firststrand' \nor 'reverse'/'fr-secondstrand'. "
                             "(Default: 'unstranded')\n\n-----------\n\n")

    # Options regarding the plot
    parser.add_argument("--biotype_filter", help=argparse.SUPPRESS)  # "Make an extra plot of categories distribution using only counts of the specified biotype." ## MB: TO DISPLAY (no suppress)
    parser.add_argument("-d", "--categories_depth", type=int, default=3, choices=range(1, 5),
                        help="Use this option to set the hierarchical level that will be considered in the GTF file (default=3): \n(1) gene,intergenic; \n(2) intron,exon,intergenic; \n(3) 5'UTR,CDS,3'UTR,intron,intergenic; \n(4) start_codon,5'UTR,CDS,3'UTR,stop_codon,intron,intergenic. \n\n")
    parser.add_argument("--pdf", nargs="?", const="ALFA_plots.pdf",
                        help="Save produced plots in PDF format at specified path ('categories_plots.pdf' if no argument provided).\n\n")
    parser.add_argument("--png", nargs="?", const="ALFA_plots.png",
                        help="Save produced plots in PNG format with provided argument as basename \nor 'categories.png' and 'biotypes.png' if no argument provided.\n\n")
    parser.add_argument("--svg", nargs="?", const="ALFA_plots.svg",
                        help="Save produced plots in SVG format with provided argument as basename \nor 'categories.svg' and 'biotypes.svg' if no argument provided.\n\n")
    parser.add_argument("-n", "--no_display", action="store_const", const=True, default=False, help="Do not display plots.\n\n") # We have to add "const=None" to avoid a bug in argparse
    parser.add_argument("-t", "--threshold", dest="threshold", nargs=2, metavar=("ymin", "ymax"), type=float,
                        help="Set axis limits for enrichment plots.\n\n")

    if len(sys.argv) == 1:
        parser.print_usage()
        sys.exit(1)

    options = parser.parse_args()

    '''
    # Booleans for steps to be executed
    make_index = False
    intersect_reads = False
    process_counts = False

    #### Check arguments conformity and define which steps have to be performed
    print "### Checking parameters"
    if options.counts:
        # Aucun autre argument requis, precise that the other won't be used (if this is true!!)
        # Vérifier extension input

        # Action: Only do the plot
        process_counts = True
    else:
        if options.annotation:
            # If '-gi' parameter is present
            if options.genome_index:
                genome_index_basename = options.genome_index
            else:
                # Otherwise the GTF filename without extension will be the basename
                genome_index_basename = options.annotation.split("/")[-1].split(".gtf")[0]
            # Check if the indexes were already created and warn the user
            if os.path.isfile(genome_index_basename + ".stranded.index"):
                if options.input:
                    print >> sys.stderr, "Warning: an index file named '%s' already exists and will be used. If you want to create a new index, please delete this file or specify an other path." % (
                    genome_index_basename + ".stranded.index")
                else:
                    sys.exit(
                        "Error: an index file named %s already exists. If you want to create a new index, please delete this file or specify an other path.\n" % (
                        genome_index_basename + ".stranded.index"))
            # Create them otherwise
            else:
                make_index = True
        # If the index is already done
        if options.input:
            # Required arguments are the input and the genome_index
            if 'genome_index_basename' not in locals():
                required_arg(options.genome_index, "-g/--genome_index")
                genome_index_basename = options.genome_index
            required_arg(options.input, "-i/--input/--bam")
            for i in xrange(0, len(options.input), 2):
                # Check whether the input file exists
                try:
                    open(options.input[i])
                except IOError:
                    sys.exit("Error: the input file " + options.input[i] + " was not found. Aborting.")
                # Check whether the input file extensions are 'bam', 'bedgraph' or 'bg' and the label extension are no
                try:
                    extension = os.path.splitext(options.input[i + 1])[1]
                    if extension in ('.bam', '.bedgraph', '.bg'):
                        sys.exit("Error: it seems input files and associated labels are not correctly provided.\n\
                        Make sure to follow the expected format : -i Input_file1 Label1 [Input_file2 Label2 ...].")
                except:
                    sys.exit(
                        "Error: it seems input files and associated labels are not correctly provided.\nMake sure to follow the expected format : -i Input_file1 Label1 [Input_file2 Label2 ...].")

            intersect_reads = True
        # Vérifier input's extension
        # TODO
    if not (options.counts or options.input or options.annotation):
        sys.exit(
            "\nError : some arguments are missing At least '-a', '-c' or '-i' is required. Please refer to help (-h/--help) and usage cases for more details.\n")
    if not options.counts:
        # Declare genome_index variables
        stranded_genome_index = genome_index_basename + ".stranded.index"
        unstranded_genome_index = genome_index_basename + ".unstranded.index"
        if options.strandness[0] == "unstranded":
            genome_index = unstranded_genome_index
        else:
            genome_index = stranded_genome_index
    '''
    #### Initialization of some variables

    # Miscellaneous variables
    reverse_strand = {"+": "-", "-": "+"}
    samples = collections.OrderedDict() # Structure: {<label>: [<filename1>(, <filename2>)]} # Example: {'Toy': ['toy.bam']}

    #Sample labels and file paths
    labels = []
    bams = []
    bedgraphs = []
    count_files = []


    # Initializing the category priority order, coding biotypes and the final list
    # prios = {"start_codon": 7, "stop_codon": 7, "five_prime_utr": 6, "three_prime_utr": 6, "UTR": 6, "CDS": 5, "exon": 4,
    #          "transcript": 3, "gene": 2, "antisense": 1, "intergenic": 0}
    prios = {"start_codon": 4, "stop_codon": 4, "five_prime_utr": 3, "three_prime_utr": 3, "UTR": 3, "CDS": 3,
             "exon": 2, "intron": 2, "transcript": 1, "gene": 1, "antisense": 0, "intergenic": -1}

    biotype_prios = None
    # biotype_prios = {"protein_coding":1, "miRNA":2}

    # Possibles groups of categories to plot
    # categs_group1 = {"start": ["start_codon"], "5UTR": ["five_prime_utr", "UTR"], "CDS": ["CDS", "exon"],
    #                  "3UTR": ["three_prime_utr"], "stop": ["stop_codon"], "introns": ["transcript", "gene"],
    #                  "intergenic": ["intergenic"], "antisense": ["antisense"]}
    # categs_group2 = {"5UTR": ["five_prime_utr", "UTR"], "CDS": ["CDS", "exon", "start_codon", "stop_codon"],
    #                  "3UTR": ["three_prime_utr"], "introns": ["transcript", "gene"], "intergenic": ["intergenic"],
    #                  "antisense": ["antisense"]}
    # categs_group3 = {"exons": ["five_prime_utr", "three_prime_utr", "UTR", "CDS", "exon", "start_codon", "stop_codon"],
    #                  "introns": ["transcript", "gene"], "intergenic": ["intergenic"], "antisense": ["antisense"]}
    # categs_group4 = {
    #     "gene": ["five_prime_utr", "three_prime_utr", "UTR", "CDS", "exon", "start_codon", "stop_codon", "transcript",
    #              "gene"], "intergenic": ["intergenic"], "antisense": ["antisense"]}


    categs_level1 = {"gene": ["five_prime_utr", "three_prime_utr", "UTR", "CDS", "exon", "intron", "start_codon",
                             "stop_codon", "transcript", "gene"],
                    "intergenic": ["intergenic"],
                    "antisense": ["antisense"]}

    categs_level2 = {"exons": ["five_prime_utr", "three_prime_utr", "UTR", "CDS", "exon", "start_codon", "stop_codon"],
                    "introns": ["intron"],
                    "undescribed_genes": ["transcript", "gene"],
                    "intergenic": ["intergenic"],
                    "antisense": ["antisense"]}

    categs_level3 = {"5UTR": ["five_prime_utr", "UTR"],
                    "CDS": ["CDS", "start_codon", "stop_codon"],
                    "3UTR": ["three_prime_utr"],
                    "undescribed_exons": ["exon"],
                    "introns": ["intron"],
                    "undescribed_genes": ["transcript", "gene"],
                    "intergenic": ["intergenic"],
                    "antisense": ["antisense"]}

    categs_level4 = {"5UTR": ["five_prime_utr", "UTR"],
                     "start": ["start_codon"],
                     "stop": ["stop_codon"],
                     "CDS_body": ["CDS"],
                     "undescribed_CDS": [], #TODO: implement CDS/undescribed_CDS distinction
                     "3UTR": ["three_prime_utr"],
                     "undescribed_exons": ["exon"],
                     "introns": ["intron"],
                     "undescribed_genes": ["transcript", "gene"],
                     "intergenic": ["intergenic"],
                     "antisense": ["antisense"]}

    # categs_groups = [categs_group4, categs_group3, categs_group2, categs_group1]  # Order and merging for the final plot
    categs_levels = [categs_level1, categs_level2, categs_level3, categs_level4]

    parent_categ_level1 = []
    parent_categ_level2 = [{"gene":[0.5,2.5]}]
    parent_categ_level3 = [{"exon":[0.5,3.5]}, {"gene":[0.5,5.5]}]
    parent_categ_level4 = [{"CDS":[1.5,4.5]}, {"exon":[0.5,6.5]},{"gene":[0.5,8.5]}]
    parent_categ_groups = [parent_categ_level1, parent_categ_level2, parent_categ_level3, parent_categ_level4]

    cat_list = ["5UTR", "start", "CDS", "CDS_body", "stop", "undescribed_CDS", "3UTR", "exons", "undescribed_exons", "introns", "gene", "undescribed_genes", "intergenic", "antisense"]

    # biotypes list
    biotypes = {"protein_coding", "polymorphic_pseudogene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene", "IG_C_gene",
                "IG_D_gene", "IG_J_gene", "IG_V_gene", "3prime_overlapping_ncrna", "lincRNA", "macro_lncRNA", "miRNA",
                "misc_RNA", "Mt_rRNA", "Mt_tRNA", "processed_transcript", "ribozyme", "rRNA", "scaRNA", "sense_intronic",
                "sense_overlapping", "snoRNA", "snRNA", "sRNA", "TEC", "vaultRNA", "antisense",
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
                       "IG": ["IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene"], \
                       "MT_RNA": ["Mt_rRNA", "Mt_tRNA"], \
                       "ncRNA": ["lincRNA", "macro_lncRNA", "3prime_overlapping_ncrna", "ncRNA"], \
                       "others": ["misc_RNA", "processed_transcript", "ribozyme", "scaRNA", "sense_intronic",
                                  "sense_overlapping", "TEC", "vaultRNA"],
                       "antisense": ["antisense"]}
    for biot in ["miRNA", "snoRNA", "snRNA", "rRNA", "sRNA", "tRNA"]:
        biotypes_group1[biot] = [biot]

    # Initializing the unknown features list
    unknown_cat = set()
    unknown_biot= set()

    # Initializing the genome category counter dict
    cpt_genome = {}
    '''
    if process_counts:
        #### If input files are the categories counts, just load them and continue to recategorization step
        cpt, cpt_genome, samples_names = read_counts_files(options.counts)
    else:
        #### Create genome index if needed and get the sizes of categories
        index_chrom_list = []
        if make_index:
            #### Get the chromosome lengths
            lengths = get_chromosome_lengths(options)
            # Generating the genome index files if the user didn't provide them
            create_genome_index(options.annotation, unstranded_genome_index, stranded_genome_index, cpt_genome, prios,
                                biotypes, lengths)
        else:
            # Retrieving chromosome names saved in index
            index_chrom_list = get_chromosome_names_in_index(genome_index)


        # print '\nChr lengths:', lengths

    if intersect_reads:
        # If the indexes already exist, read them to compute the sizes of the categories in the genome and retrieve the chromosome lengths
        if not make_index:
            print "\n### Reading genome indexes\n",
            sys.stdout.flush()
        lengths = {}
        with open(genome_index, 'r') as genome_index_file:
            for line in genome_index_file:
                if line[0] == "#":
                    lengths[line.split('\t')[0][1:]] = int(line.split('\t')[1])
                else:
                    add_info(cpt_genome, line.rstrip().split('\t')[4:], line.split('\t')[1], line.split('\t')[2],
                             biotype_prios=None, categ_prios=prios)

        #print 'Indexed chromosomes: ' + ', '.join((sorted(index_chrom_list)))
        index_chrom_list.sort(key=alphanum_key)
        print 'Indexed chromosomes: ' + ', '.join((index_chrom_list))

        #### Computing the genome intergenic count: sum of the chr lengths minus sum of the genome annotated intervals
        cpt_genome[('intergenic', 'intergenic')] = sum(lengths.itervalues()) - sum(
            [v for x, v in cpt_genome.iteritems() if x != ('antisense', 'antisense')])
        if not make_index:
            print "Done!"
        # print '\nGenome category counts:'
        # for key,val in cpt_genome.iteritems():
        # print key,"\t",val


        #### Create the Bedgraph files if needed and get the files list

        if not options.bedgraph:
            # Generating the BEDGRAPH files is the user provided BAM file(s) and get the samples labels (this names will be used in the plot legend)
            samples_files, samples_names = create_bedgraph_files(options.input, options.strandness[0])
        else:
            # Just initialize the files list with the bedgraph paths
            # samples_files = [options.input[i] for i in range(0,len(options.input),2)]
            samples_files = [re.sub('.bedgraph$', '', options.input[i]) for i in range(0, len(options.input), 2)]
            # and get the labels
            samples_names = [options.input[i] for i in range(1, len(options.input), 2)]
        #### Retrieving chromosome names saved in index
        #chrom_list = get_chromosome_names_in_index(genome_index)
        #### Processing the BEDGRAPH files: intersecting the bedgraph with the genome index and count the number of aligned positions in each category
        cpt = intersect_bedgraphs_and_index_to_counts_categories(samples_files, samples_names, prios, genome_index,
                                                                 options.strandness[0], biotype_prios=None)

        #### Write the counts on disk
        write_counts_in_files(cpt, cpt_genome)

    if not (intersect_reads or process_counts) or (options.quiet and options.pdf == False):
        quit("\n### End of program")
    print "\n### Generating plots"
    # Updating the biotypes lists (biotypes and 'biotype_group1'): adding the 'unknow biotypes' found in gtf/index
    if unknown_feature == []:  # 'unknown_feature' is define only during the index generation
        # Browse the feature to determine whether some biotypes are 'unknown'
        for sample, counts in cpt.items():
            for (cat, biot) in counts:
                if biot not in biotypes and cat not in unknown_feature:
                    unknown_feature.append(biot)
    for new_biot in unknown_feature:
        biotypes.add(new_biot)
        biotypes_group1["others"].append(new_biot)
    biotypes = sorted(biotypes)
    # move antisense categ to the end of the list
    biotypes.remove('antisense')
    biotypes.append('antisense')
    biotypes_group1 = sorted(biotypes_group1)

    # print '\nCounts for every category/biotype pair: ',cpt

    # Generating plots
    if options.pdf != False:
        if options.pdf == None:
            options.pdf = "categories_plots.pdf"
        pdf = PdfPages(options.pdf)
    else:
        pdf = False

    selected_biotype = None
    if options.biotype_filter:
        options.biotype_filter = options.biotype_filter[0]
        for sample in cpt:
            for feature in cpt[sample]:
                biotype = feature[1]
                if options.biotype_filter.lower() == biotype.lower():
                    selected_biotype = biotype
                    break
        if selected_biotype == None:
            print "\nError: biotype '" + options.biotype_filter + "' not found. Please check the biotype name and that this biotype exists in your sample(s)."
            sys.exit()

    # Print a warning message if the UTRs are not specified as 5' or 3' (they will be ploted as 5'UTR)
    if 'UTR' in [categ[0] for counts in cpt.values() for categ in counts.keys()]:
        print "\nWARNING: (some) 5'UTR/3'UTR are not precisely defined. Consequently, positions annotated as "UTR" will be counted as "5'UTR"\n"
##MB
    #### Make the plot by categories
    #### Recategorizing with the final categories
    final_cats = categs_groups[options.categories_depth - 1]
    final_cat_cpt, final_genome_cpt, filtered_cat_cpt = group_counts_by_categ(cpt, cpt_genome, final_cats, selected_biotype)
    #### Display the distribution of specified categories (or biotypes) in samples on a barplot
    # Remove the 'antisense' category if the library type is 'unstranded'
    for dic in cpt.values():
        if ('antisense', 'antisense') in dic.keys(): break
    else:
        cat_list.remove('antisense')
    make_plot(cat_list, samples_names, final_cat_cpt, final_genome_cpt, pdf, "categories", options.threshold,
              svg=options.svg, png=options.png)
    if selected_biotype:
        make_plot(cat_list, samples_names, filtered_cat_cpt, final_genome_cpt, pdf, "categories", options.threshold,
                  title="Categories distribution for '" + selected_biotype + "' biotype", svg=options.svg, png=options.png)

    #### Make the plot by biotypes
    #### Recategorizing with the final categories
    final_cat_cpt, final_genome_cpt = group_counts_by_biotype(cpt, cpt_genome, biotypes)
    #### Display the distribution of specified categories (or biotypes) in samples on a barplot
    make_plot(biotypes, samples_names, final_cat_cpt, final_genome_cpt, pdf, "biotypes", options.threshold, svg=options.svg,
              png=options.png)

    ##### Recategorizing with the final categories
    # final_cat_cpt,final_genome_cpt = group_counts_by_biotype(cpt,cpt_genome,biotypes_group1)
    ##### Display the distribution of specified categories (or biotypes) in samples on a barplot
    # make_plot(biotypes_group1,samples_names,final_cat_cpt,final_genome_cpt,pdf,"Biotype groups", options.threshold, title="Biotypes distribution in mapped reads \n(biotypes are grouped by 'family')", svg = options.svg, png = options.png)


    if options.pdf:
        pdf.close()
        print "\n### Plots saved in pdf file: %s" % options.pdf

    print "\n### End of program"
    '''
    print "### ALFA ###"
    if not (options.annotation or options.bam or options.bedgraph or options.counts):
        print >> sys.stderr, "\nError: At least one argument among '-a/--annotation', '--bam', '--bedgraph' or " \
                             "'-c/--counts' is required. Please refer to help (-h/--help) and usage cases for more " \
                             "details.\n"
        parser.print_usage()
        sys.exit()

    ## Getting the steps to execute and checking parameters

    # Getting the steps to execute
    generate_indexes = False
    generate_BedGraph = False
    intersect_indexes_BedGraph = False
    generate_plot = False

    # Checking parameters for the step: indexes generation
    if options.annotation:
        generate_indexes = True
    elif options.chr_len:
            unnecessary_param(options.chr_len, "Warning: the parameter '--chr_len' will not be used because the indexes generation step will not be performed.")

    # Checking parameters for the step: BedGraph files generation
    if options.bam:
        if options.bedgraph:
            sys.exit("\nError: parameters '--bam' and '--bedgraph' provided, only one should be.\n### End of program")
        else:
            if not generate_indexes and not options.genome_index:
               sys.exit("\nError: parameter '-g/--genome_index' should be provided.\n### End of program")
            # Checking the input BAM files
            for i in xrange(0, len(options.bam), 2):
                # Check whether the BAM file exists
                try:
                    open(options.bam[i])
                except IOError:
                    sys.exit("\nError: the BAM file " + options.bam[i] + " was not found.\n### End of program")
                except IndexError:
                    sys.exit("\nError: BAM files and associated labels are not correctly provided.\n"
                             "Make sure to follow the expected format: --bam BAM_file1 Label1 [BAM_file2 Label2 ...]."
                             "\n### End of program ###")

                # Check whether the BAM file extension in 'bam'
                if not options.bam[i].endswith(".bam"):
                    sys.exit("\nError: at least one BAM file hasn't a '.bam' extension.\n### End of program ###")
                # Check whether the labels hasn't a "bam" extension
                try:
                    if options.bam[i + 1].endswith(".bam"):
                        sys.exit("\nError: at least one label for a BAM file has a '.bam' extension.\n"
                                 "Make sure to follow the expected format: --bam BAM_file1 Label1 [BAM_file2 Label2 ...]."
                                 "\n### End of program ###")
                except IndexError:
                    sys.exit("\nError: BAM files and associated labels are not correctly provided.\n"
                             "Make sure to follow the expected format: --bam BAM_file1 Label1 [BAM_file2 Label2 ...]."
                             "\n### End of program ###")
                # Get label and replace invalid characters by "_"
                label = "_".join(re.findall(r"[\w\-']+", options.bam[i + 1]))
                # Check whether the BedGraph file(s) and counts file that will be generated already exists
                if options.strandness == "unstranded":
                    existing_file(label + ".bedgraph")
                    existing_file(label + ".unstranded.feature_counts.tsv")
                else:
                    existing_file(label + ".plus.bedgraph")
                    existing_file(label + ".minus.bedgraph")
                    existing_file(label + ".stranded.feature_counts.tsv")
                # Register the sample label and filename
                labels.append(label)
                bams.append(options.bam[i])
            # Set this step + the intersection one as tasks to process
            generate_BedGraph = True
            intersect_indexes_BedGraph = True

    # Checking parameters for the step: indexes and BedGraph files intersection
    if options.bedgraph:
        if not generate_indexes and not options.genome_index:
           sys.exit("\nError: parameter '-g/--genome_index' should be provided.\n### End of program")
        if options.strandness == "unstranded":
            sample_file_nb = 2
        else:
            sample_file_nb = 3
        # Checking the input BedGraph files
        for i in xrange(0, len(options.bedgraph), sample_file_nb):
            # Check whether the BedGraph file(s) exists
            for j in xrange(0,sample_file_nb - 1):
                try:
                    open(options.bedgraph[i + j])
                except IOError:
                    if not (options.bedgraph[i + j].endswith(".bedgraph") or options.bedgraph[i + j].endswith(".bg")):
                        sys.exit("\nError: it looks like BedGraph files and associated labels are not correctly provided.\n"
                             "Make sure to follow the expected format: --bedgraph BedGraph_file1 Label1 [BedGraph_file2 Label2 ...]."
                             "\n### End of program ###")
                    else:
                        sys.exit("\nError: the BedGraph file " + options.bedgraph[i + j] + " was not found.\n### End of program")
                except IndexError:
                    sys.exit("\nError: it looks like BedGraph files and associated labels are not correctly provided.\n"
                         "Make sure to follow the expected format: --bedgraph BedGraph_file1 Label1 [BedGraph_file2 Label2 ...]."
                         "\n### End of program ###")
                # Check whether the BedGraph file extension is "bedgraph"/"bg"
                if not (options.bedgraph[i + j].endswith(".bedgraph") or options.bedgraph[i + j].endswith(".bg")):
                    sys.exit("\nError: the BedGrapg file" + options.bedgraph[i + j] + "  hasn't a '.bedgraph'/'.bg' extension."
                             "\n### End of program ###")
            # Check whether the labels hasn't a "bedgraph"/"bg" extension
            try:
                if options.bedgraph[i  + sample_file_nb - 1].endswith('.bedgraph') or options.bedgraph[i  + sample_file_nb - 1].endswith('.bg'):
                    sys.exit("\nError: the label " + options.bedgraph[i  + sample_file_nb - 1] + " has a '.bedgraph'/'.bg' extension.\n"
                             "Make sure to follow the expected format: "
                             "--bedgraph BedGraph_file1 Label1 [BedGraph_file2 Label2 ...]."
                             "\n### End of program ###")
            except IndexError:
                sys.exit("\nError: it looks like BedGraph files and associated labels are not correctly provided.\n"
                         "Make sure to follow the expected format: --bedgraph BedGraph_file1 Label1 [BedGraph_file2 Label2 ...]."
                         "\n### End of program ###")
            # Register the sample label and filename(s)
            bedgraphs.append(re.sub("(.(plus|minus))?.((bg)|(bedgraph))", "", options.bedgraph[i]))
            label = "_".join(re.findall(r"[\w\-']+", options.bedgraph[i  + sample_file_nb - 1]))
            labels.append(label)
            # Check whether the count file(s) that will be created already exists
            if options.strandness == "unstranded":
                existing_file(label + ".unstranded.feature_counts.tsv")
            else:
                existing_file(label + ".stranded.feature_counts.tsv")
        # Set this step as a task to process
        intersect_indexes_BedGraph = True

    # Checking parameters for the step: plot generation
    if options.counts:
        # Checking unnecessary parameters
        unnecessary_param(options.annotation, "Warning: the parameter '-a/--annotation' will not be used because the counts are already provided.")
        generate_indexes = False
        unnecessary_param(options.genome_index, "Warning: the parameter '-g/--genome_index' will not be used because the counts are already provided.")
        unnecessary_param(options.bam, "Warning: the parameter '--bam' will not be used because the counts are already provided.")
        generate_BedGraph = False
        unnecessary_param(options.bedgraph, "Warning: the parameter '--bedgraph' will not be used because the counts are already provided.")
        intersect_indexes_BedGraph = False
        # Registering the sample labels and filenames
        for sample in options.counts:
            label = os.path.basename(sample)
            label = re.sub('(.(un)?stranded)?.feature_counts.tsv', '', label)
            label = "_".join(re.findall(r"[\w\-']+", label))
            count_files.append(sample)
            labels.append(label)
    else:
        # Generating the genome index filenames
        index_chrom_list = [] # Not a set because we need to sort it according to the chromosome names later
        lengths = {}
        if options.genome_index:
            genome_index_basename = options.genome_index
        elif options.annotation:
            # Otherwise the GTF filename without extension will be the basename
            genome_index_basename = options.annotation.split("/")[-1].split(".gtf")[0]
        stranded_genome_index = genome_index_basename + ".stranded.index"
        unstranded_genome_index = genome_index_basename + ".unstranded.index"
        if options.strandness == "unstranded":
            genome_index = unstranded_genome_index
        else:
            genome_index = stranded_genome_index
    # If no plot is displayed or saved, plot parameters are useless
    if (options.no_display and not (options.pdf or options.png or options.svg)) or not(intersect_indexes_BedGraph or options.counts):
        # unnecessary_param(options.categories_depth, "Warning: the parameter '-d/--categories_depth' will not be used because no plots will be displayed or saved.")
        # # Cannot be tested because options.categories_depth has always a value (default value if option not specified by user)
        unnecessary_param(options.threshold, "Warning: the parameter '-t/--threshold' will not be used because no plots will be displayed or saved.")
        unnecessary_param(options.biotype_filter, "Warning: the parameter '--biotype_filter' will not be used because no plots will be displayed or saved.")
        if options.counts:
            sys.exit("Error: there is nothing to do (counts are provided and no display or plot saving is required")
    else:
        ## [BN] X server check
        # if not X server:
            #options.no_display = True
            #Message
        try:
            x_server = os.environ['DISPLAY']
            generate_plot = True
        except:
            x_server = False
            if options.counts:
                sys.exit("Error: your current configuration does not allow graphical interface ('DISPLAY' variable is not set on your system).\nExiting")
            else:
                print >> sys.stderr, "WARNING: your current configuration does not allow graphical interface ('DISPLAY' variable is not set on your system).\nPlotting step will not be performed."

    ## Executing the step(s)

    # Indexes generation
    if generate_indexes:
        # Checking if the index files already exist
        if os.path.isfile(stranded_genome_index):
            sys.exit("\nError: index files named '" + genome_index_basename +
                     ".stranded.index' already exists but you provided a GTF file. If you want to generate a new index, "
                     "please delete these files or specify an other path.\n### End of program")
        # Running the index generation commands
        print "# Generating the genome index files"
        # Getting chromosomes lengths
        lengths = get_chromosome_lengths()
        # Generating the index files
        generate_genome_index(options.annotation, unstranded_genome_index, stranded_genome_index, lengths)
        # Displaying the list of indexed chromosomes
        index_chrom_list.sort(key=alphanum_key)
        print "Indexed chromosomes: " + ", ".join(index_chrom_list)
    if generate_indexes or not options.counts:
        # Getting index info
        read_index()
        # Computing the genome intergenic count: sum of the chr lengths minus sum of the genome annotated intervals
        cpt_genome[("intergenic", "intergenic")] = sum(lengths.values()) - sum([v for x, v in cpt_genome.iteritems() if x != ("antisense", "antisense")])

    # BedGraph files generation
    if generate_BedGraph:
        print "# Generating the BedGraph files"
        #sample_files, sample_labels = generate_bedgraph_files()
        generate_bedgraph_files(labels, bams)

    # Indexes and BedGraph files intersection
    if intersect_indexes_BedGraph:
        print "# Intersecting index and BedGraph files"
        cpt = intersect_bedgraphs_and_index_to_counts_categories(labels, bedgraphs)
        # Write the counts to an output file
        write_counts_in_files(cpt, cpt_genome)

    ## Plot generation ## MB: all the section still to review
    if generate_plot:
        print "# Generating plots"
        # If input files are the categories counts, the first step is to load them
        if options.counts:
            #cpt, cpt_genome, sample_names = read_counts(options.counts)
            cpt, cpt_genome = read_counts(labels, count_files)
        # Managing the unknown biotypes
        for sample_label, counters in cpt.items():
            for (cat, biot) in counters:
                if biot not in biotypes:
                    unknown_biot.add(biot)
        for biot in unknown_biot:
            biotypes.add(biot)
            biotypes_group1["others"].append(biot)
        biotypes = sorted(biotypes)
        # Moving antisense cat to the end of the list
        biotypes.remove("antisense")
        biotypes.append("antisense")
        biotypes_group1 = sorted(biotypes_group1)
        # Filtering biotypes if necessary
        filtered_biotype = None
        if options.biotype_filter:
            for sample_label in cpt:
                for feature in cpt[sample_label]:
                    biotype = feature[1]
                    if options.biotype_filter.lower() == biotype.lower():
                        selected_biotype = biotype
                        break
            if filtered_biotype:
                print "\nWarning: biotype '" + options.biotype_filter + "' not found. Please check the biotype name and that this biotype exists in your sample(s)."
        # Setting the plots filenames
        if options.pdf: ## MB: Do the same for svg and png??
            pdf = PdfPages(options.pdf)
        else:
            pdf = False
        ## Generate the categories plot
        # Recategorizing within the final categories and plot generation
        final_cats = categs_levels[options.categories_depth - 1]
        parent_categs = parent_categ_groups[options.categories_depth - 1]
        final_cat_cpt, final_genome_cpt, filtered_cat_cpt = group_counts_by_categ(cpt, cpt_genome, final_cats, filtered_biotype)
        # Remove the "antisense" category if the library type is "unstranded" ## MB: if options.strandness == "unstranded": cat_list.remove("antisense")??
        for dic in cpt.values():
            if ("antisense", "antisense") in dic.keys(): break
        else:
            cat_list.remove("antisense")
        #make_plot(labels, cat_list, sample_labels, final_cat_cpt, final_genome_cpt, pdf, "categories", options.threshold, svg=options.svg, png=options.png)
        make_plot(labels, cat_list, final_cat_cpt, final_genome_cpt, pdf, "Categories", options.threshold, svg=options.svg, png=options.png, categ_groups= parent_categs)
        if filtered_biotype:
            #make_plot(labels, cat_list, sample_labels, filtered_cat_cpt, final_genome_cpt, pdf, "categories", options.threshold, title="Categories distribution for '" + filtered_biotype + "' biotype", svg=options.svg, png=options.png)
            make_plot(labels, cat_list, filtered_cat_cpt, final_genome_cpt, pdf, "Categories", options.threshold, title="Categories distribution for '" + filtered_biotype + "' biotype", svg=options.svg, png=options.png, categ_groups= parent_categs)
        ## Generate the biotypes plot
        # Recategorization within the final biotypes and plot generation
        final_cat_cpt, final_genome_cpt = group_counts_by_biotype(cpt, cpt_genome, biotypes)
        #make_plot(biotypes, sample_labels, final_cat_cpt, final_genome_cpt, pdf, "biotypes", options.threshold, svg=options.svg, png=options.png)
        make_plot(labels, biotypes, final_cat_cpt, final_genome_cpt, pdf, "Biotypes", options.threshold, svg=options.svg, png=options.png)


    print "### End of program ###"