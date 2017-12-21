#!/usr/bin/python

import argparse
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile

def exit_and_explain(msg):
    logging.critical(msg)
    sys.exit(msg)

def cleanup_before_exit(tmp_dir):
    if tmp_dir and os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

def get_arg():
    parser = argparse.ArgumentParser()
    parser.add_argument('--project_name', dest='project_name', action='store', nargs=1, metavar='project_name', type=str)
    #Input 1: Annotation File
    parser.add_argument('--index', dest='indexes', action='store', nargs=2, metavar=('stranded_index_filename', 'unstranded_index_filename'), type=str)
    parser.add_argument('--bi_index', dest='bi_indexes', action='store', nargs=1, metavar='built_in_indexes_dir_path', type=str )
    parser.add_argument('--annotation', dest='annotation_file', action='store', nargs=1, metavar='annotation_gtf_file', type=str )
    #Input 2: Mapped Reads
    parser.add_argument('--reads_format', dest='reads_format', action='store', nargs=1, choices=['bam', 'bedgraph'], metavar='reads_format', type=str)
    parser.add_argument('--reads', dest='reads', action='store', nargs='+', metavar=('bam_file1 label1',""), type=str)
    parser.add_argument('--strandness', dest='strandness', action='store', nargs=1, default=['unstranded'], choices=['unstranded', 'forward', 'reverse'], metavar='strandness', type=str)
    #Output files
    parser.add_argument('--output_pdf', dest='output_pdf', action='store', nargs=1, metavar='output_pdf_filename', type=str)
    parser.add_argument('--output_svg', dest='output_svg', action='store', nargs=2, metavar=('categories_svg_filename', 'biotypes_svg_filename'), type=str)
    parser.add_argument('--output_png', dest='output_png', action='store', nargs=2, metavar=('categories_png_filename', 'biotypes_png_filename'), type=str)
    parser.add_argument('--output_count', dest='output_count', action='store', nargs=1, metavar='output_count_filename', type=str)
    parser.add_argument('--output_index', dest='output_indexes', action='store', nargs=2, metavar=('output_stranded_index_filename', 'output_unstranded_index_filename'), type=str)
    #Output Options
    parser.add_argument('--categories_depth', dest='categories_depth', action='store', nargs=1, default=[3], choices=range(1,5), metavar='categories_depth', type=int)
    parser.add_argument('--plot_format', dest='plot_format', action='store', nargs=1, choices=['pdf', 'png', 'svg'], metavar='plot_format', type=str)
    parser.add_argument('--threshold', dest='threshold', action='store', nargs=2, metavar=('yMin', 'yMax'), type=float)
    #Internal variables
    parser.add_argument('--log_report', dest='log_report', action='store', nargs=1, metavar='log_filename', type=str)
    parser.add_argument('--tool_dir', dest='GALAXY_TOOL_DIR', action='store', nargs=1, metavar='galaxy_tool_dir_path', type=str)
    args = parser.parse_args()
    return args

def symlink_user_indexes(stranded_index, unstranded_index, tmp_dir):
    index='index'
    os.symlink(stranded_index, os.path.join(tmp_dir, index + '.stranded.index'))
    os.symlink(unstranded_index, os.path.join(tmp_dir, index + '.unstranded.index'))
    return index

def get_input2_args(reads_list, format, tmp_dir):
    n = len(reads_list)
    if n%2 != 0:
        exit_and_explain('Problem with pairing reads filename and reads label')
    if format == 'bam':
        input2_args = '--bam'
    elif format == 'bedgraph':
        input2_args = '--bedgraph'
    k = 0
    reads_filenames = [''] * (n/2)
    reads_labels = [''] * (n/2)
    for i in range(0, n, 2):
        curr_filename = reads_list[i].split('__fname__')[1]
        # Alfa checks extension so the filename must end either by .bedgraph or by .bam
        # We then create a symlink from file.dat to tmp_dir/annotation_n.<format> to avoid the error message
        reads_filenames[k] = os.path.join(tmp_dir, 'annotation_' + str(k) + '.' + format)
        os.symlink(curr_filename, reads_filenames[k])
        cur_label = reads_list[i+1].split('__label__')[1]
        reads_labels[k] = re.sub(r' ', '_', cur_label)
        if not reads_labels[k]:
            reads_labels[k] = 'sample_%s' % str(k)
        input2_args='%s "%s" "%s"' % (input2_args, reads_filenames[k], reads_labels[k])
        k += 1
    return input2_args, reads_filenames, reads_labels

def redirect_errors(alfa_out, alfa_err):
    # When the option --n is enabled, alfa prints '### End of the program' in stderr even if the process worked-
    # The following lines to avoid the tool from crashing in this case
    if alfa_err and not re.search('### End of program', alfa_err):
        # When alfa prints '### End of program' in stdout, all the messages in stderr are considered
        # as warnings and not as errors.
        if re.search('### End of program', alfa_out):
            logging.warning("The script ALFA.py encountered the following warning:\n\n%s" % alfa_err)
            logging.info("\n******************************************************************\n")
        # True errors make the script exits
        else:
            exit_and_explain("The script ALFA.py encountered the following error:\n\n%s" % alfa_err)

def merge_count_files(reads_labels):
    merged_count_file = open('count_file.txt', 'wb')
    for i in range(0, len(reads_labels)):
        current_count_file = open('%s.feature_counts.tsv' % reads_labels[i], 'r')
        merged_count_file.write('##LABEL: %s\n\n' % reads_labels[i])
        merged_count_file.write(current_count_file.read())
        merged_count_file.write('__________________________________________________________________\n')
        current_count_file.close()
    merged_count_file.close()
    return 'count_file.txt'

def main():
    args = get_arg()

    if not (args.output_pdf or args.output_png or args.output_svg or args.output_indexes or args.output_count):
        exit_and_explain('Error: no output to return\nProcess Aborted\n')
    tmp_dir = tempfile.mkdtemp(prefix='tmp', suffix='')
    logging.basicConfig(level=logging.INFO, filename=args.log_report[0], filemode="a+", format='%(message)s')
    alfa_path = os.path.join(args.GALAXY_TOOL_DIR[0], 'ALFA.py')

    #INPUT1: Annotation File
    if args.indexes:
        # The indexes submitted by the user must exhibit the suffix '.(un)stranded.index' and will be called by alfa by their prefix
        index = symlink_user_indexes(args.indexes[0], args.indexes[1], tmp_dir)
        input1_args = '-g "%s"' % index
    elif args.bi_indexes:
        input1_args = '-g "%s"' % args.bi_indexes[0]
    elif args.annotation_file:
        input1_args = '-a "%s"' % args.annotation_file[0]
    else:
        exit_and_explain('No annotation file submitted !')

    #INPUT 2: Mapped Reads
    if args.reads:
        input2_args, reads_filenames, reads_labels = get_input2_args(args.reads, args.reads_format[0], tmp_dir)
        strandness = '-s %s' % args.strandness[0]
    else:
        exit_and_explain('No reads submitted !')

    ##Output options
    categories_depth = '-d %s' % args.categories_depth[0]
    if not (args.output_pdf or args.output_png or args.output_svg):
        output_args = '--n'
    else:
        plot_suffix = os.path.join(tmp_dir, "ALFA_plot");
        if args.output_pdf:
            output_args = '--pdf ' + plot_suffix + '.pdf'
        if args.output_png:
            output_args = '--png ' + plot_suffix
        if args.output_svg:
            output_args = '--svg ' + plot_suffix
        if args.threshold:
            output_args = '%s -t %.3f %.3f' % (output_args, args.threshold[0], args.threshold[1])

    ##Run alfa
    cmd = 'python %s %s %s %s %s %s' % (alfa_path, input1_args, input2_args, strandness, categories_depth, output_args)
    # Change into the tmp dir because by default, ALFA produces files in the current dir
    curr_dir = os.getcwd()
    os.chdir(tmp_dir)
    print(cmd)
    logging.info("__________________________________________________________________\n")
    logging.info("Alfa execution")
    logging.info("__________________________________________________________________\n")
    logging.info("Command Line:\n%s\n" % cmd)
    logging.info("------------------------------------------------------------------\n")
    alfa_result = subprocess.Popen(args=cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    alfa_out, alfa_err =  alfa_result.communicate()

    ##Handle stdout, warning, errors...
    redirect_errors(alfa_out, alfa_err)

    logging.info("Alfa prompt:\n%s" % alfa_out)

    ##Redirect outputs
    if args.output_pdf:
        shutil.move(plot_suffix + '.pdf', args.output_pdf[0])
    if args.output_png:
        shutil.move(plot_suffix + '.Categories.png', args.output_png[0])
        shutil.move(plot_suffix + '.Biotypes.png', args.output_png[1])
    if args.output_svg:
        shutil.move(plot_suffix + '.Categories.svg', args.output_svg[0])
        shutil.move(plot_suffix + '.Biotypes.svg', args.output_svg[1])
    if args.output_count:
        count_filename = merge_count_files(reads_labels)
        shutil.move(count_filename, args.output_count[0])
    if args.output_indexes:
        if args.annotation_file:
            indexes_regex = re.compile('.*\.index')
            indexes = filter(indexes_regex.search, os.listdir('.'))
            indexes.sort()
            shutil.move(indexes[0], args.output_indexes[0])
            shutil.move(indexes[1], args.output_indexes[1])
        if args.indexes:
            shutil.move(index + '.stranded.index', args.output_indexes[0])
            shutil.move(index + '.unstranded.index', args.output_indexes[1])
        if args.bi_indexes:
            shutil.move(args.bi_indexes[0] + '.stranded.index', args.output_index[0])
            shutil.move(args.bi_indexes[1] + '.unstranded.index', args.output_index[1])

    # Get back to the original dir and cleanup the tmp dir
    os.chdir(curr_dir)
    cleanup_before_exit(tmp_dir)
main()
