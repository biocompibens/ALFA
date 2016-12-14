#!/usr/bin/python

import sys
import shutil
import re
import urllib2
import subprocess
import gzip
import os
import tempfile
from optparse import OptionParser
from galaxy.util.json import from_json_string, to_json_string

def get_arg():
    parser = OptionParser()
    parser.add_option("-e", "--ensembl", dest = 'ensembl_info', action = "store", nargs = 2, metavar = ("kingdom", "species_name"), type = "str")
    parser.add_option("-o", "--output", dest='output_filename', action="store", nargs = 1, metavar = 'JSON_FILE')
    parser.add_option("--log", dest='log_filename', action="store", nargs=1, metavar='log_report')
    (options, args) = parser.parse_args()
    return options, args

def cleanup_before_exit(tmp_dir):
    if tmp_dir and os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

def get_page_content(url):
    req = urllib2.Request(url)
    page = urllib2.urlopen(req)
    return page.read()

def download_file(link, local_file_name):
    req = urllib2.Request(link)
    src_file = urllib2.urlopen(req)
    local_file = open(local_file_name, 'wb')
    local_file.write(src_file.read())
    local_file.close()

def uncompress_gz(gz_file_name, uncompressed_file_name):
    print("____________________________________________________________")
    print("*** Uncompressing %s" % gz_file_name)
    uncompressed_file = open(uncompressed_file_name, 'wb')
    with gzip.open(gz_file_name, 'rb') as src_file:
        uncompressed_file.write(src_file.read())
    uncompressed_file.close()
    print("-> Uncompressed !\n")

def add_data_table_entry( data_manager_dict, data_table_entry ):
    data_manager_dict['data_tables'] = data_manager_dict.get( 'data_tables', {} )
    data_manager_dict['data_tables']['alfa_indexes'] = data_manager_dict['data_tables'].get( 'alfa_indexes', data_table_entry )
    return data_manager_dict

def standardize_species_name(species_name):
    # substitute all capital letters, replace every succession of chars that are not letters to one underscore
    standard_species_name = re.sub(r'[)]$', '', species_name)
    standard_species_name = re.sub(r'[ _),-.(=]+ *', '_', standard_species_name)
    return standard_species_name.lower()

def get_ensembl_url_root(kingdom):
    print("____________________________________________________________")
    print("*** Determining Ensembl ftp root url")
    if kingdom == 'vertebrates':
        root = 'ftp://ftp.ensembl.org/pub/current_gtf/'
    else:
        root = 'ftp://ftp.ensemblgenomes.org/pub/%s/current/' % kingdom
    print("-> Determined !\n")
    return root

def test_ensembl_species_exists(kingdom, url, species_name):
    """
    Test if a species exist on the ftp & return the species name with the species_line if so.
    if the species_name matches a single string, then this string will be returned as the species name
    if the species_name matches several strings, then an error is printed with all the possible species to enter for a new run
    """
    print("____________________________________________________________")
    print ("*** Testing whether %s is referenced in Ensembl %s" % (species_name, kingdom))
    list_species_file_name = 'species_Ensembl%s%s.txt' % (kingdom[0].upper(), kingdom[1:])
    if kingdom=='vertebrates':
        download_file(url, list_species_file_name)
    else:
        download_file(url + list_species_file_name, list_species_file_name)

    grep_result = subprocess.Popen(['grep', species_name, list_species_file_name], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    species_lines_matched, grep_error = grep_result.communicate()
    if grep_error != None or species_lines_matched == "":
        msg = 'The species \'%s\' is not referenced on Ensembl (%s)' % (species_name, kingdom)
        sys.exit(msg)

    species_lines = species_lines_matched.split('\n')
    del species_lines[-1]
    nb_lines = len(species_lines)

    if nb_lines == 1:
        if kingdom == 'vertebrates':
            fields = species_lines[0].split(' ')
            columns = fields[-1].split('\r')
            found_species_name = columns[0]
        else:
            columns = species_lines[0].split('\t')
            found_species_name = columns[1]
        if species_name != found_species_name:
            print('-> \'%s\' has been replace with the complete species name \'%s\'' % (species_name, found_species_name))
            return found_species_name, species_lines_matched
        print("-> Referenced !\n")
        return species_name, species_lines_matched
    else:
        list_species = [''] * nb_lines
        for i in range(0, nb_lines):
            if kingdom == 'vertebrates':
                fields = species_lines[i].split(' ')
                columns = fields[-1].split('\r')
                list_species[i] = columns[0]
            else:
                columns = species_lines[i].split('\t')
                list_species[i] = columns[1]
            exact_match = re.search('^%s$' % species_name, list_species[i])
            if exact_match:
                print("-> Referenced !\n")
                return species_name, species_lines[i]
        msg = ("The string \'%s\' has been matched against the list of Ensembl Species but is not a complete species name.\n"
                "Please retry with one of these following species names:\n" % species_name)
        for s in list_species:
            msg = ("%s- %s\n" % (msg, s))
        sys.exit(msg)

def get_ensembl_collection(kingdom, species_line):
    print("*** Extracting the %s_collection of the species" % kingdom)
    collection_regex = re.compile('%s_.+_collection' % kingdom.lower())
    collection_match = re.search(collection_regex, species_line)
    if not collection_match:
        print("-> Skiped: this species is not classified in a Ensembl %s collection\n" % kingdom)
        return None
    print("-> Extracted !\n")
    return collection_match.group(0)

def get_ensembl_gtf_archive_name(url_dir, species_name):
    print("____________________________________________________________")
    print("*** Extracting the gtf archive name of %s" % species_name)
    gtf_archive_regex = re.compile('%s\..*\.[0-9]+\.gtf\.gz' % species_name, flags = re.IGNORECASE)
    dir_content = get_page_content(url_dir)
    gtf_archive_match = re.search(gtf_archive_regex, dir_content)
    if not gtf_archive_match:
        sys.exit('The species is referenced on Ensembl but error of nomenclature led to download failure')
    gtf_archive_name = gtf_archive_match.group(0)
    print("-> Extracted !\n")
    return gtf_archive_name

def get_ensembl_gtf_archive(kingdom, url, species_name, species_line):
    if kingdom != 'vertebrates':
        url = url + 'gtf/' 
        if kingdom == 'bacteria' or kingdom == 'protists' or kingdom == 'fungi':
            collection = get_ensembl_collection(kingdom, species_line)
            if collection != None:
                url = url + "%s/" % collection
    final_url = url + species_name + '/'
    gtf_archive_name = get_ensembl_gtf_archive_name(final_url, species_name)
    print("____________________________________________________________")
    print("*** Download the gtf archive of %s" % species_name)
    download_file(final_url + gtf_archive_name, gtf_archive_name)
    print("-> Downloaded !\n")
    return gtf_archive_name

def generate_alfa_indexes(path_to_alfa, gtf_file_name):
    print("____________________________________________________________")
    print("*** Generating alfa indexes from %s" % gtf_file_name)
    alfa_result = subprocess.Popen(['python', path_to_alfa, '-a',  gtf_file_name], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    alfa_out, alfa_err = alfa_result.communicate()
    if alfa_err != None and not re.search('### End of program', alfa_err):
        msg = 'Generation Failed due an alfa error: %s' % (alfa_err)
        sys.exit(msg)
    print("Alfa prompt:\n%s" % alfa_out)
    print("-> Generated !\n")

def get_data_table_new_entry(gtf_archive_name):
    info_list = gtf_archive_name.split('.')
    species = info_list[0]
    version = info_list[1]
    release = info_list[2]
    value = '%s_%s_%s' % (species, version, release)
    dbkey = value
    name = '%s: %s (release %s)' % (species, version, release)
    prefix = '%s.%s.%s' % (species, version, release)
    entry_dict = { 'species': species, 'version': version, 'release': release, 'value': value, 'dbkey': dbkey, 'name': name, 'prefix': prefix }
    return entry_dict

def main():
    options, args = get_arg()
    tool_dir = args[0]

    path_to_alfa = os.path.join(tool_dir, 'ALFA.py')

    if options.output_filename == None:
        msg = 'No json output file specified'
        sys.exit(msg)
    output_filename = options.output_filename

    # Interestingly the output file to return is not empty initially.
    # it contains a dictionary, with notably the path to the dir where the alfa_indexes
    # are expected to be found
    params = from_json_string(open(output_filename).read())
    target_directory = params['output_data'][0]['extra_files_path']
    os.mkdir(target_directory)

    tmp_dir = tempfile.mkdtemp(prefix='tmp', suffix='')
    os.chdir(tmp_dir)

    data_manager_dict = {}

    if options.ensembl_info:
        kingdom, species_name = options.ensembl_info
        species_name = standardize_species_name(species_name)
        url = get_ensembl_url_root(kingdom)
        species_name, species_line = test_ensembl_species_exists(kingdom, url, species_name)
        gtf_archive_name = get_ensembl_gtf_archive(kingdom, url, species_name, species_line)
        data_table_entry = get_data_table_new_entry(gtf_archive_name)
        gtf_file_name = '%s.gtf' % data_table_entry['prefix']
        uncompress_gz(gtf_archive_name, gtf_file_name)
        generate_alfa_indexes(path_to_alfa, gtf_file_name)
        stranded_index_name = '%s.stranded.index' % data_table_entry['prefix']
        unstranded_index_name = '%s.unstranded.index' % data_table_entry['prefix']
        add_data_table_entry(data_manager_dict, data_table_entry)

    print("____________________________________________________________")
    print("*** General Info")
    print("URL ROOT:\t%s" % url)
    print("SPECIES:\t%s" % data_table_entry['species'])
    print("VERSION:\t%s" % data_table_entry['version'])
    print("RELEASE:\t%s" % data_table_entry['release'])
    print("VALUE:\t%s" % data_table_entry['value'])
    print("DBKEY:\t%s" % data_table_entry['dbkey'])
    print("NAME:\t%s" % data_table_entry['name'])
    print("PREFIX:\t%s" % data_table_entry['prefix'])

    shutil.copyfile(stranded_index_name, os.path.join(target_directory, stranded_index_name))
    shutil.copyfile(unstranded_index_name, os.path.join(target_directory, unstranded_index_name))

    cleanup_before_exit(tmp_dir)

    open(output_filename, 'wb').write(to_json_string(data_manager_dict))
main()
