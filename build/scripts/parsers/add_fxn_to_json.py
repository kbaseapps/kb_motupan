#!/usr/bin/python3
'''
Copyright 2022-2023 Dylan Chivian

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
'''

import sys
import os
import argparse
import gzip
import re
import json


# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="add functions to pangenome json")

    parser.add_argument("-i", "--input_clades_file", help="input clades file")
    parser.add_argument("-b", "--base_dir", help="base directory with pangenomes")
    parser.add_argument("-d", "--domain", help="Bacteria or Archaea")
    parser.add_argument("-p", "--prefered_genomes_file", help="list of prefered species genomes")
    parser.add_argument("-u", "--upa_mapping_file", help="file mapping from genome ID to UPA")
    parser.add_argument("-f", "--function_dir", help="directory where eggNOG functions in mapping format are found")
    args = parser.parse_args()

    args_pass = True

    if len(sys.argv) < 6:
        parser.print_help()
        sys.exit (-1)

    if args.input_clades_file is None:
        print ("must specify --{}\n".format('input_clades_file'))
        args_pass = False
    elif not os.path.exists(args.input_clades_file) or \
         not os.path.isfile(args.input_clades_file) or \
         not os.path.getsize(args.input_clades_file) > 0:
        print ("--{} {} must exist and not be empty\n".format('input_clades_file', args.input_clades_infile))
        args_pass = False
    if args.domain is None:
        print ("must specify --{}\n".format('domain'))
        args_pass = False
    if args.base_dir is None:
        print ("must specify --{}\n".format('base_dir'))
        args_pass = False
    elif not os.path.exists(args.base_dir) or \
         not os.path.isdir(args.base_dir) or \
         not len(os.listdir(args.base_dir)) > 0:
        print ("--{} {} must exist and not be empty\n".format('base_dir', args.base_dir))
        args_pass = False
    if args.prefered_genomes_file is None:
        print ("must specify --{}\n".format('prefered_genomes_file'))
        args_pass = False
    elif not os.path.exists(args.prefered_genomes_file) or \
         not os.path.isfile(args.prefered_genomes_file) or \
         not os.path.getsize(args.prefered_genomes_file) > 0:
        print ("--{} {} must exist and not be empty\n".format('prefered_genomes_file', args.prefered_genomes_file))
        args_pass = False
    if args.upa_mapping_file is None:
        print ("must specify --{}\n".format('upa_mapping_file'))
        args_pass = False
    elif not os.path.exists(args.upa_mapping_file) or \
         not os.path.isfile(args.upa_mapping_file) or \
         not os.path.getsize(args.upa_mapping_file) > 0:
        print ("--{} {} must exist and not be empty\n".format('upa_mapping_file', args.upa_mapping_file))
        args_pass = False
    if args.function_dir is None:
        print ("must specify --{}\n".format('function_dir'))
        args_pass = False
    elif not os.path.exists(args.function_dir) or \
         not os.path.isdir(args.function_dir):
        print ("--{} {} must exist\n".format('function_dir', args.function_dir))
        args_pass = False
        
    if not args_pass:
        parser.print_help()
        sys.exit (-1)
        
    return args


# create_json_paths ()
#
def create_json_paths (base_dir, input_clades_file):
    print ("creating json paths from {} ...".format(input_clades_file))
    input_json_files = []
    output_json_files = []

    with open (input_clades_file, 'r') as clades_h:
        for line in clades_h:
            if line.startswith('#'):
                continue
            #clade = line.rstrip()
            [count, lineage] = line.rstrip().split()
            lineage_list = lineage.split(';')
            clade = lineage_list[-1]
            input_json_file = os.path.join(base_dir, clade, clade+'-mOTUpan-pangenome.json')
            output_json_file = os.path.join(base_dir, clade, clade+'-mOTUpan-pangenome-fxn.json')
            input_json_files.append(input_json_file)
            output_json_files.append(output_json_file)
            
    return (input_json_files, output_json_files)


# read_pangenome_json ()
#
def read_pangenome_json (input_json_file):
    print ("reading pangenome json file {} ...".format(input_json_file))
    with open (input_json_file, 'r') as json_h:
        pangenome_obj = json.load (json_h)

    return pangenome_obj


# read_prefered_genomes ()
#
def read_prefered_genomes (prefered_genomes_file):
    prefered_genomes = dict()

    print ("reading prefered genomes file {} ...".format(prefered_genomes_file))
    if prefered_genomes_file.lower().endswith('.gz'):
        f = gzip.open(prefered_genomes_file, 'rt')
    else:
        f = open(prefered_genomes_file, 'r')
    for line in f:
        if line.startswith('#'):
            continue
        genome_id = line.rstrip()
        prefered_genomes[genome_id] = True
    f.close()

    return prefered_genomes


# read_genome_UPA_to_ID ()
#
def read_genome_UPA_to_ID (upa_mapping_file, target_genome_IDs):
    genome_UPAs_to_IDs = dict()

    print ("reading genome IDs and UPAs file {} ...".format(upa_mapping_file))
    if upa_mapping_file.lower().endswith('.gz'):
        f = gzip.open(upa_mapping_file, 'rt')
    else:
        f = open(upa_mapping_file, 'r')
    for line in f:
        if line.startswith('#'):
            continue
        (genome_id, upa) = line.rstrip().split("\t")
        if genome_id in target_genome_IDs:
            upa = '/'.join(upa.split('/')[0:2])
            #print ("FOUND TARGET UPA for GENOME_ID {}: {}".format(genome_id, upa))  # DEBUG
            genome_UPAs_to_IDs[upa] = genome_id
    f.close()
    
    # DEBUG
    #for upa in sorted(genome_UPAs_to_IDs.keys()):
    #    print ("TARGET UPA {} GENOME_ID {}".format(upa, genome_UPAs_to_IDs[upa]))
    
    return genome_UPAs_to_IDs

    
# get_target_genome_IDs ()
#
def get_target_genome_IDs (input_json_files):
    target_genome_IDs = dict()

    for input_json_file in input_json_files:
        this_pangenome_obj = read_pangenome_json (input_json_file)

        for genome_name in this_pangenome_obj['genome_names']:
            genome_id = genome_name.replace('_protein','')
            #print ("TARGET_GENOME_ID: {}".format(genome_id))  # DEBUG
            target_genome_IDs[genome_id] = True

    return target_genome_IDs


# read_gene_functions ()
#    
def read_gene_functions (function_dir, genome_UPAs_to_IDs, domain):
    gene_names = dict()
    gene_fxns = dict()

    # DEBUG
    #for upa in genome_UPAs_to_IDs.keys():
    #    print ("UPA: {} {}".format(upa, genome_UPAs_to_IDs[upa]))
    
    print ("reading gene functions from dir {} ...".format(function_dir))

    if domain.upper().startswith('A'):
        filename_prefix = 'GTDB_Arc'
    else:
        filename_prefix = 'GTDB_Bac'
    
    UPA_done = dict()
    
    #for f_name in os.listdir (function_dir):
    for subdir, dirs, files in os.walk(function_dir, topdown=True, followlinks=True):
        for f_name in files:

            if not f_name.startswith(filename_prefix):
                continue
            #print ("FILE: {}".format(f_name))

            f_path = os.path.join (subdir,f_name)
            if not os.path.isfile (f_path):
                continue
            print ("\treading gene functions from file {} ...".format(f_name))

            all_done = False
            last_upa = None
        
            if f_path.lower().endswith('.gz'):
                f = gzip.open(f_path_file, 'rt')
            else:
                f = open(f_path, 'r')

            for line in f:
                if line.startswith('#'):
                    continue

                (upa, gene_id, aliases_str, functions_str, inference_str) = line.split("\t")
                upa = '/'.join(upa.split('/')[0:2])
                #print ("UPA {}".format(upa))  # DEBUG

                if last_upa is not None and upa != last_upa and last_upa in genome_UPAs_to_IDs:
                    UPA_done[last_upa] = True
                    print ("UPA DONE {}".format(last_upa))  # DEBUG
                    last_upa = upa
                    all_done = True
                    for this_upa in genome_UPAs_to_IDs.keys():
                        if this_upa not in UPA_done:
                            all_done = False
                            break
                    if all_done:
                        break
                elif last_upa is None or upa != last_upa:
                    last_upa = upa

                if upa not in genome_UPAs_to_IDs:
                    continue

                gene_id = re.sub (r'^gene-', '', gene_id)
                gene_id = re.sub (r'\.CDS$', '', gene_id)

                aliases_list = json.loads(re.sub (r'^"aliases":', '', aliases_str))
                functions_list = json.loads(re.sub (r'^"functions":', '', functions_str))

                these_gene_names = []
                for alias in aliases_list:
                    if alias[0] == 'gene':
                        these_gene_names.append(alias[1])

                genome_id = genome_UPAs_to_IDs[upa]
                    
                if genome_id not in gene_names:
                    gene_names[genome_id] = dict()
                gene_names[genome_id][gene_id] = these_gene_names

                if genome_id not in gene_fxns:
                    gene_fxns[genome_id] = dict()
                gene_fxns[genome_id][gene_id] = functions_list
                    
            f.close()        

            if all_done:
                break
        
            if last_upa is not None and last_upa in genome_UPAs_to_IDs:
                UPA_done[last_upa] = True
                all_done = True
                for this_upa in genome_UPAs_to_IDs.keys():
                    if this_upa not in UPA_done:
                        all_done = False
                        break
            if all_done:
                break

    # DEBUG
    '''
    for genome_id in sorted(gene_names.keys()):
        print ("GENE_NAMES: GENOME_ID: {}".format(genome_id))
        #for gene_id in gene_names[genome_id].keys():
        #    print ("GENE_NAMES: GENE_ID: {}".format(gene_id))
    for genome_id in sorted(gene_fxns.keys()):
        print ("GENE_FXNS: GENOME_ID: {}".format(genome_id))
        #for gene_id in gene_fxns[genome_id].keys():
        #    print ("GENE_FXNS: GENE_ID: {}".format(gene_id))
    # END DEBUG
    '''
            
    return (gene_names, gene_fxns)


# add_gene_fxn_to_pangenome ()
#
def add_gene_fxn_to_pangenome (pangenome_obj, gene_names, gene_fxns, prefered_genomes, genome_UPAs_to_IDs):

    print ("adding gene functions to pangenome obj ...")

    for cluster_i,cluster in enumerate(pangenome_obj['orthologs']):
        prefered_genome_seen = False
        functions_seen = dict()
        these_functions = []
        gene_names_seen = dict()
        these_gene_names = []
        these_function_sources = []
        
        for ortholog in cluster['orthologs']:
            #(gene_id, gene_order, genome_id) = ortholog
            #genome_id = genome_id.replace('_protein','')
            (gene_id, gene_order, genome_upa) = ortholog
            these_function_sources.append([gene_id,genome_upa])
            upa = "/".join(genome_upa.split('/')[0:2])
            genome_id = genome_UPAs_to_IDs[upa]
            #print ("GENOME_ID {}".format(genome_id))  # DEBUG
            if genome_id in prefered_genomes:
                prefered_genome_seen = True
                if genome_id in gene_fxns and \
                   gene_id in gene_fxns[genome_id]:

                    for fxn in gene_fxns[genome_id][gene_id]:
                        if fxn not in functions_seen:
                            functions_seen[fxn] = True
                            these_functions.append(fxn)
                            
                if genome_id in gene_names and \
                   gene_id in gene_names[genome_id]:

                    for gene_name in gene_names[genome_id][gene_id]:
                        if gene_name not in gene_names_seen:
                            gene_names_seen[gene_name] = True
                            these_gene_names.append(gene_name)

        if not prefered_genome_seen:
            pangenome_obj['orthologs'][cluster_i]['gene_name'] = []
            pangenome_obj['orthologs'][cluster_i]['function_sources'] = []
            pangenome_obj['orthologs'][cluster_i]['function'] = "NA"
        else:
            if len(these_gene_names) > 0:
                pangenome_obj['orthologs'][cluster_i]['gene_name'] = these_gene_names
            else:
                pangenome_obj['orthologs'][cluster_i]['gene_name'] = []
            if len(these_functions) > 0:
                pangenome_obj['orthologs'][cluster_i]['function'] = ';'.join(these_functions)
                pangenome_obj['orthologs'][cluster_i]['function_sources'] = these_function_sources
            else:
                pangenome_obj['orthologs'][cluster_i]['gene_name'] = []
                pangenome_obj['orthologs'][cluster_i]['function'] = "hypothetical protein"
                pangenome_obj['orthologs'][cluster_i]['function_sources'] = these_function_sources

    return pangenome_obj

    
# write_pangenome_json_file ()
#
def write_pangenome_json_file (pangenome_outfile, pangenome_obj):
    print ("writing pangenome as json {} ...".format(pangenome_outfile))

    with open(pangenome_outfile, 'w', encoding='utf-8') as f:
        json.dump(pangenome_obj, f, ensure_ascii=False, indent=4)

    return pangenome_outfile


# main()
#
def main() -> int:
    args = getargs()

    # create lists of input and output json pangenome files
    (input_json_files, output_json_files) = create_json_paths (args.base_dir, args.input_clades_file)

    # read prefered genomes
    prefered_genomes = read_prefered_genomes (args.prefered_genomes_file)

    # get target genome IDs
    target_genome_IDs = get_target_genome_IDs (input_json_files)

    # read genome UPA to genome ID mapping
    genome_UPAs_to_IDs = read_genome_UPA_to_ID (args.upa_mapping_file, target_genome_IDs)

    # read annotations and save those for genomes OI
    (gene_names, gene_fxns) = read_gene_functions (args.function_dir, genome_UPAs_to_IDs, args.domain)

    # process each clade
    for clade_i,input_json_file in enumerate(input_json_files):

        output_json_file = output_json_files[clade_i]
        
        # read pangenome obj from json file
        pangenome_obj = read_pangenome_json (input_json_file)
    
        # add annotations to pangenome obj from prefered genomes and other sp reps
        pangenome_obj = add_gene_fxn_to_pangenome (pangenome_obj, gene_names, gene_fxns, prefered_genomes, genome_UPAs_to_IDs)

        # write updated pangenome json file
        write_pangenome_json_file (output_json_file, pangenome_obj)

    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())
