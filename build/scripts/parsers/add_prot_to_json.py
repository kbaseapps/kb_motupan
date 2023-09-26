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
import hashlib


# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="add functions to pangenome json")

    parser.add_argument("-i", "--input_clades_file", help="input clades file")
    parser.add_argument("-b", "--base_dir", help="base directory with pangenomes")
    parser.add_argument("-u", "--upa_mapping_file", help="file mapping from genome ID to UPA")
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
    if args.base_dir is None:
        print ("must specify --{}\n".format('base_dir'))
        args_pass = False
    elif not os.path.exists(args.base_dir) or \
         not os.path.isdir(args.base_dir) or \
         not len(os.listdir(args.base_dir)) > 0:
        print ("--{} {} must exist and not be empty\n".format('base_dir', args.base_dir))
        args_pass = False
    if args.upa_mapping_file is None:
        print ("must specify --{}\n".format('upa_mapping_file'))
        args_pass = False
    elif not os.path.exists(args.upa_mapping_file) or \
         not os.path.isfile(args.upa_mapping_file) or \
         not os.path.getsize(args.upa_mapping_file) > 0:
        print ("--{} {} must exist and not be empty\n".format('upa_mapping_file', args.upa_mapping_file))
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
            input_json_file = os.path.join(base_dir, clade, clade+'-mOTUpan-pangenome-fxn.json')
            output_json_file = os.path.join(base_dir, clade, clade+'-mOTUpan-pangenome-fxn-prot.json')
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


# read_genome_ID_to_UPA ()
#
def read_genome_ID_to_UPA (upa_mapping_file, target_genome_IDs):
    genome_IDs_to_UPAs = dict()

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
            #print ("FOUND TARGET UPA for GENOME_ID {}: {}".format(genome_id, upa))  # DEBUG
            genome_IDs_to_UPAs[genome_id] = upa
    f.close()
    
    # DEBUG
    #for upa in sorted(genome_UPAs_to_IDs.keys()):
    #    print ("TARGET UPA {} GENOME_ID {}".format(upa, genome_UPAs_to_IDs[upa]))
    
    return genome_IDs_to_UPAs

    
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


# get_clust_rep_seqs ()
#
def get_clust_rep_seqs (input_json_file):
    clust_rep_seqs = dict()

    clust_rep_seq_file = input_json_file.replace('-mOTUpan-pangenome-fxn.json', '-clust_rep_seq.fasta')

    with open (clust_rep_seq_file, 'r') as clust_rep_seq_h:
        last_clust_id = None
        seq = ''
        for clust_rep_seq_line in clust_rep_seq_h:
            clust_rep_seq_line = clust_rep_seq_line.rstrip()
            if clust_rep_seq_line.startswith('>'):
                clust_id = clust_rep_seq_line.lstrip('>').split()[0]
                if last_clust_id and seq:
                    seq = seq.rstrip('*')
                    clust_rep_seqs[last_clust_id] = seq
                seq = ''
                last_clust_id = clust_id
            else:
                seq += clust_rep_seq_line
        if last_clust_id and seq:
            seq = seq.rstrip('*')
            clust_rep_seqs[last_clust_id] = seq
            seq = ''
            last_clust_id = None
        
    return clust_rep_seqs
    

# get_gene_id_map ()
#
def get_gene_id_map (input_json_file):
    gene_id_map = dict()

    gene_id_map_file = input_json_file.replace('-mOTUpan-pangenome-fxn.json', '.gene_id_map')

    with open (gene_id_map_file, 'r') as gene_id_map_h:
        for gene_id_map_line in gene_id_map_h:
            [genome_based_gene_id, scaffold_based_gene_id] = gene_id_map_line.rstrip().split()
            gene_id_map[genome_based_gene_id] = scaffold_based_gene_id

    return gene_id_map
    

# add_prot_seqs_to_pangenome ()
#
def add_prot_seqs_to_pangenome (pangenome_obj, cluster_rep_seqs, gene_id_map, genome_IDs_to_UPAs):

    print ("adding gene protein translations to pangenome obj ...")

    for cluster_i,cluster in enumerate(pangenome_obj['orthologs']):
        cluster_id = cluster['id']
        
        if cluster_id not in cluster_rep_seqs:
            raise ValueError ("Missing cluster rep seq for cluster id {} in pangenome {}".format(cluster_id, pangeome_obj['name']))

        genome_id = re.sub('_\d+$', '', cluster_id)
        rep_seq_source_gene_id = gene_id_map[cluster_id]
        rep_seq_source_upa = genome_IDs_to_UPAs[genome_id]
        cluster_rep_seq_source = [rep_seq_source_gene_id, rep_seq_source_upa]
        
        pangenome_obj['orthologs'][cluster_i]['protein_translation'] = cluster_rep_seqs[cluster_id]
        pangenome_obj['orthologs'][cluster_i]['protein_translation_source'] = cluster_rep_seq_source
        pangenome_obj['orthologs'][cluster_i]['md5'] = hashlib.md5(cluster_rep_seqs[cluster_id].encode('utf-8')).hexdigest()

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

    # get target genome IDs
    target_genome_IDs = get_target_genome_IDs (input_json_files)

    # read genome IDs to UPAs mapping
    genome_IDs_to_UPAs = read_genome_ID_to_UPA (args.upa_mapping_file, target_genome_IDs)

    # process each clade
    for clade_i,input_json_file in enumerate(input_json_files):

        output_json_file = output_json_files[clade_i]

        # read protein seqs
        cluster_rep_seqs = get_clust_rep_seqs (input_json_file)
        
        # read gene id map
        gene_id_map = get_gene_id_map (input_json_file)
        
        # read pangenome obj from json file
        pangenome_obj = read_pangenome_json (input_json_file)

        # add protein seqs to pangenome obj
        pangenome_obj = add_prot_seqs_to_pangenome (pangenome_obj, cluster_rep_seqs, gene_id_map, genome_IDs_to_UPAs)

        # write updated pangenome json file
        write_pangenome_json_file (output_json_file, pangenome_obj)

    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())
