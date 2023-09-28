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
    parser = argparse.ArgumentParser(description="rm function source field to pangenome json")

    parser.add_argument("-i", "--input_json_file", help="input json file")
    args = parser.parse_args()

    args_pass = True

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit (-1)

    if args.input_json_file is None:
        print ("must specify --{}\n".format('input_json_file'))
        args_pass = False
    elif not os.path.exists(args.input_json_file) or \
         not os.path.isfile(args.input_json_file) or \
         not os.path.getsize(args.input_json_file) > 0:
        print ("--{} {} must exist and not be empty\n".format('input_json_file', args.input_json_infile))
        args_pass = False
        
    if not args_pass:
        parser.print_help()
        sys.exit (-1)
        
    return args


# read_pangenome_json ()
#
def read_pangenome_json (input_json_file):
    print ("reading pangenome json file {} ...".format(input_json_file))
    with open (input_json_file, 'r') as json_h:
        pangenome_obj = json.load (json_h)

    return pangenome_obj


# get_pg_obj_with_cluster_subset ()
#
def get_pg_obj_with_cluster_subset (pangenome_obj, target_cat):
    new_pangenome_obj = dict()
    for field in (list(pangenome_obj.keys())):
        new_pangenome_obj[field] = pangenome_obj[field]
    new_pangenome_obj['id'] += '-'+target_cat
    new_pangenome_obj['name'] += '-'+target_cat
        
    new_clusters = []
    for cluster_i,cluster in enumerate(new_pangenome_obj['orthologs']):
        if cluster['cat'] == target_cat:
            new_clusters.append(cluster)
    new_pangenome_obj['orthologs'] = new_clusters

    return new_pangenome_obj


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

    # read pangenome obj from json file
    pangenome_obj = read_pangenome_json (args.input_json_file)
    
    # get core and write
    output_core_json_file = re.sub('\.json$', '-core.json', args.input_json_file)
    print ("CORE: {}".format(output_core_json_file))
    pangenome_obj_core = get_pg_obj_with_cluster_subset (pangenome_obj, 'core')
    write_pangenome_json_file (output_core_json_file, pangenome_obj_core)
    pangenome_obj_core = dict()
    
    # get accessory and write
    output_acc_json_file = re.sub('\.json$', '-acc.json', args.input_json_file)
    print ("ACC: {}".format(output_acc_json_file))
    pangenome_obj_acc = get_pg_obj_with_cluster_subset (pangenome_obj, 'accessory')
    write_pangenome_json_file (output_acc_json_file, pangenome_obj_acc)

    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())
