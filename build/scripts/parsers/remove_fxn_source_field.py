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
    parser.add_argument("-o", "--output_json_file", help="output json file")
    args = parser.parse_args()

    args_pass = True

    if len(sys.argv) < 4:
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
    if args.output_json_file is None:
        print ("must specify --{}\n".format('output_json_file'))
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

# rm_fxn_sources ()
#
def rm_fxn_sources (pangenome_obj):

    for cluster_i,cluster in enumerate(pangenome_obj['orthologs']):
        pangenome_obj['orthologs'][cluster_i]['function_sources'] = []

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

    # read pangenome obj from json file
    pangenome_obj = read_pangenome_json (args.input_json_file)
    
    # remove fxn sources fields
    pangenome_obj = rm_fxn_sources (pangenome_obj)

    # write updated pangenome json file
    write_pangenome_json_file (args.output_json_file, pangenome_obj)

    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())
