#!/usr/bin/python3
'''
Copyright 2022 Dylan Chivian

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
'''

import sys
import argparse
import gzip
import re


# getargs()
#
def getargs():
    default_len = 25
    parser = argparse.ArgumentParser(description="format protein fastas and checkM files for mOTUpan")
    parser.add_argument("-f", "--faalistfile", help="input file with list of protein fasta files")
    parser.add_argument("-c", "--checkminfile", help="input file with checkm scores")
    parser.add_argument("-m", "--motupanfaafile", help="output fasta file for mOTUpan")
    parser.add_argument("-q", "--qualitymotupanfile", help="output quality checkm file for mOTUpan")
    parser.add_argument("-g", "--geneidmappingfile", help="output gene id mapping file")
    args = parser.parse_args()

    if len(sys.argv) < 10:
        parser.print_help()
        sys.exit(-1)
    if not args.faalistfile:
        print ("must specify --faalistfile\n")
        parser.print_help()
        sys.exit (-1)
    if not args.checkminfile:
        print ("must specify --checkminfile\n")
        parser.print_help()
        sys.exit (-1)
    if not args.motupanfaafile:
        print ("must specify --motupanfaafile\n")
        parser.print_help()
        sys.exit (-1)
    if not args.qualitymotupanfile:
        print ("must specify --qualitymotupanfile\n")
        parser.print_help()
        sys.exit (-1)
    if not args.geneidmappingfile:
        print ("must specify --geneidmappingfile\n")
        parser.print_help()
        sys.exit (-1)
        
    return args


# get_fasta_files_from_listfile()
#
def get_fasta_files_from_listfile (listfile):
    print ("reading input file {} ...".format(listfile))
    input_faa_files = dict()
    
    with open(listfile, 'r') as f:
        for line in f:
            filepath = line.strip()
            genome_id = re.sub (r'^.*\/', '', filepath)
            genome_id = re.sub (r'\.faa', '', genome_id)
            genome_id = re.sub (r'\.fasta', '', genome_id)
            genome_id = re.sub (r'_protein', '', genome_id)

            input_faa_files[genome_id] = filepath

    return input_faa_files


# get_genome_ids_from_faas ()
#
def get_genome_ids_from_faas (input_faa_files):
    genome_ids = dict()
    for genome_id in input_faa_files.keys():
        print ("GENOME_ID: "+genome_id)  # DEBUG
        genome_ids[genome_id] = True

    return genome_ids


# get_checkm_scores ()
#
def get_checkm_scores (checkminfile):
    checkm_scores = dict()

    with open(checkminfile, 'r') as f:
        for line in f:
            if line.startswith ('accession'):
                continue
            rec = line.strip()
            (accession, completeness, contamination) = rec.split("\t")
            checkm_scores[accession] = dict()
            checkm_scores[accession]['comp'] = completeness
            checkm_scores[accession]['cont'] = contamination

    return checkm_scores


# write_checkm_file ()
#
def write_checkm_file (qualitymotupanfile, checkm_scores, genome_ids):
    
    with open(qualitymotupanfile, 'w') as f:
        f.write("\t".join(['Bin Id', 'Completeness', 'Contamination'])+"\n")
        for genome_id in sorted(genome_ids.keys()):
            f.write("\t".join([genome_id, checkm_scores[genome_id]['comp'], checkm_scores[genome_id]['cont']])+"\n")
            
    return

                    
# write_faa_file ()
#
def write_faa_file (motupanfaafile, input_faa_files):
    gene_id_mapping = dict()

    with open (motupanfaafile, 'w') as faa_out:
        
        for genome_id in sorted(input_faa_files.keys()):
            gene_cnt = 0
            if input_faa_files[genome_id].lower().endswith('.gz'):
                faa_in = gzip.open(input_faa_files[genome_id], 'rt')
            else:
                faa_in = open(input_faa_files[genome_id], 'r')

            for faa_line in faa_in:
                if faa_line.startswith('>'):
                    gene_cnt += 1
                    old_gene_id = faa_line.split()[0].replace('>','')
                    new_gene_id = genome_id+'_'+str(gene_cnt)
                    gene_id_mapping[new_gene_id] = old_gene_id
                    new_faa_line = faa_line.replace(old_gene_id, new_gene_id)
                    faa_out.write(new_faa_line)
                else:
                    faa_out.write(faa_line)
                    
        faa_in.close()

    return gene_id_mapping


# write_gene_id_mapping_file ()
#
def write_gene_id_mapping_file (geneidmappingfile, gene_id_mapping):

    with open (geneidmappingfile, 'w') as f_out:
        for new_gene_id in sorted(gene_id_mapping.keys()):
            f_out.write("\t".join([new_gene_id, gene_id_mapping[new_gene_id]])+"\n")
            

# main()
#
def main() -> int:
    args = getargs()

    # get faa files
    input_faa_files = get_fasta_files_from_listfile (args.faalistfile)

    # get genome ids
    genome_ids = get_genome_ids_from_faas (input_faa_files)
    
    # get checkm scores
    checkm_scores = get_checkm_scores (args.checkminfile)

    # write new checkm file
    write_checkm_file (args.qualitymotupanfile, checkm_scores, genome_ids)

    # write new faa file
    gene_id_mapping = write_faa_file (args.motupanfaafile, input_faa_files)

    # write gene id mapping file
    write_gene_id_mapping_file (args.geneidmappingfile, gene_id_mapping)
    
    print ("DONE")
    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())
