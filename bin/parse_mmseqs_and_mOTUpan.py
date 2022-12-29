#!/opt/conda3/bin/python3
'''
Copyright 2022 Dylan Chivian

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
    parser = argparse.ArgumentParser(description="parse mOTUpan output into JSON")

    parser.add_argument("-m", "--mOTUpan_infile", help="mOTUpan output file to reformat")
    parser.add_argument("-g", "--genefamily_mmseqs_infile", help="mmseqs2 gene family clusters file")
    parser.add_argument("-i", "--id_map_file", help="file with gene id mapping")
    parser.add_argument("-p", "--pangenome_outfile", help="json pangenome out file")
    parser.add_argument("-c", "--completeness_outfile", help="posterior completeness scores calculated by mOTUpan")
    args = parser.parse_args()

    args_pass = True

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit (-1)

    if args.mOTUpan_infile is None:
        print ("must specify --{}\n".format('mOTUpan_infile'))
        args_pass = False
    elif not os.path.exists(args.mOTUpan_infile) or \
         not os.path.isfile(args.mOTUpan_infile) or \
         not os.path.getsize(args.mOTUpan_infile) > 0:
        print ("--{} {} must exist and not be empty\n".format('mOTUpan_infile', args.mOTUpan_infile))
        args_pass = False
    if args.genefamily_mmseqs_infile is None:
        print ("must specify --{}\n".format('genefamily_mmseqs_infile'))
        args_pass = False
    elif not os.path.exists(args.genefamily_mmseqs_infile) or \
         not os.path.isfile(args.genefamily_mmseqs_infile) or \
         not os.path.getsize(args.genefamily_mmseqs_infile) > 0:
        print ("--{} {} must exist and not be empty\n".format('genefamily_mmseqs_infile', args.genefamily_mmseqs_infile))
        args_pass = False
    if args.id_map_file is None:
        print ("must specify --{}\n".format('id_map_file'))
        args_pass = False
    elif not os.path.exists(args.id_map_file) or \
         not os.path.isfile(args.id_map_file) or \
         not os.path.getsize(args.id_map_file) > 0:
        print ("--{} {} must exist and not be empty\n".format('id_map_file', args.id_map_file))
        args_pass = False
    if args.pangenome_outfile is None:
        args.pangenome_outfile = args.mOTUpan_infile+'.json'
    if args.completeness_outfile is None:
        args.pangenome_outfile = args.mOTUpan_infile+'.completeness'
        
    if not args_pass:
        parser.print_help()
        sys.exit (-1)
        
    return args


# get_gene2gene_map ()
#
def get_gene2gene_map (id_map_file):

    print ("reading id map file {} ...".format(id_map_file))

    gene2gene_map = dict()
    
    if id_map_file.lower().endswith('.gz'):
        f = gzip.open(id_map_file, 'rt')
    else:
        f = open(id_map_file, 'r')

    for line in f:
        line = line.rstrip()

        (genome_based_gene_id, scaffold_based_gene_id) = line.split("\t")
        gene2gene_map[genome_based_gene_id] = scaffold_based_gene_id
    f.close()

    return gene2gene_map


# get_cluster_genes ()
#
def get_cluster_genes (mmseqs_file):

    print ("reading cluster members file {} ...".format(mmseqs_file))

    cluster_genes = dict()
    
    if mmseqs_file.lower().endswith('.gz'):
        f = gzip.open(mmseqs_file, 'rt')
    else:
        f = open(mmseqs_file, 'r')

    for line in f:
        line = line.rstrip()

        (cluster_id, genome_based_gene_id) = line.split("\t")
        if cluster_id not in cluster_genes:
            cluster_genes[cluster_id] = []
        cluster_genes[cluster_id].append(genome_based_gene_id)
    f.close()

    return cluster_genes


# get_completeness_scores ()
#
def get_completeness_scores (mOTUpan_infile):
    completeness_scores = dict()
    print ("reading mOTUpan file {} for completeness scores ...".format(mOTUpan_infile))
    if mOTUpan_infile.lower().endswith('.gz'):
        f = gzip.open(mOTUpan_infile, 'rt')
    else:
        f = open(mOTUpan_infile, 'r')

    for line in f:
        if line.startswith('#genomes='):
            genome_line = line.rstrip().replace('#genomes=', '')
            for genome_info in genome_line.split(';'):
                [genome_id, prior, posterior] = genome_info.split(':')
                prior_comp = prior.replace('prior_complete=', '')
                posterior_comp = posterior.replace('posterior_complete=', '')
                completeness_scores[genome_id] = posterior_comp
            break
    f.close()

    return completeness_scores


# get_pangenome_obj ()
#
def get_pangenome_obj (mOTUpan_infile, gene2gene_map, cluster_genes, completeness_scores):
    print ("reading mOTUpan file {} for pangenome_obj ...".format(mOTUpan_infile))

    # init structs
    pangenome_obj = dict()
    orthologs = []
    prior_genome_completeness = dict()
    posterior_genome_completeness = dict()

    # get pangenome name
    pangenome_name = re.sub(r'^.*/', '', mOTUpan_infile)
    pangenome_name += '.Pangenome'

    # get genome names
    genome_names = sorted(completeness_scores.keys())

    # assign type
    pangenome_type = 'mOTUpan'
    
    # get ortholog clusters
    if mOTUpan_infile.lower().endswith('.gz'):
        f = gzip.open(mOTUpan_infile, 'rt')
    else:
        f = open(mOTUpan_infile, 'r')

    for line in f:
        line = line.rstrip()
        if line == '':
            continue
        elif line.startswith('#'):

            if line.startswith('#mOTUlizer:mOTUpan:'):
                type_ver = line.replace('#mOTUlizer:mOTUpan:', '')
            elif line.startswith('#run_name='):
                pangenome_id = line.replace('#run_name=', '')
            elif line.startswith('#genome_count='):
                genome_count = line.replace('#genome_count=', '')
            elif line.startswith('#core_length='):
                core_length = line.replace('#core_length=', '')
            elif line.startswith('#mean_est_genome_size='):
                mean_est_genome_size = line.replace('#mean_est_genome_size=', '')
                mean_est_genome_size = mean_est_genome_size.replace(';traits_per_genome', '')
            elif line.startswith('#genomes='):
                genome_line = line.rstrip().replace('#genomes=', '')
                for genome_info in genome_line.split(';'):
                    [genome_id, prior, posterior] = genome_info.split(':')
                    prior_comp = prior.replace('prior_complete=', '')
                    posterior_comp = posterior.replace('posterior_complete=', '')
                    prior_genome_completeness[genome_id] = float(prior_comp)
                    posterior_genome_completeness[genome_id] = float(posterior_comp)
        elif line.startswith('trait_name'):
            continue
        else:
            [cluster_id, cat_acc_core, genome_occurences, log_likelihood_to_be_core, mean_copy_per_genome, genomes_in_clust, genes_in_clust] = line.split("\t")

            # note: genes_in_clust should be 'NA'
            this_cluster = dict()
            this_cluster['function'] = ''
            this_cluster['id'] = cluster_id
            this_cluster['genome_occ'] = int(genome_occurences)
            this_cluster['cat'] = cat_acc_core  # either 'accessory' or 'core'
            this_cluster['core_log_likelihood'] = float(log_likelihood_to_be_core)
            this_cluster['mean_copies'] = float(mean_copy_per_genome) * len(cluster_genes[cluster_id])
            these_genes = []
            #for gene_id in genes_in_clust.split(';'):
            #    these_genes.append([gene_id, gene2order[gene_id], gene2genome_map[gene_id]])
            for genome_based_gene_id in cluster_genes[cluster_id]:
                scaffold_based_gene_id = gene2gene_map[genome_based_gene_id]
                gene_order = re.sub('^GC[AF]_\d+\.\d+_', '', genome_based_gene_id)
                genome_id = re.sub('_\d+$', '', genome_based_gene_id)
                these_genes.append([scaffold_based_gene_id,
                                    gene_order,
                                    genome_id])
            
            this_cluster['orthologs'] = these_genes
            orthologs.append(this_cluster)
            
    f.close()

    # build pangenome_obj
    pangenome_obj['name'] = pangenome_name
    pangenome_obj['id'] = pangenome_id
    pangenome_obj['genome_names'] = genome_names
    pangenome_obj['type'] = 'mOTUpan'
    pangenome_obj['type_ver'] = type_ver
    pangenome_obj['genome_count'] = int(genome_count)
    pangenome_obj['core_length'] = int(core_length)
    pangenome_obj['mean_est_genome_size'] = float(mean_est_genome_size)
    pangenome_obj['prior_genome_completeness'] = prior_genome_completeness
    pangenome_obj['posterior_genome_completeness'] = posterior_genome_completeness
    pangenome_obj['orthologs'] = orthologs
    
    return pangenome_obj
    

# write_completeness_file()
#
def write_completeness_file (completeness_file, completeness_scores):
    print ("writing completeness {} ...".format(completeness_file))
    if completeness_file.lower().endswith('.gz'):
        f = gzip.open(completeness_file, 'wt')
    else:
        f = open(completeness_file, 'w')

    outbuf = []
    outbuf.append("\t".join(['Bin Id','Completeness', 'Contamination']))
    for genome_id in sorted (completeness_scores.keys()):
        outbuf.append("\t".join([genome_id, completeness_scores[genome_id], '-']))
    f.write("\n".join(outbuf)+"\n")
    f.close()

    return completeness_file


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

    # read gene id to gene id mapping
    gene2gene_map = get_gene2gene_map (args.id_map_file)

    # read cluster gene id members
    cluster_genes = get_cluster_genes (args.genefamily_mmseqs_infile)
    
    # parse out posterior completenss scores and write file
    completeness_scores = get_completeness_scores (args.mOTUpan_infile)
    write_completeness_file (args.completeness_outfile, completeness_scores)

    # parse out clusters and gene ids and write json file
    pangenome_obj = get_pangenome_obj (args.mOTUpan_infile, gene2gene_map, cluster_genes, completeness_scores)
    write_pangenome_json_file (args.pangenome_outfile, pangenome_obj)

    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())
