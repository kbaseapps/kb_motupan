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
import hashlib


# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="parse mOTUpan output into JSON")

    parser.add_argument("-m", "--mOTUpan_infile", help="mOTUpan output file to reformat")
    parser.add_argument("-r", "--reference_map_infile", help="genome name to kbase object reference map file")
    parser.add_argument("-j", "--json_genome_obj_paths_file", help="genome json objs paths mapping file")
    parser.add_argument("-g", "--genefamily_mmseqs_infile", help="mmseqs2 gene family clusters file")
    parser.add_argument("-i", "--id_map_file", help="file with gene id mapping")
    parser.add_argument("-p", "--pangenome_outfile", help="json pangenome out file")
    parser.add_argument("-c", "--completeness_outfile", help="posterior completeness scores calculated by mOTUpan")
    parser.add_argument("-v", "--version_mmseqs2", help="version of MMseqs2 binary")
    parser.add_argument("-a", "--cluster_method_params", help="command line params for clustering")
    parser.add_argument("-b", "--pangenome_method_params", help="command line params for pangenome calc")
    parser.add_argument("-f", "--force_oldfields", help="don't write newer pangenome typedef fields (True/False)")

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
    if args.json_genome_obj_paths_file is not None:
        if not os.path.exists(args.json_genome_obj_paths_file) or \
           not os.path.isfile(args.json_genome_obj_paths_file) or \
           not os.path.getsize(args.json_genome_obj_paths_file) > 0:
            print ("{} {} must exist and not be empty\n".format('json_genome_obj_paths_file', args.json_genome_obj_paths_file))
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

    if args.force_oldfields is None or args.force_oldfields.upper().startswith('F'):
        args.force_oldfields = False
    else:
        args.force_oldfields = True
        
    if not args_pass:
        parser.print_help()
        sys.exit (-1)
        
    return args


# get_genome_name2ref_map ()
#
def get_genome_name2ref_map (reference_map_infile):

    print ("reading genome name to reference map file {} ...".format(reference_map_infile))

    genome_name2ref_map = dict()
    
    if reference_map_infile.lower().endswith('.gz'):
        f = gzip.open(reference_map_infile, 'rt')
    else:
        f = open(reference_map_infile, 'r')

    for line in f:
        line = line.rstrip()

        (genome_name, genome_ref) = line.split("\t")
        genome_name2ref_map[genome_name] = genome_ref
    f.close()

    return genome_name2ref_map


# get_genome_objs ()
#
def get_genome_objs (json_genome_obj_paths_file):

    genome_objs = dict()
    if json_genome_obj_paths_file is not None:
        with open (json_genome_obj_paths_file, 'r') as jgopf:
            for jgopf_line in jgopf:
                jgopf_line = jgopf_line.rstrip()
                [genome_name, json_genome_obj_path] = jgopf_line.split("\t")

                print ("reading genome obj {} from file {} ...".format(genome_name, json_genome_obj_path))
                with open (json_genome_obj_path, 'r') as json_genome_obj_file:
                    genome_objs[genome_name] = json.loads(json_genome_obj_file.read())
        
    return genome_objs


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


# build_pangenome_obj ()
#
def build_pangenome_obj (mOTUpan_infile,
                         genome_name2ref_map,
                         genome_objs,
                         gene2gene_map,
                         cluster_genes,
                         completeness_scores,
                         version_mmseqs2,
                         cluster_method_params_str,
                         pangenome_method_params_str,
                         force_oldfields):
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

    # get genome refs
    genome_refs = []
    genome_ref2name_map = dict()
    if genome_name2ref_map:
        for genome_name in genome_names:
            genome_ref = genome_name2ref_map[genome_name]
            genome_refs.append(genome_ref)
            genome_ref2name_map[genome_ref] = genome_name

    # prep for functions and protein_translation
    gene_names = dict()
    gene_functions = dict()
    protein_translations = dict()
    if not force_oldfields and genome_objs:
        for genome_name in genome_names:
            if genome_name not in genome_objs:
                raise ValueError ("Missing genome {} in genome_objs".format(genome_name))
            genome_obj = genome_objs[genome_name]
            gene_names[genome_name] = dict()
            gene_functions[genome_name] = dict()
            protein_translations[genome_name] = dict()
            for feature in genome_obj['features']:
                fid = feature['id']
                gene_names[genome_name][fid] = []
                gene_functions[genome_name][fid] = []
                protein_translations[genome_name][fid] = ''
                if 'aliases' in feature:
                    for alias in feature['aliases']:
                        [alias_type, alias_val] = alias
                        if alias_type == 'gene':
                            gene_names[genome_name][fid].append(alias_val)
                if 'functions' in feature:
                    gene_functions[genome_name][fid] = feature['functions']
                if 'protein_translation' in feature:
                    protein_translations[genome_name][fid] = feature['protein_translation']
            
    # assign pangenome type and params
    pangenome_type = 'mOTUpan'
    pangenome_method_params = dict()
    if pangenome_method_params_str:
        pg_args = pangenome_method_params_str.split(';')
        for arg in pg_args:
            (pg_key,pg_val) = arg.split('=')
            pangenome_method_params[pg_key] = pg_val
               
    # clustering method, ver, and params
    clustering_method = 'MMseqs2'
    clustering_method_ver = 'bb0a1b3569b9fe115f3bf63e5ba1da234748de23'
    if version_mmseqs2:
        clustering_method_ver = version_mmseqs2
    clustering_method_params = dict()
    if cluster_method_params_str:
        cl_args = cluster_method_params_str.split(';')
        for arg in cl_args:
            (cl_key,cl_val) = arg.split('=')
            clustering_method_params[pg_key] = pg_val

    # assign cluster cats mapping
    cluster_cats = { 'mOTUpan': { 'core': 'core', 'accessory': 'flexible' } }

    
    # get ortholog clusters
    #
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
            this_cluster['id'] = cluster_id
            this_cluster['function'] = ''
            this_cluster['protein_translation'] = ''
            this_cluster['md5'] = ''
            if not force_oldfields:
                this_cluster['genome_occ'] = int(genome_occurences)
                this_cluster['cat'] = cat_acc_core  # either 'accessory' or 'core'
                this_cluster['core_log_likelihood'] = float(log_likelihood_to_be_core)
                this_cluster['mean_copies'] = float(mean_copy_per_genome) * len(cluster_genes[cluster_id])
                this_cluster['function_sources'] = []
                this_cluster['function_logic'] = ''
                this_cluster['protein_translation_source'] = None

            these_genes = []
            this_longest_protein_translation = ''
            this_protein_translation_source = None
            these_gene_names = []
            this_function_logic = 'union'
            these_functions_order = []
            these_functions = dict()
            these_functions_sources = []
            #for gene_id in genes_in_clust.split(';'):
            #    these_genes.append([gene_id, gene2order[gene_id], gene2genome_map[gene_id]])
            for genome_based_gene_id in cluster_genes[cluster_id]:
                scaffold_based_gene_id = re.sub('^.*\.f:', '', gene2gene_map[genome_based_gene_id])
                gene_order = int(re.sub('^.+_', '', genome_based_gene_id))
                genome_name = re.sub('_\d+$', '', genome_based_gene_id)
                genome_ref = genome_name2ref_map[genome_name]
                genome_id = genome_ref
                #if not force_oldfields:
                #    genome_id = genome_name
                these_genes.append([scaffold_based_gene_id,
                                    gene_order,
                                    genome_id])

                if not force_oldfields and genome_objs and genome_name in genome_objs:
                    # gene names
                    if scaffold_based_gene_id in gene_names[genome_name]:
                        for gene_name in gene_names[genome_name][scaffold_based_gene_id]:
                            if gene_name not in these_gene_names:
                                these_gene_names.append(gene_name)
                                
                    # functions
                    if scaffold_based_gene_id in gene_functions[genome_name]:
                        if len(gene_functions[genome_name][scaffold_based_gene_id]) > 0:
                            these_functions_sources.append((scaffold_based_gene_id,genome_ref))
                        for each_function in gene_functions[genome_name][scaffold_based_gene_id]:
                            if each_function not in these_functions:
                                these_functions[each_function] = True
                                these_functions_order.append(each_function)

                    # protein translation
                    if scaffold_based_gene_id in protein_translations[genome_name]:
                        if not this_longest_protein_translation or \
                           len(this_longest_protein_translation) < len(protein_translations[genome_name][scaffold_based_gene_id]):
                            this_longest_protein_translation = protein_translations[genome_name][scaffold_based_gene_id]
                            this_protein_translation_source = (scaffold_based_gene_id,genome_ref)
                            
                    
            this_cluster['orthologs'] = these_genes

            if not force_oldfields:
                this_cluster['gene_name'] = these_gene_names
                this_cluster['function'] = ';'.join(these_functions_order)
                this_cluster['function_sources'] = these_functions_sources
                this_cluster['function_logic'] = this_function_logic
                this_cluster['protein_translation'] = this_longest_protein_translation
                this_cluster['protein_translation_source'] = this_protein_translation_source
                this_cluster['md5'] = hashlib.md5(this_longest_protein_translation.encode('utf-8')).hexdigest()

            orthologs.append(this_cluster)
            
    f.close()

    # build pangenome_obj
    pangenome_obj['name'] = pangenome_name
    pangenome_obj['id'] = pangenome_id
    pangenome_obj['type'] = pangenome_type
    pangenome_obj['orthologs'] = orthologs
    pangenome_obj['genome_refs'] = genome_refs
    if not force_oldfields:
        pangenome_obj['genome_names'] = genome_names
        pangenome_obj['genome_name_to_ref'] = genome_name2ref_map
        pangenome_obj['genome_ref_to_name'] = genome_ref2name_map
        
        pangenome_obj['type_ver'] = type_ver
        pangenome_obj['cluster_cats'] = cluster_cats
        pangenome_obj['clustering_method'] = clustering_method
        pangenome_obj['clustering_method_ver'] = clustering_method_ver
        pangenome_obj['clustering_method_params'] = clustering_method_params
        pangenome_obj['pangenome_method_params'] = pangenome_method_params
        
        pangenome_obj['genome_count'] = int(genome_count)
        pangenome_obj['core_length'] = int(core_length)
        pangenome_obj['mean_est_genome_size'] = float(mean_est_genome_size)
        pangenome_obj['prior_genome_completeness'] = prior_genome_completeness
        pangenome_obj['posterior_genome_completeness'] = posterior_genome_completeness
    
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

    # read genome_name to ref mapping
    genome_name2ref_map = dict()
    if args.reference_map_infile:
        genome_name2ref_map = get_genome_name2ref_map (args.reference_map_infile)

    # read genome objs
    genome_objs = None
    if args.json_genome_obj_paths_file:
        genome_objs = get_genome_objs (args.json_genome_obj_paths_file)
        
    # read gene id to gene id mapping
    gene2gene_map = get_gene2gene_map (args.id_map_file)

    # read cluster gene id members
    cluster_genes = get_cluster_genes (args.genefamily_mmseqs_infile)
    
    # parse out posterior completenss scores and write file
    completeness_scores = get_completeness_scores (args.mOTUpan_infile)
    write_completeness_file (args.completeness_outfile, completeness_scores)

    # parse out clusters and gene ids and write json file
    pangenome_obj = build_pangenome_obj (args.mOTUpan_infile,
                                         genome_name2ref_map,
                                         genome_objs,
                                         gene2gene_map,
                                         cluster_genes,
                                         completeness_scores,
                                         args.version_mmseqs2,
                                         args.cluster_method_params,
                                         args.pangenome_method_params,
                                         args.force_oldfields)
    write_pangenome_json_file (args.pangenome_outfile, pangenome_obj)

    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())
