#!/usr/bin/python3
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
import shutil
import json


# getargs()
#
def getargs():
    parser = argparse.ArgumentParser(description="set up for mOTUpan runs with docker kb_motupan module")

    parser.add_argument("-t", "--target_clades_file", help="file with list of target clades")
    parser.add_argument("-s", "--set_name", help="name for this set of targets")
    parser.add_argument("-d", "--docker_base_dir", help="base dir for docker runs - e.g. docker_exec (must be absolute path)")
    parser.add_argument("-i", "--id2ref_genome_map_file", help="genome id to obj ref file")
    parser.add_argument("-g", "--gtdb_metadata_file", help="file with GTDB metadata (checkm scores, lineages, species reps, etc")
    parser.add_argument("-f", "--faa_dir", default="/kbase/ke/db/gtdb_r214/all_faas", help="location of the faa files (def: /kbase/ke/db/gtdb_r214/all_faas)")
    parser.add_argument("-r", "--run_dir", help="where to build the run files - e.g. docker_work (must be relative path)")
    parser.add_argument("-m", "--mount_path", default="/pangenome", help="docker path to run_dir (def: /pangenome)")
    args = parser.parse_args()

    if len(sys.argv) < 5:
        parser.print_help()
        sys.exit(-1)
    args_pass = True
    
    if not args.target_clades_file:
        print ("must specify --{}\n".format('target_clades_file'))
        args_pass = False
    elif not os.path.exists(args.target_clades_file) or \
         not os.path.isfile(args.target_clades_file) or \
         not os.path.getsize(args.target_clades_file) > 0:
        print ("--{} {} must exist and not be empty\n".format('target_clades_file', args.target_clades_file))
        args_pass = False
    if not args.set_name:
        print ("must specify --{}\n".format('set_name'))
        args_pass = False
    if not args.id2ref_genome_map_file:
        print ("must specify --{}\n".format('id2ref_genome_map_file'))
        args_pass = False
    elif not os.path.exists(args.id2ref_genome_map_file) or \
         not os.path.isfile(args.id2ref_genome_map_file) or \
         not os.path.getsize(args.id2ref_genome_map_file) > 0:
        print ("--{} {} must exist and not be empty\n".format('id2ref_genome_map_file', args.id2ref_genome_map_file))
        args_pass = False
    if not args.gtdb_metadata_file:
        print ("must specify --{}\n".format('gtdb_metadata_file'))
        args_pass = False
    elif not os.path.exists(args.gtdb_metadata_file) or \
         not os.path.isfile(args.gtdb_metadata_file) or \
         not os.path.getsize(args.gtdb_metadata_file) > 0:
        print ("--{} {} must exist and not be empty\n".format('gtdb_metadata_file', args.gtdb_metadata_file))
        args_pass = False
    if not args.docker_base_dir:
        print ("must specify --{}\n".format('docker_base_dir'))
        args_pass = False
    elif not os.path.exists(args.docker_base_dir) or \
         not os.path.isdir(args.docker_base_dir) or \
         not len(os.listdir(args.docker_base_dir)) > 0:
        print ("--{} {} must exist and not be empty\n".format('docker_base_dir', args.docker_base_dir))
        args_pass = False
    elif not args.docker_base_dir.startswith('/'):
        print ("--{} must be absolute path, not ".format('docker_base_dir', args.docker_base_dir))
        args_pass = False
    if not args.run_dir:
        print ("must specify --{}\n".format('run_dir'))
        args_pass = False
    elif not os.path.exists(args.run_dir) or \
         not os.path.isdir(args.run_dir):
        print ("--{} {} must exist\n".format('run_dir', args.run_dir))
        args_pass = False
    elif args.run_dir.startswith('/'):
        print ("--{} must be relative path to cwd, not ".format('run_dir', args.run_dir))
        args_pass = False
    if not args.faa_dir:
        print ("must specify --{}\n".format('faa_dir'))
        args_pass = False
    elif not os.path.exists(args.faa_dir) or \
         not os.path.isdir(args.faa_dir) or \
         not len(os.listdir(args.faa_dir)) > 0:
        print ("--{} {} must exist and not be empty\n".format('faa_dir', args.faa_dir))
        args_pass = False

    if not args_pass:
        parser.print_help()
        sys.exit (-1)
        
    return args


# get_target_clades()
#
def get_target_clades (target_clades_file):
    these_target_clades = []
    print ("reading target clades file {} ...".format(target_clades_file))
    if target_clades_file.lower().endswith('.gz'):
        f = gzip.open(clades_file, 'rt')
    else:
        f = open(target_clades_file, 'r')

    this_clade = None
    for line in f:
        line = line.rstrip()
        (this_cnt, this_clade) = line.split("\t")
        these_target_clades.append(this_clade)
    f.close()

    return these_target_clades

# get_genome_name2ref_map ()
#
def get_genome_name2ref_map (id2ref_genome_map_file):
    all_genome_name2ref_map = dict()
    with open (id2ref_genome_map_file, 'r') as id2ref_file:
        for id2ref_line in id2ref_file:
            id2ref_line = id2ref_line.rstrip()
            [genome_id, genome_ref] = id2ref_line.split("\t")

            all_genome_name2ref_map[genome_id] = genome_ref

    return all_genome_name2ref_map

# get_member_genome_ids()
#
def get_member_genome_ids (gtdb_metadata_file, sp_reps_only_for_genus=True):
    print ("reading clade genomes from file {} ...".format(gtdb_metadata_file))
    genus_members = dict()
    species_members = dict()

    GENOME_ID_I            = 0
    #CHECKM_COMPLETENESS_I  = 2
    #CHECKM_CONTAMINATION_I = 3
    GTDB_GENOME_REP_I      = 14
    GTDB_REP_FLAG_I        = 15
    GTDB_TAX_I             = 16
    
    if gtdb_metadata_file.lower().endswith('.gz'):
        f = gzip.open(gtdb_metadata_file, 'rt')
    else:
        f = open(gtdb_metadata_file, 'r')

    for line in f:
        if line.startswith ('accession'):
            continue
        line = line.rstrip()
        metadata = line.split("\t")
        
        genome_id = metadata[GENOME_ID_I]
        genome_id = re.sub(r'^GB_', '', genome_id)
        genome_id = re.sub(r'^RS_', '', genome_id)

        sp_rep_genome_id = metadata[GTDB_GENOME_REP_I]
        sp_rep_genome_id = re.sub(r'^GB_', '', sp_rep_genome_id)
        sp_rep_genome_id = re.sub(r'^RS_', '', sp_rep_genome_id)
        
        sp_rep_flag = metadata[GTDB_REP_FLAG_I]

        species_lineage = metadata[GTDB_TAX_I]
        genus_lineage_list = species_lineage.split(';')
        genus_lineage = ";".join(genus_lineage_list[0:len(genus_lineage_list)-1])

        if species_lineage not in species_members:
            species_members[species_lineage] = []
        species_members[species_lineage].append(genome_id)

        if sp_reps_only_for_genus and sp_rep_flag == 'f':
            continue
        if genus_lineage not in genus_members:
            genus_members[genus_lineage] = []
        genus_members[genus_lineage].append(genome_id)

        
    f.close()

    return (genus_members, species_members)


# get_checkm_scores()
#
def get_checkm_scores (gtdb_metadata_file):
    these_checkm_scores = dict()
    GENOME_ID_I            = 0
    CHECKM_COMPLETENESS_I  = 2
    CHECKM_CONTAMINATION_I = 3
    #GTDB_GENOME_REP_I      = 14
    #GTDB_REP_FLAG_I        = 15
    #GTDB_TAX_I             = 16

    print ("reading CheckM scores from file {} ...".format(gtdb_metadata_file))
    if gtdb_metadata_file.lower().endswith('.gz'):
        f = gzip.open(gtdb_metadata_file, 'rt')
    else:
        f = open(gtdb_metadata_file, 'r')

    for line in f:
        if line.startswith ('accession'):
            continue
        line = line.rstrip()
        metadata = line.split("\t")
        
        genome_id = metadata[GENOME_ID_I]
        genome_id = re.sub(r'^GB_', '', genome_id)
        genome_id = re.sub(r'^RS_', '', genome_id)

        completeness = metadata[CHECKM_COMPLETENESS_I]
        contamination = metadata[CHECKM_CONTAMINATION_I]
        
        these_checkm_scores[genome_id] = (completeness, contamination)
    f.close()

    return these_checkm_scores


# create_run_dir ()
#
def create_run_dir (run_dir, short_clade):
    if short_clade.startswith('g__'):
        folder = 'genus-spreps'
    else:
        folder = 'species'
    this_run_dir = os.path.join (run_dir, folder, short_clade)
    print ("creating run dir {} ...".format(this_run_dir))

    if not os.path.exists (this_run_dir):
        os.makedirs (this_run_dir, mode=0o777, exist_ok=False)
    return this_run_dir


# create_faa_file ()
#
def create_faa_file (faa_dir,
                     this_run_dir,
                     full_clade,
                     short_clade,
                     genome_members):
    faa_out_file = os.path.join (this_run_dir, short_clade+'.faa')
    id_map_file = os.path.join (this_run_dir, short_clade+'.gene_id_map')
    print ("creating faa file {} ...".format(faa_out_file))

    faa_buf = []
    id_map_buf = []

    for genome_id in genome_members[full_clade]:
        db_src = genome_id[0:3]
        f_1 = genome_id[4:7]
        f_2 = genome_id[7:10]
        f_3 = genome_id[10:13]

        faa_in_file = genome_id+'_protein.faa.gz'
        faa_path = os.path.join (faa_dir, db_src, f_1, f_2, f_3, faa_in_file)
        if not os.path.exists (faa_path) or \
           not os.path.getsize (faa_path) > 0:
                
            print ("faa file for {} is missing or empty\n".format(genome_id))
            sys.exit (-2)

        # rewrite gene ids to match genome_id as base and store old genome id
        gene_cnt = 0
        if faa_path.lower().endswith('.gz'):
            faa_in = gzip.open(faa_path, 'rt')
        else:
            faa_in = open(faa_path, 'r')

        for faa_line in faa_in:
            if faa_line.startswith('>'):
                gene_cnt += 1
                old_gene_id = faa_line.split()[0].replace('>','')
                new_gene_id = genome_id+'_'+str(gene_cnt)
                id_map_buf.append("\t".join([new_gene_id,old_gene_id]))
                new_faa_line = faa_line.replace(old_gene_id, new_gene_id)
                faa_buf.append(new_faa_line)
            else:
                faa_buf.append(faa_line)
        faa_in.close()        
            
    # write files
    with open (faa_out_file, 'w') as faa_path_handle:
        faa_path_handle.write("".join(faa_buf))

    with open (id_map_file, 'w') as id_map_path_handle:
        id_map_path_handle.write("\n".join(id_map_buf)+"\n")
                
    return (faa_out_file, id_map_file)


# create_checkm_file ()
#
def create_checkm_file (this_run_dir,
                        target_clade,
                        short_clade,
                        genome_members,
                        all_checkm_scores):
    checkm_file = os.path.join (this_run_dir, short_clade+'.checkm')
    print ("creating checkm file {} ...".format(checkm_file))

    checkm_buf = []
    header = "\t".join(['Bin Id', 'Completeness', 'Contamination'])
    checkm_buf.append (header)
    for genome_id in sorted(genome_members[target_clade]):
        completeness = all_checkm_scores[genome_id][0]
        contamination = all_checkm_scores[genome_id][1]
        checkm_buf.append ("\t".join([genome_id, completeness, contamination]))

    with open (checkm_file, 'w') as checkm_handle:
        checkm_handle.write ("\n".join(checkm_buf)+"\n")
        
    return checkm_file


# create_genome_name2ref_file ()
#
def create_genome_name2ref_file (this_run_dir,
                                 target_clade,
                                 short_clade,
                                 genome_members,
                                 all_genome_name2ref_map):
    genome_name2ref_file = os.path.join (this_run_dir, short_clade+'.genomeid2ref.map')
    print ("creating genome_name2ref file {} ...".format(genome_name2ref_file))

    name2ref_buf = []
    for genome_id in sorted(genome_members[target_clade]):
        genome_ref = all_genome_name2ref_map[genome_id]
        name2ref_buf.append ("\t".join([genome_id, genome_ref]))

    with open (genome_name2ref_file, 'w') as name2ref_handle:
        name2ref_handle.write ("\n".join(name2ref_buf)+"\n")
        
    return genome_name2ref_file


# create_image_dir ()
#
def create_image_dir (docker_base_dir,
                      short_clade):
    image_dir = os.path.join(docker_base_dir, 'images', short_clade)
    print ("creating image dir {} ...".format(image_dir))
    src_sdk_cfg_file = os.path.join(docker_base_dir,'sdk.cfg')
    dst_sdk_cfg_file = os.path.join(image_dir,'sdk.cfg')

    if not os.path.exists (image_dir):
        os.makedirs (image_dir, mode=0o777, exist_ok=False) 

    shutil.copy (src_sdk_cfg_file, dst_sdk_cfg_file)

    return image_dir

    
# create_params_file ()
#
def create_params_file (docker_base_dir,
                        docker_mount_path,
                        short_clade,
                        this_run_dir,
                        faa_path,
                        id_map_path,
                        checkm_path,
                        genome_name2ref_path):

    params_dir = os.path.join (docker_base_dir, 'params')
    if not os.path.exists (params_dir):
        os.makedirs (params_dir, mode=0o777, exist_ok=False) 
    
    params_path = os.path.join (params_dir, short_clade+'.json')
    print ("writing params as json to {} ...".format(params_path))


    output_pg_path = os.path.join(this_run_dir, short_clade+'-mOTUpan-pangenome.json')
    params_obj = {
        'input_faa_path': os.path.join(docker_mount_path, faa_path),
        'input_qual_path': os.path.join(docker_mount_path, checkm_path),
        'genome_name2ref_path': os.path.join(docker_mount_path, genome_name2ref_path),
        'input_gene_id_map_path': os.path.join(docker_mount_path, id_map_path),
        'run_dir': os.path.join(docker_mount_path,this_run_dir),
        'output_pangenome_json_path': os.path.join(docker_mount_path,output_pg_path),
        'mmseqs_cluster_mode': 'easy-cluster',
        'mmseqs_min_seq_id': 0.0,
        'mmseqs_min_coverage': 0.8,
        'motupan_max_iter': 1,
        'force_redo': 0
        #'json_genome_obj_paths_file': motupan_input_files['json_genome_obj_paths_file'],
    }

    with open(params_path, 'w', encoding='utf-8') as f:
        json.dump(params_obj, f, ensure_ascii=False, indent=4)    
    
    return params_path
    

# make_runner_cmd ()
#
def make_runner_cmd (clade_i,
                     short_clade,
                     docker_mount_path,
                     image_dir,
                     params_file):
    cmd = []
    full_mount_path = os.getcwd()
    full_image_dir_path = os.path.join(os.getcwd(), image_dir)
    full_params_file_path = os.path.join(os.getcwd(), params_file)
    method = 'kb_motupan.run_mmseqs2_and_mOTUpan_files'
    
    cmd.append('kb-sdk run')
    cmd.append('-t beta')
    cmd.append('--input '+full_params_file_path)
    cmd.append('--mount-points '+full_mount_path+':'+docker_mount_path)
    cmd.append('--sdk-home '+full_image_dir_path)
    cmd.append(method)

    report_cmd = 'echo; echo RUNNING CLADE NUMBER '+str(clade_i+1)+'; ' + 'echo '+" ".join(cmd)+'; '
    
    return report_cmd+" ".join(cmd)
    

# create_run_script_file ()
#
def create_run_script_file (docker_base_dir,
                            set_name,
                            super_runner_buf):
    scripts_dir = os.path.join(docker_base_dir, 'scripts')
    if not os.path.exists (scripts_dir):
        os.makedirs (scripts_dir, mode=0o777, exist_ok=False) 
    run_script_file = os.path.join(scripts_dir, set_name+'.sh')

    with open (run_script_file, 'w') as script_h:
        script_h.write("\n".join(super_runner_buf)+"\n")

    return run_script_file


# main()
#
def main() -> int:
    args = getargs()

    # gather needed info
    target_clades = get_target_clades (args.target_clades_file)
    all_checkm_scores = get_checkm_scores (args.gtdb_metadata_file)
    (all_genus_members, all_species_members) = get_member_genome_ids (args.gtdb_metadata_file)
    all_genome_name2ref_map = get_genome_name2ref_map (args.id2ref_genome_map_file)
    

    # DEBUG
    """
    for target_clade in target_clades:
        print ("TARGET_CLADE: '{}'".format(target_clade))

    #for genome_id in sorted (all_checkm_scores.keys()):
    #    print ("{} CHECKM_SCORES: {} {}".format(genome_id, all_checkm_scores[genome_id][0], all_checkm_scores[genome_id][1]))

    #for genus in sorted (all_genus_members.keys()):
    #    print ("{} GENOME MEMBERS: {}".format(genus, ", ".join(all_genus_members[genus])))
        
    for species in sorted (all_species_members.keys()):
        print ("{} GENOME MEMBERS: {}".format(species, ", ".join(all_species_members[species])))
        
    return 0
    """
    
    
    # init runner buf
    super_runner_buf = ['#!/bin/sh']
    
    # add each target clade
    for clade_i,target_clade in enumerate(target_clades):
        short_clade = target_clade.split(';')[-1]
        short_clade = short_clade.replace(' ', '_')

        print ("processing {}".format(short_clade))
        
        # determine mode
        if short_clade.startswith ('g__'):
            genome_members = all_genus_members
        else:
            genome_members = all_species_members
        
        # create run_dir
        this_run_dir = create_run_dir (args.run_dir, short_clade)
        
        # create faa file and id map file
        (faa_file, id_map_file) = create_faa_file (args.faa_dir,
                                                   this_run_dir,
                                                   target_clade,
                                                   short_clade,
                                                   genome_members)

        # create checkm file
        checkm_file = create_checkm_file (this_run_dir,
                                          target_clade,
                                          short_clade,
                                          genome_members,
                                          all_checkm_scores)

        # create genome_name2ref file
        genome_name2ref_file = create_genome_name2ref_file (this_run_dir,
                                                            target_clade,
                                                            short_clade,
                                                            genome_members,
                                                            all_genome_name2ref_map)

        # create image dir
        image_dir = create_image_dir (args.docker_base_dir,
                                      short_clade)

        # create params
        params_file = create_params_file (args.docker_base_dir,
                                          args.mount_path,
                                          short_clade,
                                          this_run_dir,
                                          faa_file,
                                          id_map_file,
                                          checkm_file,
                                          genome_name2ref_file)

        # write run mOTUpan script
        runner_cmd = make_runner_cmd (clade_i,
                                      short_clade,
                                      args.mount_path,
                                      image_dir,
                                      params_file)
        super_runner_buf.append (runner_cmd)


    # create run script
    run_script_file = create_run_script_file (args.docker_base_dir,
                                              args.set_name,
                                              super_runner_buf)

    return 0


# exec()
#
if __name__ == '__main__':
    sys.exit(main())
