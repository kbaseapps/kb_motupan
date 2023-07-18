# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import sys
import shutil
import re
import json
import subprocess
import traceback
import uuid
import gzip
from datetime import datetime
from pprint import pprint, pformat

# KBase libs
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace as workspaceService
from installed_clients.DataFileUtilClient import DataFileUtil

# KBase modules
from installed_clients.kb_MsuiteClient import kb_Msuite
from installed_clients.kb_phylogenomicsClient import kb_phylogenomics
#END_HEADER


class kb_motupan:
    '''
    Module Name:
    kb_motupan

    Module Description:
    A KBase module: kb_motupan
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.2"
    GIT_URL = "https://github.com/kbaseapps/kb_motupan.git"
    GIT_COMMIT_HASH = "6add5ae5a2c5d1c5fd97654f2047fa7884edd6b6"

    #BEGIN_CLASS_HEADER
    workspaceURL = None
    shockURL = None
    handleURL = None
    callbackURL = None
    scratch = None

    MMSEQS_BINDIR = "/kb/module/mmseqs/bin"
    MMSEQS_BIN = MMSEQS_BINDIR+"/mmseqs"
    MOTUCONVERT_BIN = "/opt/conda3/bin/mOTUconvert.py"
    MOTUPAN_BIN = "/opt/conda3/bin/mOTUpan.py"
    PARSE_MOTUPAN_BIN = "/kb/module/bin/parse_mmseqs_and_mOTUpan.py"

    
    ### now_ISO()
    #
    def now_ISO(self):
        now_timestamp = datetime.now()
        now_secs_from_epoch = (now_timestamp - datetime(1970,1,1)).total_seconds()
        now_timestamp_in_iso = datetime.fromtimestamp(int(now_secs_from_epoch)).strftime('%Y-%m-%d_%T')
        return now_timestamp_in_iso

    
    ### log()
    #
    def log(self, target, message):
        message = '['+self.now_ISO()+'] '+message
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()

        
    ### check_params ()
    #
    def check_params (self, params, required_params):
        missing_params = []
        for param in required_params:
            if params.get(param,'') == '':
                missing_params.append(param)
        if len(missing_params):
            raise ValueError("Missing required param(s):\n" + "\n".join(missing_params))
        return True

    
    ### set_default_params()
    #
    def set_default_params(self, params, default_vals, console):
        for param in sorted(default_vals.keys()):
            if not params.get(param):
                if param in params and params[param] == 0:
                    continue
                else:
                    self.log(console, "Setting param {} to default val {}".format(param, default_vals[param]))
                    params[param] = default_vals[param]
                
        return params
                

    ### run_subprocess ()
    #
    def run_subprocess (self, run_cmd, run_dir, console=None):
        p = subprocess.Popen(run_cmd,
                             cwd=run_dir,
                             stdin=subprocess.DEVNULL,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             shell=False)
                             #shell=True)

        while True:
            line = p.stdout.readline()
            if not line:
                break
            if console is not None:
                self.log(console, str(line).replace('\n', ''))

        p.stdout.close()
        p.wait()
        if console is not None:
            self.log(console, 'return code: ' + str(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running CMD: {}, return code: '.format(" ".join(run_cmd)) + str(p.returncode))

        return p.returncode
    

    ### getUPA_fromInfo ()
    #
    def getUPA_fromInfo (self,obj_info):
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
        return '/'.join([str(obj_info[WSID_I]),
                         str(obj_info[OBJID_I]),
                         str(obj_info[VERSION_I])])

    
    ### validate_and_default_params ()
    #
    def validate_and_default_params (self, params, console):
        required_params = [ 'workspace_name',
                            'input_ref',
                            'output_pangenome_name'
                            ]
        if self.check_params (params, required_params):
            self.log(console, 'All required params met');

        default_vals = {'checkm_version': 'CheckM-1',
                        'mmseqs_cluster_mode': 'easy-cluster',
                        'mmseqs_min_seq_id': 0.0,                        
                        'mmseqs_min_coverage': 0.8,
                        'motupan_max_iter': 1
                        }
        params = self.set_default_params(params, default_vals, console)
        return params


    ### get_genome_objs ()
    #
    def get_genome_objs (self, input_ref, console):
        genome_refs = []
        genome_objs = []
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        top_obj = self.dfuClient.get_objects({'object_refs':[input_ref]})['data'][0]
        top_type = top_obj['info'][TYPE_I].split('-')[0]
        top_data = top_obj['data']
        
        if top_type == 'KBaseTrees.Tree':
            tree_obj = top_data
            for node_id in sorted(tree_obj['ws_refs'].keys()):
                for ref_type in tree_obj['ws_refs'][node_id].keys():
                    genome_refs.extend(tree_obj['ws_refs'][node_id][ref_type])
        elif top_type == 'KBaseSets.GenomeSet':
            genome_refs = [gsi['ref'] for gsi in top_data['items']]
        elif top_type == 'KBaseSearch.GenomeSet':
            genome_refs = [gse['ref'] for gse in top_data['elements'].values()]
        else:
            raise ValueError(f'{top_type} type is not supported')

        # get objs
        self.log(console, "Getting genome objects")
        for genome_ref in genome_refs:
            self.log(console, "Getting genome object for ref {}".format(genome_ref))
            genome_objs.append(self.dfuClient.get_objects({'object_refs':[genome_ref]})['data'][0])
            
        return (genome_refs, genome_objs)


    ### get_genome_qual_scores()
    #
    def get_genome_qual_scores (self, workspace_name, genome_refs, genome_objs, checkm_version, console):
        genome_qual_scores = dict()
        needing_checkm_run = dict()
        genome_names = []

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        # where possible, get completeness from genome obj
        need_to_run_checkm = False
        for genome_i,genome_obj in enumerate(genome_objs):
            genome_name = genome_obj['info'][NAME_I]
            genome_names.append(genome_name)

            score_found = False
            if genome_obj['data'].get('quality_scoress'):
                for qual_score in genome_obj['data']['quality_scores']:
                    if 'method' in qual_score and 'score_interpretation' in qual_score and 'score' in qual_score:
                        if (checkm_version == 'CheckM-1' and qual_score['method'] == 'CheckM') or \
                           (checkm_version == 'CheckM-2' and qual_score['method'] == 'CheckM2'):
                            if qual_score['score_interpretation'] == 'percent_completeness':
                                self.log(console,"found completeness score for {}".format(genome_name))
                                score_found = True
                                if genome_name not in genome_qual_scores:
                                    genome_qual_scores[genome_name] = dict()
                                    genome_qual_scores[genome_name]['contamination'] = 'N/A'
                                genome_qual_scores[genome_name]['completeness'] = qual_score['score']
                            elif qual_score['score_interpretation'] == 'percent_contamination':
                                if genome_name not in genome_qual_scores:
                                    genome_qual_scores[genome_name] = dict()
                                genome_qual_scores[genome_name]['contamination'] = qual_score['score']

            if not score_found:
                need_to_run_checkm = True
                needing_checkm_run[genome_refs[genome_i]] = True

        # run CheckM or CheckM2 for those without scores
        if need_to_run_checkm:

            # DEBUG since can't run CheckM without refdata
            #TEST_MODE = False
            TEST_MODE = True
            if TEST_MODE:
                for genome_i,genome_name in enumerate(genome_names):
                    genome_qual_scores[genome_name] = dict()
                    genome_qual_scores[genome_name]['completeness'] = 90.0
                    genome_qual_scores[genome_name]['contamination'] = 2.0
                return genome_qual_scores
                
            else:
                try:
                    checkM_Client = kb_Msuite(self.callbackURL, token=self.token, service_ver=self.SERVICE_VER)
                except Exception as e:
                    raise ValueError("unable to instantiate checkM_Client. "+str(e))            


            # TODO: this should be run on top obj but currently CheckM doesn't support trees.
            # best solution is to save hidden genomeset of those genomes that need to be scored
            # also, get rid of passing top_obj.  not necessary in this solution
            #
            for genome_i,genome_ref in enumerate(genome_refs):
                if genome_ref not in needing_checkm_run:
                    continue
                
                checkM_params = {'workspace_name': workspace_name,
                                 'input_ref': genome_ref,
                                 'reduced_tree': 1,
                                 'save_output_dir': '0',
                                 'save_plots_dir': '0',
                                 'threads': 4
                                 }                
                
                if checkm_version == 'CheckM-2':
                    raise ValueError ("CheckM2 version not implemented yet")
                else:
                    sub_method = 'CheckM'
                    try:
                        self.log(console, 'RUNNING CheckM for {}'.format(genome_names[genome_i]))
                        this_retVal = checkM_Client.run_checkM_lineage_wf(checkM_params)
                    except Exception as e:
                        raise ValueError ("unable to run "+sub_method+". "+str(e))

                try:
                    this_report_obj = self.wsClient.get_objects2({'objects':[{'ref':this_retVal['report_ref']}]})['data'][0]['data']
                except Exception as e:
                    raise ValueError("unable to fetch "+sub_method+" report: " + this_retVal['report_ref']+". "+str(e))
            
                # retrieve CheckM TSV file
                checkM_outdir = os.path.join(self.scratch, 'checkM.'+str(datetime.now()))
                if not os.path.exists(checkM_outdir):
                    os.makedirs(checkM_outdir)
                checkM_tsv_basefile = 'CheckM_summary_table.tsv'
                checkM_tsv_path = os.path.join(checkM_outdir, checkM_tsv_basefile)
                found_checkM_summary = False
                if len(this_report_obj.get('file_links',[])) > 0:
                    for file_link in this_report_obj['file_links']:
                        if 'name' in file_link and file_link['name'] == checkM_tsv_basefile+'.zip':
                            self.log(console, "CheckM FILE_LINK contents")
                            for key in file_link.keys():
                                self.log(console, "FILE_LINK "+key+": "+str(file_link[key]))
                                
                            download_ret = self.dfuClient.shock_to_file({'handle_id': file_link['handle'],
                                                                         'file_path': checkM_tsv_path+'.zip',
                                                                         'unpack': 'unpack'})
                            # DEBUG
                            #for key in download_ret.keys():
                            #    self.log(console, "DOWNLOAD "+str(key)+": "+str(download_ret[key]))
                            #for inode in os.listdir(checkM_outdir):
                            #    print ("INODE: "+str(inode))
                                
                            found_checkM_summary = True
                            break
                if not found_checkM_summary:
                    raise ValueError ("Failure retrieving CheckM summary TSV file")
                [GENOME_I, LINEAGE_I, GENOME_CNT_I, MARKER_CNT_I, MARKER_SET_I, CNT_0, CNT_1, CNT_2, CNT_3, CNT_4, CNT_5plus, COMPLETENESS_I, CONTAMINATION_I] = range(13)
                self.log(console, "CheckM TSV:")
                with open (checkM_tsv_path, 'r') as checkM_tsv_handle:
                    for checkM_line in checkM_tsv_handle.readlines():
                        checkM_line = checkM_line.rstrip()
                        self.log(console, checkM_line)
                        checkM_info = checkM_line.split("\t")
                        genome_name = checkM_info[GENOME_I]
                        if genome_name == 'Bin Name':
                            continue
                        genome_qual_scores[genome_name] = dict()

                        genome_qual_scores[genome_name]['completeness'] = str(checkM_info[COMPLETENESS_I])  # percent
                        genome_qual_scores[genome_name]['contamination'] = str(checkM_info[CONTAMINATION_I])  # percent
            
        return genome_qual_scores


    ### run_pangenome_circle_plot()
    #
    def run_pangenome_circle_plot (self, pangenome_upa, calling_params, console):
        objects_created = []
        file_links = []
        html_links = []
        
        try:
            phylogenomics_Client = kb_phylogenomics(self.callbackURL, token=self.token, service_ver=self.SERVICE_VER)
        except Exception as e:
            raise ValueError("unable to instantiate phylogenomics_Client. "+str(e))

        if not calling_params.get('pcp_input_genome_ref'):
            self.log(console, "DETERMINING CENTROID GENOME TO USE AS BASE")
            base_genome_ref = self.get_base_genome_ref (pangenome_upa, console)
        else:
            self.log(console, "USING REQUESTED GENOME {} AS BASE".format(calling_params['pcp_input_genome_ref']))
            base_genome_ref = calling_params['pcp_input_genome_ref']
            
        circle_plot_params = {
            'workspace_name': calling_params['workspace_name'],
            'input_genome_ref': base_genome_ref,
            'input_pangenome_ref': pangenome_upa,
            'save_featuresets': calling_params['pcp_save_featuresets'],
            'genome_disp_name_config': calling_params['pcp_genome_disp_name_config']
            }
        if calling_params.get('input_compare_genome_refs'):
            circle_plot_params['input_compare_genome_refs'] = calling_params['input_compare_genome_refs']
        if calling_params.get('input_outgroup_genome_refs'):
            circle_plot_params['input_outgroup_genome_refs'] = calling_params['pcp_input_outgroup_genome_refs']

        this_retVal = phylogenomics_Client.view_pan_circle_plot (circle_plot_params)

        try:
            this_report_obj = self.wsClient.get_objects2({'objects':[{'ref':this_retVal['report_ref']}]})['data'][0]['data']
        except Exception as e:
            raise ValueError("unable to fetch "+sub_method+" report: " + this_retVal['report_ref']+". "+str(e))
        if 'objects_created' in this_report_obj:
            objects_created = this_report_obj['objects_created']
        if 'file_links' in this_report_obj:
            file_links = this_report_obj['file_links']
        if 'html_links' in this_report_obj:
            html_links = this_report_obj['html_links']

        return (objects_created, file_links, html_links)


    ### get_base_genome_ref ()
    #
    def get_base_genome_ref (self, pangenome_upa, console):
        base_genome_ref = None

        pg_obj = self.dfuClient.get_objects({'object_refs':[pangenome_upa]})['data'][0]
        pg_obj_data = pg_obj['data']
        
        centroid_score = dict()

        for cluster in pg_obj_data['orthologs']:
            genome_refs_hit = dict()
            for homolog in cluster['orthologs']:
                (gene_id, gene_order, genome_ref) = homolog
                genome_refs_hit[genome_ref] = True
            this_cluster_genome_refs = sorted(genome_refs_hit.keys())
            num_genomes_in_cluster = len(this_cluster_genome_refs)
            if num_genomes_in_cluster > 1:
                for genome_ref in this_cluster_genome_refs:
                    if genome_ref not in centroid_score:
                        centroid_score[genome_ref] = 0
                    centroid_score[genome_ref] += num_genomes_in_cluster - 1

        high_score = 0
        high_score_genome = None
        for genome_ref in centroid_score.keys():
            if centroid_score[genome_ref] > high_score:
                high_score_genome_ref = genome_ref
                high_score = centroid_score[genome_ref]
                
        base_genome_ref = high_score_genome_ref
        return base_genome_ref
    
            
    ### prepare_motupan_files ()
    #
    def prepare_motupan_files (self, genome_objs, genome_qual_scores, console):
        motupan_input_files = dict()
        stamp = str(int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds() * 1000))

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple


        ### get genome_names and refs
        #
        genome_names = []
        genome_refs = []
        for genome_obj in genome_objs:
            genome_names.append(genome_obj['info'][NAME_I])
            genome_refs.append(self.getUPA_fromInfo(genome_obj['info']))

            
        ### create run directory
        #
        folder = 'mOTUpan_run.'+stamp
        this_run_dir = os.path.join (self.scratch, folder)
        self.log(console,"creating run dir {} ...".format(this_run_dir))

        # create dir and save path to return
        if not os.path.exists (this_run_dir):
            os.makedirs (this_run_dir, mode=0o777, exist_ok=False)
        motupan_input_files['run_dir'] = this_run_dir


        ### create genome name to ref map file
        #
        name2ref_map_buf = []
        name2ref_map_file = os.path.join (this_run_dir, stamp+'-genome_name2ref.map')
        for genome_i,genome_name in enumerate(genome_names):
            name2ref_map_buf.append("\t".join([genome_name, genome_refs[genome_i]]))
            with open (name2ref_map_file, 'w') as name2ref_map_path_handle:
                name2ref_map_path_handle.write("\n".join(name2ref_map_buf)+"\n")
        motupan_input_files['genome_name2ref_path'] = name2ref_map_file

                           
        ### create faa and id_map files
        #
        faa_out_file = os.path.join (this_run_dir, stamp+'.faa')
        id_map_file = os.path.join (this_run_dir, stamp+'.gene_id_map')
        self.log (console,"creating faa file {} ...".format(faa_out_file))

        faa_buf = []
        id_map_buf = []
        for genome_i,genome_obj in enumerate(genome_objs):

            genome_name = genome_names[genome_i]

            # rewrite gene ids to match genome_id as base and store old genome id
            gene_cnt = 0

            #for cds in genome_obj['data']['cdss']:
            for feature in genome_obj['data']['features']:
                if feature.get('protein_translation'):
                    gene_cnt += 1
                    old_gene_id = genome_name+'.f:'+feature['id']
                    new_gene_id = genome_name+'_'+str(gene_cnt)
                    id_map_buf.append("\t".join([new_gene_id,old_gene_id]))
                    faa_buf.append('>'+new_gene_id)
                    faa_buf.append(feature['protein_translation'])
            
        # write files and save path to return
        with open (faa_out_file, 'w') as faa_path_handle:
            faa_path_handle.write("\n".join(faa_buf)+"\n")
        with open (id_map_file, 'w') as id_map_path_handle:
            id_map_path_handle.write("\n".join(id_map_buf)+"\n")
        motupan_input_files['input_faa_path'] = faa_out_file
        motupan_input_files['input_gene_id_map_path'] = id_map_file
            
        
        ### create genome qual file
        #
        checkm_file = os.path.join (this_run_dir, stamp+'.checkm')
        self.log (console,"creating checkm file {} ...".format(checkm_file))

        checkm_buf = []
        header = "\t".join(['Bin Id', 'Completeness', 'Contamination'])
        checkm_buf.append (header)
        for genome_name in genome_names:
            completeness = genome_qual_scores[genome_name]['completeness']
            contamination = genome_qual_scores[genome_name]['contamination']
            checkm_buf.append ("\t".join([genome_name, str(completeness), str(contamination)]))

        # write file and save path to return
        with open (checkm_file, 'w') as checkm_handle:
            checkm_handle.write ("\n".join(checkm_buf)+"\n")
        motupan_input_files['input_qual_path'] = checkm_file


        ### set path for output pangenome json file
        #
        motupan_input_files['output_pangenome_json_path'] = os.path.join (this_run_dir, stamp+'-mOTUpan.json')
        

        return motupan_input_files


    ### save_pangenome_obj ()
    #
    def save_pangenome_obj (self, ctx, input_ref, workspace_name, pangenome_json_file, output_pangenome_name, console):
        pangenome_upa = None

        # set provenance
        self.log(console, "SETTING PROVENANCE")  # DEBUG
        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        provenance[0]['input_ws_objects'] = [input_ref]
        provenance[0]['service'] = 'kb_motupan'
        provenance[0]['method'] = 'run_kb_motupan'        

        # load pg data
        pg_json_str = ''
        with open (pangenome_json_file, 'r') as pg_json_file_h:
            for line in pg_json_file_h:
                pg_json_str += line

        pg_data = json.loads(pg_json_str)
        if 'id' not in pg_data:
            pg_data['id'] = output_pangenome_name
        
        self.log(console, "saving pangenome object {}".format(output_pangenome_name))
        try:
            pg_obj_info = self.wsClient.save_objects(
                { 'workspace': workspace_name,
                  'objects': [{ 'type': 'KBaseGenomes.Pangenome',
                                'name': output_pangenome_name,
                                'data': pg_data,
                                'provenance': provenance
                                }]
                })[0]
        except:
            raise ValueError ("error saving pangenome object")

        pangenome_upa = self.getUPA_fromInfo(pg_obj_info)
        return pangenome_upa


    ### create_motupan_report ()
    #
    def create_motupan_report (self,
                               workspace_name,
                               pangenome_upa,
                               run_dir,
                               pcp_objects_created,
                               pcp_file_links,
                               pcp_html_links,
                               console):

        objects_created = []
        file_links = []
        html_links = []

        report_name = 'kb_motupan_report_' + str(uuid.uuid4())
        
        # archive run
        upload_ret = self.dfuClient.file_to_shock({'file_path': str(run_dir),
                                                   'make_handle': 0,
                                                   'pack': 'zip'})
        output_file_archive = {
            'shock_id': upload_ret['shock_id'],
            'name': 'mOTUpan_run_archive.zip',
            'description': 'mOTUpan run output'
        }
        file_links.append (output_file_archive)

        # put pangenome into objects created
        pg_desc = 'Calculated Pangenome'
        objects_created.append({'ref': pangenome_upa, 'description': pg_desc})

        # add pangenome circle plot output
        objects_created.extend (pcp_objects_created)

        for file_link_item in pcp_file_links:  # file links can't just be extended
            #this_shock_id = file_link_item['URL']
            this_shock_id = re.sub('^.*/', '', file_link_item['URL'])
            new_file_link_item = {'shock_id': this_shock_id,
                                  'name': file_link_item['name'],
                                  'label': file_link_item['label']
            }
            file_links.append(new_file_link_item)
            
        for html_link_item in pcp_html_links:  # html links can't just be extended
            #this_shock_id = html_link_item['URL']
            this_shock_id = re.sub('^.*/', '', html_link_item['URL'])
            new_html_link_item = {'shock_id': this_shock_id,
                                  'name': html_link_item['name'],
                                  'label': html_link_item['label']
            }
            html_links.append(new_html_link_item)            
            

        # create report
        report_params = {
            'direct_html_link_index': 0,
            'html_links': html_links,
            'file_links': file_links,
            'report_object_name': report_name,
            'workspace_name': workspace_name
        }
        if objects_created is not None:
            report_params['objects_created'] = objects_created
            
        report_info = self.reportClient.create_extended_report(report_params)        
        return report_info

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)

        self.token = os.environ['KB_AUTH_TOKEN']
        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.serviceWizardURL = config['srv-wiz-url']

        self.callbackURL = os.environ.get('SDK_CALLBACK_URL')
        if self.callbackURL == None:
            raise ValueError("SDK_CALLBACK_URL not set in environment")

        self.SERVICE_VER = 'release'
        try:
            self.wsClient = workspaceService(self.workspaceURL, token=self.token)
        except:
            raise ValueError ("failed to get wsClient")
        try:
            self.dfuClient = DataFileUtil(self.callbackURL, token=self.token, service_ver=self.SERVICE_VER)
        except:
            raise ValueError ("failed to get dfuClient")
        try:
            self.reportClient = KBaseReport(self.callbackURL, token=self.token, service_ver=self.SERVICE_VER)
        except:
            raise ValueError ("failed to get reportClient")
        
        self.shared_folder = config['scratch']
        self.scratch = os.path.abspath(config['scratch'])
        if self.scratch == None:
            self.scratch = os.path.join('/kb', 'module', 'local_scratch')
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)

        # set i/o dirs
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds() * 1000)
        self.input_dir = os.path.join(self.scratch, 'input.' + str(timestamp))
        self.output_dir = os.path.join(self.scratch, 'output.' + str(timestamp))
        if not os.path.exists(self.input_dir):
            os.makedirs(self.input_dir)
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            
        #END_CONSTRUCTOR
        pass


    def run_mmseqs2_and_mOTUpan_files(self, ctx, params):
        """
        :param params: instance of type
           "run_mmseqs2_and_mOTUpan_files_Params"
           (run_mmseqs2_and_mOTUpan_files() ** **  Method for running mmseqs2
           and mOTUpan from files) -> structure: parameter "input_faa_path"
           of type "file_path", parameter "input_qual_path" of type
           "file_path", parameter "input_gene_id_map_path" of type
           "file_path", parameter "genome_name2ref_path" of type "file_path",
           parameter "run_dir" of type "file_path", parameter
           "output_pangenome_json_path" of type "file_path", parameter
           "mmseqs_cluster_mode" of String, parameter "mmseqs_min_seq_id" of
           Double, parameter "mmseqs_min_coverage" of Double, parameter
           "motupan_max_iter" of Long
        :returns: instance of type "run_mmseqs2_and_mOTUpan_files_Output" ->
           structure: parameter "pangenome_json" of type "file_path"
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_mmseqs2_and_mOTUpan_files
        console = []
        invalid_msgs = []
        tool_name = 'run_mmseqs2_and_motupan_files'
        self.log(console, 'Running ' + tool_name + '_Search with params='
)
        self.log(console, "\n" + pformat(params))
        report = ''
#        report = 'Running '+tool_name+'_Search with params='
#        report += "\n"+pformat(params)

        #### do some basic checks
        #
        required_params =  ['input_faa_path',
                            'input_qual_path',
                            'input_gene_id_map_path',
                            'run_dir',
                            'output_pangenome_json_path',
                            'mmseqs_cluster_mode',
                            'mmseqs_min_seq_id',
                            'mmseqs_min_coverage',
                            'motupan_max_iter'
                            ]
        self.check_params (params, required_params)


        # workflow
        #
        # 1. calculate mmseqs2 clusters
        # 2. format genome clusters for mOTUpan
        # 3. run mOTUpan
        # 4. parse mOTUpan to JSON and add genes in each cluster from mmseqs
        # (5. optionally add functions to pangenome clusters)

        
        # 1. calculate mmseqs2 clusters
        #    Note: subprocess shell must be False.  I think bourne shell messes up mmseqs
        #
        cluster_basename = os.path.basename (params['input_faa_path'])
        cluster_basename = re.sub(r'\.faa', '', cluster_basename)
        cluster_basename = cluster_basename + '-clust'
        mmseqs_cluster_outfile = os.path.join (params['run_dir'], cluster_basename + '_cluster.tsv')

        if int(params.get('force_redo',0)) != 0 or \
           not os.path.isfile (mmseqs_cluster_outfile) or \
           not os.path.getsize (mmseqs_cluster_outfile) > 0:

            cov_mode = "0"
            mmseqs_workdir = 'mmseqs_work'
            mmseqs_cmd = [self.MMSEQS_BIN]
            mmseqs_cmd += [params['mmseqs_cluster_mode']]
            mmseqs_cmd += [params['input_faa_path']]
            mmseqs_cmd += [cluster_basename]
            mmseqs_cmd += [mmseqs_workdir]
            mmseqs_cmd += ['--min-seq-id']
            mmseqs_cmd += [str(params['mmseqs_min_seq_id'])]
            mmseqs_cmd += ['--cov-mode']
            mmseqs_cmd += [str(cov_mode)]
            mmseqs_cmd += ['-c']
            mmseqs_cmd += [str(params['mmseqs_min_coverage'])]
            
            self.log(console, "RUN: "+" ".join(mmseqs_cmd))
            #self.run_subprocess (mmseqs_cmd, params['run_dir'], console)  # too verbose
            self.run_subprocess (mmseqs_cmd, params['run_dir'])

            # clean up mmseqs_workdir because symlinks to nothng mess up archive
            (rundir_path, faa_file) = os.path.split (params['input_faa_path'])
            mmseqs_workdir_path = os.path.join (rundir_path, mmseqs_workdir)
            shutil.rmtree (mmseqs_workdir_path)
            
        
        # 2. format genome clusters for mOTUpan
        #
        motupan_genome_cluster_file = os.path.join (params['run_dir'], cluster_basename+'-motupan_in.json')

        if int(params.get('force_redo',0)) != 0 or \
           not os.path.isfile (motupan_genome_cluster_file) or \
           not os.path.getsize (motupan_genome_cluster_file) > 0:

            mOTUconvert_cmd = [self.MOTUCONVERT_BIN]
            mOTUconvert_cmd += ['--in_type']
            mOTUconvert_cmd += ['mmseqs2']
            mOTUconvert_cmd += ['-o']
            mOTUconvert_cmd += [motupan_genome_cluster_file]
            mOTUconvert_cmd += [mmseqs_cluster_outfile]
            
            self.log(console, "RUN: "+" ".join(mOTUconvert_cmd))
            self.run_subprocess (mOTUconvert_cmd, params['run_dir'], console)

            
        # 3. run mOTUpan
        #
        pangenome_basename = os.path.basename (params['input_faa_path'])
        pangenome_basename = re.sub(r'\.faa', '', pangenome_basename)
        motupan_outfile = os.path.join (params['run_dir'], pangenome_basename+'-pangenome.mOTUpan')
        
        if int(params.get('force_redo',0)) != 0 or \
           not os.path.isfile (motupan_outfile) or \
           not os.path.getsize (motupan_outfile) > 0:

            mOTUpan_cmd = [self.MOTUPAN_BIN]
            mOTUpan_cmd += ['--gene_clusters_file']
            mOTUpan_cmd += [motupan_genome_cluster_file]
            mOTUpan_cmd += ['--checkm']
            mOTUpan_cmd += [params['input_qual_path']]
            mOTUpan_cmd += ['--max_iter']
            mOTUpan_cmd += [str(params['motupan_max_iter'])]
            mOTUpan_cmd += ['--output']
            mOTUpan_cmd += [motupan_outfile]
            
            self.log(console, "RUN: "+" ".join(mOTUpan_cmd))
            self.run_subprocess (mOTUpan_cmd, params['run_dir'], console)

            
        # 4. parse mOTUpan to JSON and add genes in each cluster from mmseqs
        #
        posterior_qual_file = os.path.basename (params['input_qual_path'])
        posterior_qual_file = re.sub(r'\.checkm', '', posterior_qual_file)
        posterior_qual_path = os.path.join (params['run_dir'], posterior_qual_file+'-mOTUpan.qual')
        motupan_json_outfile = motupan_outfile + '.json'
        
        if int(params.get('force_redo',0)) != 0 or \
           not os.path.isfile (motupan_json_outfile) or \
           not os.path.getsize (motupan_json_outfile) > 0 or \
           not os.path.isfile (posterior_qual_path) or \
           not os.path.getsize (posterior_qual_path) > 0:

            parse_mOTUpan_cmd = [self.PARSE_MOTUPAN_BIN]
            parse_mOTUpan_cmd += ['--mOTUpan_infile']
            parse_mOTUpan_cmd += [motupan_outfile]
            parse_mOTUpan_cmd += ['--genefamily_mmseqs_infile']
            parse_mOTUpan_cmd += [mmseqs_cluster_outfile]
            parse_mOTUpan_cmd += ['--id_map_file']
            parse_mOTUpan_cmd += [params['input_gene_id_map_path']]
            parse_mOTUpan_cmd += ['--pangenome_outfile']
            parse_mOTUpan_cmd += [params['output_pangenome_json_path']]
            parse_mOTUpan_cmd += ['--completeness_outfile']
            parse_mOTUpan_cmd += [posterior_qual_path]
            parse_mOTUpan_cmd += ['--force_oldfields']
            parse_mOTUpan_cmd += ['True']  # change to 'False' when pangenome typedef updated
            if params.get('genome_name2ref_path'):
                parse_mOTUpan_cmd += ['--reference_map_infile']
                parse_mOTUpan_cmd += [params['genome_name2ref_path']]
                                    
            self.log(console, "RUN: "+" ".join(parse_mOTUpan_cmd))
            self.run_subprocess (parse_mOTUpan_cmd, params['run_dir'], console)

            
        # Return
        #
        output = { 'pangenome_json': params['output_pangenome_json_path'] }
        self.log(console, "run_mmseqs2_and_motupan_files DONE")

        #END run_mmseqs2_and_mOTUpan_files

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_mmseqs2_and_mOTUpan_files return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def run_kb_motupan(self, ctx, params):
        """
        :param params: instance of type "run_kb_motupan_Params"
           (run_kb_motupan() ** **  Method for running mOTUpan from kbase
           widget) -> structure: parameter "workspace_name" of Long,
           parameter "input_ref" of type "data_obj_ref", parameter
           "output_pangenome_name" of type "data_obj_name", parameter
           "checkm_version" of String, parameter "mmseqs_cluster_mode" of
           String, parameter "mmseqs_min_seq_id" of Double, parameter
           "mmseqs_min_coverage" of Double, parameter "motupan_max_iter" of
           Long, parameter "pcp_input_genome_ref" of type "data_obj_ref",
           parameter "pcp_input_compare_genome_refs" of type "data_obj_ref",
           parameter "pcp_input_outgroup_genome_refs" of type "data_obj_ref",
           parameter "pcp_save_featuresets" of type "bool", parameter
           "pcp_genome_disp_name_config" of String
        :returns: instance of type "ReportResults" (Report results **   
           report_name: The name of the report object in the workspace. **   
           report_ref: The UPA of the report object, e.g. wsid/objid/ver.) ->
           structure: parameter "report_name" of type "data_obj_name",
           parameter "report_ref" of type "data_obj_ref"
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_kb_motupan

        #### STEP 0: init
        console = []
        invalid_msgs = []
        objects_created = []
        file_links = []
        html_dir = os.path.join(self.output_dir, 'html')
        if not os.path.exists(html_dir):
            os.makedirs(html_dir)
        self.log(console, 'Running run_kb_motupan() with params=')
        self.log(console, "\n" + pformat(params))

        #### STEP 1: check and default input params
        self.log(console, "VALIDATING AND DEFAULTING INPUT PARAMS")
        params = self.validate_and_default_params (params, console)

        
        #### STEP 2: get genomes
        self.log(console, "GETTING INPUT GENOME OBJECTS")
        (genome_refs, genome_objs) = self.get_genome_objs (params['input_ref'], console)
        

        #### STEP 3: get completeness scores
        self.log(console, "GETTING COMPLETENESS SCORES")
        genome_qual_scores = self.get_genome_qual_scores (params['workspace_name'],
                                                          genome_refs,
                                                          genome_objs,
                                                          params['checkm_version'],
                                                          console)
        

        ### STEP 4: prepare files
        self.log(console, "PREPARING FILES")
        motupan_input_files = self.prepare_motupan_files (genome_objs, genome_qual_scores, console)
        

        ### STEP 5: run MMseqs2 and mOTUpan on files
        self.log(console, "RUNNING MMSEQS2 and MOTUPAN")
        motupan_files_params = {
            'input_faa_path': motupan_input_files['input_faa_path'],
            'input_qual_path': motupan_input_files['input_qual_path'],
            'genome_name2ref_path': motupan_input_files['genome_name2ref_path'],
            'input_gene_id_map_path': motupan_input_files['input_gene_id_map_path'],
            'run_dir': motupan_input_files['run_dir'],
            'output_pangenome_json_path': motupan_input_files['output_pangenome_json_path'],
            'mmseqs_cluster_mode': params['mmseqs_cluster_mode'],
            'mmseqs_min_seq_id': params['mmseqs_min_seq_id'],
            'mmseqs_min_coverage': params['mmseqs_min_coverage'],
            'motupan_max_iter': params['motupan_max_iter']
        }
        motupan_output_files = self.run_mmseqs2_and_mOTUpan_files (ctx, motupan_files_params)[0]

        
        ### STEP 6: save pangenome object
        self.log(console, "SAVING PANGENOME OUTPUT OBJECT")
        pangenome_upa = self.save_pangenome_obj (ctx,
                                                 params['input_ref'],
                                                 params['workspace_name'],
                                                 motupan_output_files['pangenome_json'],
                                                 params['output_pangenome_name'],
                                                 console)
        

        ### STEP 7: run pangenome circle plot
        self.log(console, "GETTING PANGENOME CIRCLE PLOT")
        circle_plot_limit = 20
        if len (genome_refs) > circle_plot_limit:
            self.log(console, "TOO MANY GENOMES TO PLOT")
        else:
            (pcp_objects_created,
             pcp_file_links,
             pcp_html_links) = self.run_pangenome_circle_plot (pangenome_upa, params, console)

            
        ### STEP 8: make report
        self.log(console, "CREATING REPORT")
        report_info = self.create_motupan_report (params['workspace_name'],
                                                  pangenome_upa,
                                                  motupan_input_files['run_dir'],
                                                  pcp_objects_created,
                                                  pcp_file_links,
                                                  pcp_html_links,
                                                  console)
        
        output = {'report_name': report_info['name'],
                  'report_ref': report_info['ref']
                  }

        self.log(console, "run_kb_motupan() DONE")
        #END run_kb_motupan

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_kb_motupan return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
