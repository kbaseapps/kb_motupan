OB# -*- coding: utf-8 -*-
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
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/kbaseapps/kb_motupan.git"
    GIT_COMMIT_HASH = "282847aabbd7fcf518b3dbabad935a9ed5121575"

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
    
    def now_ISO(self):
        now_timestamp = datetime.now()
        now_secs_from_epoch = (now_timestamp - datetime(1970,1,1)).total_seconds()
        now_timestamp_in_iso = datetime.fromtimestamp(int(now_secs_from_epoch)).strftime('%Y-%m-%d_%T')
        return now_timestamp_in_iso

    def log(self, target, message):
        message = '['+self.now_ISO()+'] '+message
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()

    def check_params (self, params, required_params):
        missing_params = []
        for param in required_params:
            if params.get(param) is None:
                missing_params.append(param)
        if len(missing_params):
            raise ValueError("Missing required param(s):\n" + "\n".join(missing_params))

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
    
    def getUPA_fromInfo (self,obj_info):
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
        return '/'.join([str(obj_info[WSID_I]),
                         str(obj_info[OBJID_I]),
                         str(obj_info[VERSION_I])])

    
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)

        self.workspaceURL = config['workspace-url']
        self.shockURL = config['shock-url']
        self.handleURL = config['handle-service-url']
        self.serviceWizardURL = config['srv-wiz-url']

        self.callbackURL = os.environ.get('SDK_CALLBACK_URL')
        if self.callbackURL == None:
            raise ValueError("SDK_CALLBACK_URL not set in environment")

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
           "file_path", parameter "run_dir" of type "file_path", parameter
           "output_pangenome_json_path" of type "file_path", parameter
           "mmseqs_cluster_mode" of String, parameter "mmseqs_min_seq_id" of
           Double, parameter "mmseqs_min_coverage" of Double, parameter
           "motupan_max_iter" of Long
        :returns: instance of type "run_mmseqs2_and_mOTUpan_files_Output" ->
           structure: parameter "pangenome_json" of type "file_path"
        """
        # ctx is the context object
        # return variables are: returnVal
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
            
            self.log(console, "RUN: "+" ".join(parse_mOTUpan_cmd))
            self.run_subprocess (parse_mOTUpan_cmd, params['run_dir'], console)

            
        # Return
        #
        returnVal = { 'pangenome_json': params['output_pangenome_json_path'] }
        self.log(console, "run_mmseqs2_and_motupan_files DONE")

        #END run_mmseqs2_and_mOTUpan_files

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method run_mmseqs2_and_mOTUpan_files return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
