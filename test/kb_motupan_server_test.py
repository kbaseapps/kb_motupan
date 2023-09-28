# -*- coding: utf-8 -*-
import os
import time
import unittest
import shutil
import json
from shutil import copyfile
from pathlib import Path
from pprint import pprint
from configparser import ConfigParser

from kb_motupan.kb_motupanImpl import kb_motupan
from kb_motupan.kb_motupanServer import MethodContext
from installed_clients.authclient import KBaseAuth as _KBaseAuth
from installed_clients.GenomeFileUtilClient import GenomeFileUtil

from installed_clients.WorkspaceClient import Workspace


class kb_motupanTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_motupan'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(cls.token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': cls.token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_motupan',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        cls.wsClient = Workspace(cls.wsURL)
        cls.gfu = GenomeFileUtil(cls.callback_url, token=cls.token)
        cls.serviceImpl = kb_motupan(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_mOTUpan_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

        cls.prepare_data()
        

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')
        if hasattr(cls, 'ws2Name'):
            cls.wsClient.delete_workspace({'workspace': cls.ws2Name})
            print('Test workspace2 was deleted')
        if hasattr(cls, 'shock_ids'):
            for shock_id in cls.shock_ids:
                print('Deleting SHOCK node: ' + str(shock_id))
                cls.delete_shock_node(shock_id)

    @classmethod
    def delete_shock_node(cls, node_id):
        header = {'Authorization': 'Oauth {0}'.format(cls.token)}
        requests.delete(cls.shockURL + '/node/' + node_id, headers=header,
                        allow_redirects=True)
        print('Deleted shock node ' + node_id)

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_kb_SetUtilities_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})
        self.__class__.wsName = wsName
        return wsName

    def getWs2Name(self):
        if hasattr(self.__class__, 'ws2Name'):
            return self.__class__.ws2Name
        suffix = int(time.time() * 1000)
        ws2Name = "test_kb_SetUtilities_" + str(suffix)+'-2'
        ret = self.getWsClient().create_workspace({'workspace': ws2Name})
        self.__class__.ws2Name = ws2Name
        return ws2Name

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    @classmethod
    def ref_from_info(cls, obj_info):
        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple
        return "/".join([str(obj_info[WSID_I]), str(obj_info[OBJID_I]), str(obj_info[VERSION_I])])

    def isUpa (self, candidate_upa):
        legit_upa = True
        for upa_element in candidate_upa.split('/'):
            if not upa_element.isdigit():
                print ("Error: not UPA element: "+upa_element)
                legit_upa = False
        return legit_upa
        
    @classmethod
    def prepare_data(cls):
        tempdir = Path(os.path.join(cls.scratch, 'tempstuff'))
        tempdir.mkdir(parents=True, exist_ok=True)

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        # prepare genome set and tree
        cls.genomes = []
        gs_genome_elements = {}
        tree_ws_refs = dict()
        tree_default_node_labels = dict()
        tree_leaf_list = []
        tree_newick = '("GCF_014107475.1_ASM1410747v1":0.0655035,("GCF_014771645.1_ASM1477164v1":0.010225,("GCF_016584425.1_ASM1658442v1":0.00859808,"GCF_017896245.1_ASM1789624v1":0.010)0.0519334)0.0191353);'
        for this_genome_id in [
                # wolbachia
                'GCF_014107475.1_ASM1410747v1',
                'GCF_014771645.1_ASM1477164v1',
                'GCF_016584425.1_ASM1658442v1',
                'GCF_017896245.1_ASM1789624v1']:
            this_gbff_filename = this_genome_id + '_genomic.gbff.gz'

            # genome
            gbfffile = tempdir / this_gbff_filename
            copyfile(Path(__file__).parent / 'data' / this_gbff_filename, gbfffile)
            genome_ref = cls.gfu.genbank_to_genome (
                { "workspace_name": cls.wsName,
                  "genome_name": this_genome_id,
                  "file": {"path": str(gbfffile)},
                  "source": "Genbank",
                  "scientific_name": "Genus_foo species_bar",
                  "generate_missing_genes": "True"                
                })['genome_ref']
            cls.genomes.append(genome_ref)

            # prep container obj data
            gs_genome_elements[this_genome_id] = { 'ref': genome_ref }
            tree_ws_refs[this_genome_id] = {'g': [genome_ref]}
            tree_default_node_labels[this_genome_id] = this_genome_id
            tree_leaf_list.append(this_genome_id)


        # save GenomeSet
        genomeSet_name = 'Wolbachia.GenomeSet'
        genomeSet_obj_data = {'description': 'Test GS', 'elements': gs_genome_elements }
        try:
            genomeSet_info = cls.wsClient.save_objects(
                {'workspace': cls.wsName,
                 'objects': [{
                     'type': 'KBaseSearch.GenomeSet',
                     'data': genomeSet_obj_data,
                     'name': genomeSet_name
                     }]})[0]
        except Exception as e:
            raise ValueError ("ABORT: unable to save GenomeSet object.\n"+str(e))
        cls.genomeSet = cls.ref_from_info(genomeSet_info)


        # save Tree
        tree_name = 'Wolbachia.Tree'
        tree_obj_data = { 'name': tree_name,
                          'description': 'test tree',
                          'type': 'SpeciesTree',
                          'tree': tree_newick,
                          'default_node_labels': tree_default_node_labels,
                          'ws_refs': tree_ws_refs,
                          'leaf_list': tree_leaf_list
                        }

        try:
            tree_info = cls.wsClient.save_objects(
                {'workspace': cls.wsName,
                 'objects': [{
                     'type': 'KBaseTrees.Tree',
                     'data': tree_obj_data,
                     'name': tree_name
                     }]})[0]
        except Exception as e:
            raise ValueError ("ABORT: unable to save Tree object.\n"+str(e))
        cls.tree = cls.ref_from_info(tree_info)


    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    # Prepare test objects in workspace if needed using
    # self.getWsClient().save_objects({'workspace': self.getWsName(),
    #                                  'objects': []})
    #
    # Run your method by
    # ret = self.getImpl().your_method(self.getContext(), parameters...)
    #
    # Check returned data with
    # self.assertEqual(ret[...], ...) or other unittest methods


    #### test_run_mmseqs2_and_mptupan_files_01 ():
    #
    # HIDE @unittest.skip("skipped test_run_mmseqs2_and_mptupan_files_01()")  # uncomment to skip
    def test_run_mmseqs2_and_mptupan_files_01 (self):
        method = 'test_run_mmseqs2_and_mptupan_files_01'
        msg = "RUNNING: " + method + "()"
        print("\n\n" + msg)
        print("=" * len(msg) + "\n\n")

        # put test data somewhere to run
        run_dir = os.path.join(self.scratch, 'motupan_test')
        if not os.path.exists(run_dir):
            os.makedirs(run_dir)

        target = 'g__Archaeoglobus'
        output_pg_path = os.path.join (run_dir, target+'-mOTUpan-pangenome.json')
        faa_file = target+'.faa'
        qual_file = target+'.checkm'
        id_map_file = target+'.gene_id_map'
        genome_id_map_file = target+'.genome_name_to_ref.map'
        faa_path = os.path.join (run_dir, faa_file)
        qual_path = os.path.join (run_dir, qual_file)
        id_map_path = os.path.join (run_dir, id_map_file)
        genome_id_map_path = os.path.join (run_dir, genome_id_map_file)
        shutil.copy (os.path.join('data',faa_file), faa_path)
        shutil.copy (os.path.join('data',qual_file), qual_path)
        shutil.copy (os.path.join('data',id_map_file), id_map_path)
        shutil.copy (os.path.join('data',genome_id_map_file), genome_id_map_path)


        params = { 'input_faa_path': faa_path,
                   'input_qual_path': qual_path,
                   'input_gene_id_map_path': id_map_path,
                   'genome_name2ref_path': genome_id_map_path,
                   'run_dir': run_dir,
                   'output_pangenome_json_path': output_pg_path,
                   'mmseqs_cluster_mode': 'easy-cluster',
                   'mmseqs_min_seq_id': 0.0,
                   'mmseqs_min_coverage': 0.8,
                   'motupan_max_iter': 1,
                   'force_redo': 1
        }
        
        ret = self.serviceImpl.run_mmseqs2_and_mOTUpan_files (self.ctx, params)

        print('RESULT:')
        pprint(ret)

        pass


    #### test_run_kb_motupan_genomeset_02 ():
    #
    # HIDE @unittest.skip("skipped test_run_kb_motupan_genomeset_02()")  # uncomment to skip
    def test_run_kb_motupan_genomeset_02 (self):
        method = 'test_run_kb_motupan_genomeset_02'
        msg = "RUNNING: " + method + "()"
        print("\n\n" + msg)
        print("=" * len(msg) + "\n\n")


        input_ref = self.genomeSet
        params = { 'workspace_name': self.wsName,
                   'input_ref': input_ref,
                   'output_pangenome_name': 'foo.Pangenome',
                   'checkm _version': 'CheckM-1',
                   'mmseqs_cluster_mode': 'easy-cluster',
                   'mmseqs_min_seq_id': 0.0,
                   'mmseqs_min_coverage': 0.8,
                   'motupan_max_iter': 1,
                   'pcp_input_genome_ref': self.genomes[0],
                   #'pcp_input_compare_genome_refs': [],
                   #'pcp_input_compare_genome_refs': None,
                   'pcp_input_compare_genome_refs': [self.genomes[0], self.genomes[1], self.genomes[2], self.genomes[3]],
                   'pcp_input_outgroup_genome_refs': [self.genomes[3]],
                   'pcp_save_featuresets': 1,
                   'pcp_genome_disp_name_config': 'obj_name_ver_sci_name',
                   'run_as_test_mode': 1
        }
        
        ret = self.serviceImpl.run_kb_motupan (self.ctx, params)

        print('RESULT:')
        pprint(ret)

        pass


    #### test_run_kb_motupan_tree_03 ():
    #
    # HIDE @unittest.skip("skipped test_run_kb_motupan_tree_03()")  # uncomment to skip
    def test_run_kb_motupan_tree_03 (self):
        method = 'test_run_kb_motupan_tree_03'
        msg = "RUNNING: " + method + "()"
        print("\n\n" + msg)
        print("=" * len(msg) + "\n\n")

        pg_ws_name = 'dylan:narrative_1689301418337'

        
        input_ref = self.tree
        params = { 'workspace_name': self.wsName,
                   'input_ref': input_ref,
                   'output_pangenome_name': 'foo.Pangenome',
                   'checkm _version': 'CheckM-1',
                   'mmseqs_cluster_mode': 'easy-cluster',
                   'mmseqs_min_seq_id': 0.0,
                   'mmseqs_min_coverage': 0.8,
                   'motupan_max_iter': 1,
                   'pcp_input_genome_ref': self.genomes[0],
                   #'pcp_input_compare_genome_refs': [],
                   #'pcp_input_compare_genome_refs': None,
                   'pcp_input_compare_genome_refs': [self.genomes[0], self.genomes[1], self.genomes[2], self.genomes[3]],
                   'pcp_input_outgroup_genome_refs': [self.genomes[3]],
                   'pcp_save_featuresets': 1,
                   'pcp_genome_disp_name_config': 'obj_name_ver_sci_name',
                   'run_as_test_mode': 1

        }
        
        ret = self.serviceImpl.run_kb_motupan (self.ctx, params)

        print('RESULT:')
        pprint(ret)

        pass

    
    #### test_upload_pangenomes ():
    #
    # SKIP THIS!!!  JUST USED FOR UPLOADING PG OBJECTS
    @unittest.skip("skipped test_upload_pangenomes()")  # uncomment to skip
    def test_upload_pangenomes (self):
        method = 'test_upload_pangenomes'
        msg = "RUNNING: " + method + "()"
        print("\n\n" + msg)
        print("=" * len(msg) + "\n\n")

        [OBJID_I, NAME_I, TYPE_I, SAVE_DATE_I, VERSION_I, SAVED_BY_I, WSID_I, WORKSPACE_I, CHSUM_I, SIZE_I, META_I] = range(11)  # object_info tuple

        # workspace
        pg_workspace_name = 'dylan:narrative_1689301418337'

        # config
        obj_prefix = 'GTDB_Bac'
        domain = 'bacteria'
        genome_size = '100-948'
        chunk = '0001-0001'

        target_file = domain+'-species_by_genus-'+genome_size+'.cnt.'+chunk
        target_path = os.path.join(os.sep, 'pangenome', 'lists', target_file)
        work_path = os.path.join(os.sep, 'pangenome', 'docker_work-'+domain, 'genus', 'chunk_'+genome_size, 'genus-spreps')

        genus_names = []
        species_counts = []
        with open (target_path, 'r') as targ_h:
            for targ_line in targ_h:
                [genome_cnt, lineage_str] = targ_line.rstrip().split()
                genus = lineage_str.split(';')[-1]
                genus_names.append (genus)
                species_counts.append (genome_cnt)

        # make sure we have all of them first
        for genus_i,genus in enumerate(genus_names):
            pg_json_path = os.path.join (work_path, genus, genus+'-mOTUpan-pangenome-fxn-prot.json')
            if not os.path.exists(pg_json_path):
                raise ValueError ("no such file {}".format(pg_json_path))

        # now upload them
        print ("UPLOADING PANGENOMES")
        for genus_i,genus in enumerate(genus_names):
            print ("================================\n{} UPLOADING {} {}\n================================".format ((genus_i+1), species_counts[genus_i], genus), flush=True)

            #chunk = 'core'
            #chunk = 'acc-part1'
            #chunk = 'acc-part2'
            chunk = 'acc-part3'
            pg_json_path = os.path.join (work_path, genus, genus+'-mOTUpan-pangenome-fxn-prot-'+chunk+'.json')
            #pg_json_path = os.path.join (work_path, genus, genus+'-mOTUpan-pangenome-fxn-prot.json')
            with open (pg_json_path, 'r') as pg_json_h:
                pg_data = json.load(pg_json_h)

            obj_name = obj_prefix+'-'+genus+'-'+'r214'+'-'+str(species_counts[genus_i])+'genomes'+'-subset_'+chunk+'.Pangenome'
            print ("UPLOADING PG TO OBJ {}".format(obj_name), flush=True)
                
            try:
                pg_info = self.wsClient.save_objects(
                    {'workspace': pg_workspace_name,
                     'objects': [{
                         'type': 'KBaseGenomes.Pangenome',
                         'data': pg_data,
                         'name': obj_name
                     }]})[0]
            except Exception as e:
                raise ValueError ("ABORT: unable to save Pangenome object.\n"+str(e))
            
        pass
    
