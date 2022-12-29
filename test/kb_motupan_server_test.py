# -*- coding: utf-8 -*-
import os
import time
import unittest
import shutil
from pprint import pprint
from configparser import ConfigParser

from kb_motupan.kb_motupanImpl import kb_motupan
from kb_motupan.kb_motupanServer import MethodContext
from installed_clients.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace


class kb_motupanTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_motupan'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_motupan',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = kb_motupan(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_ContigFilter_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa


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
        faa_path = os.path.join (run_dir, faa_file)
        qual_path = os.path.join (run_dir, qual_file)
        id_map_path = os.path.join (run_dir, id_map_file)
        shutil.copy (os.path.join('data',faa_file), faa_path)
        shutil.copy (os.path.join('data',qual_file), qual_path)
        shutil.copy (os.path.join('data',id_map_file), id_map_path)


        params = { 'input_faa_path': faa_path,
                   'input_qual_path': qual_path,
                   'input_gene_id_map_path': id_map_path,
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
