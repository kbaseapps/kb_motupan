/*
A KBase module: kb_motupan
*/

module kb_motupan {

    /* 
    ** The workspace object refs are of form:
    **
    **    objects = ws.get_objects([{'ref': params['workspace_id']+'/'+params['obj_name']}])
    **
    ** "ref" means the entire name combining the workspace id and the object name
    ** "id" is a numerical identifier of the workspace or object, and should just be used for workspace
    ** "name" is a string identifier of a workspace or object.  This is received from Narrative.
    */
    typedef string workspace_name;
    typedef string file_path;
    typedef string data_obj_name;
    typedef string data_obj_ref;
    typedef int    bool;

    typedef structure {
        data_obj_name report_name;
        data_obj_ref  report_ref;
    } ReportResults;


    /* run_mmseqs2_and_mOTUpan_files()
    **
    **  Method for running mmseqs2 and mOTUpan from files
    */
    typedef structure {
	file_path input_faa_path;
	file_path input_qual_path;
	file_path input_gene_id_map_path;
	file_path run_dir;
	file_path output_pangenome_json_path;

	string mmseqs_cluster_mode;
	float  mmseqs_min_seq_id;
	/*int    mmseqs_cov_mode;*/
	float  mmseqs_min_coverage;
	int    motupan_max_iter;
    } run_mmseqs2_and_mOTUpan_files_Params;

    typedef structure {
	file_path pangenome_json;
    } run_mmseqs2_and_mOTUpan_files_Output;

    funcdef run_mmseqs2_and_mOTUpan_files (run_mmseqs2_and_mOTUpan_files_Params params)  returns (run_mmseqs2_and_mOTUpan_files_Output) authentication required;

};
