/*
Reference to a Genome object in the workspace
@id ws KBaseGenomes.Genome KBaseGenomeAnnotations.GenomeAnnotation
*/
typedef string Genome_ref;
typedef string Genome_name;

/*
OrthologFamily object: this object holds all data for a single ortholog family in a metagenome.

Fields:
    id - string - group identifier
    type - string - ... (may duplicate "cat"? leave alone for backwards)
    gene_name - list<string> - name as described in KBaseGenomes.Genome
    function - string - functions as described in KBaseGenomes.Genome
    function_sources - list<tuple<string,string>> - list of tuples of:
        (0) string - gene identifier (ID in gff file)
        (1) string - genome workspace reference
    function_logic - string - controlled vocab of:
        (1) 'single_source'
        (2) 'union'
        (3) 'intersection'
    md5 - string - md5 encoded string of protein_translation
    protein_translation - string - protein translation string
    protein_translation_source - tuple<string,string> - tuple of:
        (0) string - gene identifier (ID in gff file)
        (1) string - genome workspace reference
    orthologs - list<tuple<string,float,string>> - list of tuples of:
        (0) string - gene identifier (ID in gff file)
        (1) float - numerical order in gff file OR gene order in BLAST
        (2) string - genome workspace reference
    genome_occ - int - number of input genomes in cluster
    cat - string - controlled vocabulary of cluster classification -
        depends on method (see Pangenome struct)
    core_log_likelihood - float - mOTUpan-specific score
    mean_copies - float - avg number of copies per genome with cluster

@optional type gene_name function function_sources function_logic md5 protein_translation protein_translation_source genome_occ cat core_log_likelihood mean_copies
*/
typedef structure {
  string id;
  string type;
  list<string> gene_name;
  string function;
  list<tuple<string,string>> function_sources;
  string function_logic;
  string md5;
  string protein_translation;
  tuple<string,string> protein_translation_source;
  list<tuple<string, float, string>> orthologs;
  int genome_occ;
  string cat;
  float core_log_likelihood;           
  float mean_copies;

} OrthologFamily;


/*
Pangenome object: this object holds all data regarding a pangenome

@searchable ws_subset id name
@metadata ws type as Type
@metadata ws name as Name
@metadata ws length(orthologs) as Number orthologs
@metadata ws length(genome_refs) as Number genomes

Fields:
    id - string - pangenome identifier
    name - string - pangenome name
    type - string - method used to calculate pangenome -
        controlled vocab:
        (1) 'OrthoMCL'
        (2) 'mOTUpan'
    type_ver - string - version of pangenome calc method
    cluster_cats - mapping<string,mapping<string,string>> -
        controlled vocab keyed by pg calc type to standard terms
        current standard terms: 'core','flexible','singleton'
        e.g.:
        cluster_cats = { 'mOTUpan': { 'core': 'core',
                                      'accessory': 'flexible'
                                    },
                       { 'OrthoMCL': { 'core': 'core',
                                       'peripheral': 'flexible',
                                       'singleton': 'singleton'
                                     }
                       }
    clustering_method - method for clustering genes -
        controlled vocab:
        (1) OrthoMCL
        (2) MMseqs2
    clustering_method_ver - string - version of clustering method
    clustering_method_params - mapping<string,string> - clust params
    pangenome_method_params - mapping<string,string> - pg calc params
    genome_count - int - number of input genomes
    core_length - int - number of clusters considered 'core'
    mean_est_genome_size - float - estimate of avg genes per genome
    prior_genome_completeness - mapping<Genome_name,float> -
        percentage value dict with genome_ref as key
    posterior_genome_completeness - mapping<Genome_name,float> -
        percentage value dict with genome_ref as key
    genome_refs - list<string> - input genome refs
    genome_names - list<string> - input genome names
    genome_ref_to_name - mapping<Genome_ref,Genome_name> - lookup dict
    genome_name_to_ref - mapping<Genome_name,Genome_ref> - lookup dict
    orthologs - list<OrthologFamily> - pangenome clusters

@optional type_ver cluster_cats clustering_method clustering_method_ver clustering_method_params pangenome_method_params genome_count core_length mean_est_genome_size prior_genome_completeness posterior_genome_completeness genome_names genome_ref_to_name genome_name_to_ref 
*/
typedef structure {
  string id;
  string name;
  string type;
  string type_ver;
  mapping<string,mapping<string,string>> cluster_cats;
  string clustering_method;
  string clustering_method_ver;
  mapping <string,string> clustering_method_params;
  mapping <string,string> pangenome_method_params;
  int genome_count;
  int core_length;
  float mean_est_genome_size;
  mapping<Genome_name,float> prior_genome_completeness;
  mapping<Genome_name,float> posterior_genome_completeness;
  list<Genome_ref> genome_refs;
  list<Genome_name> genome_names;
  mapping<Genome_ref,Genome_name> genome_ref_to_name;
  mapping<Genome_name,Genome_ref> genome_name_to_ref;
  list<OrthologFamily> orthologs;
} Pangenome;
