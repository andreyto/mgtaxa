from MGT.Common import *

def test_gi_taxonomy():
    inp_gi_list = pjoin(options.testDataDir,"gbank","ncbi_ids.tab")
    #taxonomy_dir = options.taxaDumpDir
    taxonomy_dir = '/usr/local/projects/HMP/atovtchi/taxonomy.2013-04-29'
    cmd = ("$MGT_HOME/bin/mgtaxa bin/gi_taxonomy.py "+\
        "--out-dir gi_taxonomy.results "+\
        "--ncbi-taxonomy-dir {0} "+\
        "{1}").\
        format(taxonomy_dir,inp_gi_list)
    run(cmd,debug=True)

if __name__ == "__main__":
    test_gi_taxonomy()

