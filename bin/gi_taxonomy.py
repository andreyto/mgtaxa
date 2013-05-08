"""Convert an input file with NCBI GI sequence IDs into corresponding taxonomy and aggregate reports"""
from MGT.Common import *
from MGT.Taxa import *
from MGT.TaxaPred import TaxaPred
from MGT.ImmClassifierApp import ImmClassifierApp

def gi_list_to_taxa_pred(gi_list,gi_format,ncbi_taxonomy_dir,pred_out_taxa):
    import re
    re_gi = re.compile(r"[^\t]+\t[^\t]*\bgi\|([0-9]+\b)")
    giToTaxa = loadGiTaxBin(ncbiTaxaDumpDir=ncbi_taxonomy_dir)
    inp = openCompressed(gi_list,"r")
    if gi_format == "ncbi_pipes":
        gi_extractor = lambda line: line.split("gi|")[1].split("|")[0].strip()
    elif gi_format == "ncbi_pipes_blast_m8":
        gi_extractor = lambda line: line.split("\t",2)[1].split("gi|")[1].split("|")[0].strip()
        #gi_extractor = lambda line: re_gi.match(line).group(1)
    else:
        gi_extractor = lambda line: line.strip()
    makeFilePath(pred_out_taxa)
    pred = TaxaPred(pred_out_taxa,mode="w")
    pred.initForAppend(idSampDtype=int,expectedrows=10**6)
    data = pred.getData()
    gi_rej = []
    n_unresolved = 0
    gi_rej_buf_size = 10**0
    def _resolve_gi_to_taxid_entrez(gi_rej,predTaxid):
        print ("Attempting to resolve taxids for {0} missing records by calling NCBI EUtils Web server. "+\
                "This may fail or block with a potentially long time out depending on a quality of your Internet "+\
                "connection or NCBI server load.").format(len(gi_rej))
        from MGT.Entrez import EzRequest 
        id_rej = set( (g_r[1] for g_r in gi_rej) )
        req = EzRequest(rettype="xml")
        gi_rej_to_taxid = dict( (rec["Gi"],rec["TaxId"]) for rec in req.summary(ids=id_rej) )
        print "DEBUG: gi_rej_to_taxid = ", gi_rej_to_taxid
        n_unresolved = 0
        for (i_samp,gi_samp) in gi_rej:
            if gi_samp in gi_rej_to_taxid:
                predTaxid[i_samp] = gi_rej_to_taxid[gi_samp]
                print "DEBUG: predTaxid[{0}] = gi_rej_to_taxid[{1}]".format(i_samp,gi_samp)
            else:
                n_unresolved += 1
        return n_unresolved

    idSamp = []
    lenSamp = []
    predTaxid = []
    predScore = []
    i_samp_batch = 0
    flush_buff_size = 10**5
    for i_samp,line in enumerate(inp):
        gi = int(gi_extractor(line))
        
        taxid = rejTaxid
        
        try:
            taxid = giToTaxa[gi]
            if not taxid:
                taxid = rejTaxid
                print "WARNING: GI corresponds to undefined taxid: {0}".format(gi)
        except:
            print "WARNING: GI is larger than the max GI in our taxonomy DB:{0}".format(gi)
        
        if taxid == rejTaxid:
            gi_rej.append((i_samp_batch,gi))

        idSamp.append(i_samp)
        lenSamp.append(1)
        predTaxid.append(taxid)
        predScore.append(1)
        
        i_samp_batch += 1
        
            
        if i_samp % flush_buff_size == 0:
            if len(gi_rej):
                n_unresolved += _resolve_gi_to_taxid_entrez(gi_rej=gi_rej,predTaxid = predTaxid)
                gi_rej = []
            print "Read {0} input records. Taxid not found for {1} records".format(i_samp,n_unresolved)
            data.idSamp.append(idSamp)
            data.lenSamp.append(lenSamp)
            data.predTaxid.append(predTaxid)
            data.predScore.append(predScore)
            idSamp = []
            lenSamp = []
            predTaxid = []
            predScore = []
            i_samp_batch = 0
            pred.flush()
    
    if len(gi_rej):
        n_unresolved += _resolve_gi_to_taxid_entrez(gi_rej=gi_rej,predTaxid = predTaxid)
        gi_rej = []
    
    data.idSamp.append(idSamp)
    data.lenSamp.append(lenSamp)
    data.predTaxid.append(predTaxid)
    data.predScore.append(predScore)
    idSamp = []
    lenSamp = []
    predTaxid = []
    predScore = []
    i_samp_batch = 0
    pred.flush()

    pred.close()
    inp.close()

    print "Read {0} input records. Taxid not found for {1} records".format(i_samp,n_unresolved)


def gi_taxonomy(gi_list,out_dir="results",gi_format="ncbi_pipes_blast_m8",ncbi_taxonomy_dir=None):
    import time
    out_dir = os.path.abspath(out_dir)
    pred_out_taxa=pjoin(out_dir,"taxa_samp")
    start_time = time.time()  
    gi_list_to_taxa_pred(gi_list,gi_format,ncbi_taxonomy_dir,pred_out_taxa)
    end_time = time.time()
    print "DEBUG: loaded hit list in {0} sec.".format(end_time-start_time)
    opt = Struct()
    opt.runMode = "inproc"
    opt.web = False
    opt.needTerminator = False
    opt.mode = "export-predictions"
    opt.predOutDir = out_dir
    opt.predOutTaxa = pred_out_taxa
    opt.predMinLenSamp = 1
    opt.taxaTreeNcbiDir = os.path.abspath(ncbi_taxonomy_dir)
    opt.predOutStatsKronaEmbed = "krona"

    ImmClassifierApp.fillWithDefaultOptions(opt)
    imm = ImmClassifierApp(opt=opt)
    imm.run()

def main():
    import argh
    argh.dispatching.dispatch_command(gi_taxonomy)

if __name__ == '__main__':
    main()
