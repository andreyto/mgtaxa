from MGT.Taxa import *
from MGT.Kmers import *


class CountAggregateVisitor:
    
    def __call__(self,node):
        for child in node.children:
            node.data.kmerCnt += child.data.kmerCnt



class TaxaSampler(Options):
    def __init__(self):
        PhyOptions.__init__(self)
        #options
        self.kmerSize = 6 #7
        self.kmerWindow = 1000
        self.kmerPath = '6mers_1K' #'7mers_1K'
        self.maxKmerRankTrain = 500
        self.minTestSimSeqLen = 750
        self.maxTestSimSeqLen = 3000
        self.kmerPathVir = self.kmerPath+'.vir'
        self.kmerTxtPathVir = self.kmerPathVir + '.gz'        
        self.kmerPathNonVir = self.kmerPath+'.nvir'
        self.kmerPathVirNonVir = self.kmerPath+'.vnv'
        self.trainSetPathVirNonVir = self.kmerPath+'.vnv-tr'
        self.trainSvmPathVirNonVir = self.trainSetPathVirNonVir + '.svm.inp'
        self.testSvmPathVirNonVir = self.kmerPath+'.vnv-ts.svm.inp'
        self.labeledSetPathVirFamily = self.kmerPath+'.v_f-lb'
        self.trainSetPathVirFamily = self.kmerPath+'.v_f-tr'
        self.testSetPathVirFamily = self.kmerPath+'.v_f-ts'
        self.testIndPathVirFamily = self.kmerPath+'.v_f-ts.ind'
        self.testSeqSimPathVirFamily = self.kmerPath+'.v_f-ts.sim.fasta.gz'
        self.testKmerTxtSimPathVirFamily = self.kmerPath+'.v_f-ts.sim.gz'
        self.testSimPathVirFamily = self.kmerPath+'.v_f-ts.sim'
        self.testSvmSimPathVirFamily = self.testSimPathVirFamily + '.svm.inp'
        self.knownSeqFasta = "phyla_sel.fasta.gz"
        self.knownSeqVirFasta = "phyla_sel.vir.fasta.gz"
        self.knownSeqVirKmerFasta = "phyla_sel.vir.kmer.fasta.gz"
        self.knownSeqSampFasta = "phyla_sel.samp.fasta.gz"
        self.knownSeqSampKmerFasta = "phyla_sel.samp.kmer.fasta.gz"
        self.kmerPathSamp = self.kmerPath+'.samp'
        self.kmerTxtPathSamp = self.kmerPathSamp + '.gz'
        self.viralRootTaxid = viralRootTaxid
        #data
        self.taxaTree = None


    def loadTaxaTree(self):
        if self.taxaTree is None:
            print "Loading taxonomy tree"
            self.taxaTree = TaxaTree(ncbiDumpFile=self.taxaNodesFile,save=False,load=False)
            print "Finished loading taxonomy tree"
    
    def extractViralKmers(self):
        self.loadTaxaTree()
        mask = self.taxaTree.buildSubtreeMask(ids=(self.viralRootTaxid,), values=(1,))
        kmers = KmerBinReader(self.kmerPath)
        kmersSel = KmerBinWriter(self.kmerPathVir)
        for batch in kmers.readBatches(batchSize=10000):
            sel = batch[numpy.where(mask.take(batch['taxid']) > 0)]
            kmersSel.writeBatch(sel)
            print "Read %s k-mer records out of %s, selected %s" % (kmers.posRec(),kmers.numRec(),kmersSel.numRec())
        kmersSel.close()
        kmers.close()
        kmersSel.saveCounts()
        
    def extractViralSequence(self):
        #Load kmers taxids for sanity check only, see assertion below
        kmers = KmerBinReader(self.kmerPathVir)
        counts = kmers.readCounts()
        taxidKmers = set(numpy.where(counts>0)[0])
        self.loadTaxaTree()
        maskVir = self.taxaTree.buildSubtreeMask(ids=(self.viralRootTaxid,), values=(1,))
        taxidTree = set(numpy.where(maskVir==1)[0])
        print "#taxidKmers = %s, #taxidTree = %s" % (len(taxidKmers),len(taxidTree))
        assert taxidKmers <= taxidTree, "kmer viral taxids must be a subset of all viral taxids"
        print "Extracting all sequence records for %s taxonomy ids." % (len(maskVir[maskVir==1]),)
        reader = fastaReaderGzip(self.knownSeqFasta)
        writer = openGzip(self.knownSeqVirFasta,'w')
        iRecIn = 0
        iRecOut = 0
        for rec in reader.records():
            hdr = rec.header()
            taxid = taxidFromPhyFastaHeader(hdr)
            if maskVir[taxid] == 1:
                writer.write(hdr)
                for line in rec.seqLines():
                    writer.write(line)
                iRecOut += 1
            iRecIn += 1
            if iRecIn % 10000 == 0:
                print "Read %s fasta records, wrote %s." % (iRecIn,iRecOut)
        reader.close()
        writer.close()
        
    def extractSampleSequence(self):
        reader = fastaReaderGzip(self.knownSeqFasta)
        writer = openGzip(self.knownSeqSampFasta,'w')
        iRecIn = 0
        iRecOut = 0
        for rec in reader.records():
            hdr = rec.header()
            taxid = taxidFromPhyFastaHeader(hdr)
            if nrnd.ranf() > 0.9:
                writer.write(hdr)
                seqLen = 0
                maxSeqLen = 50000
                for line in rec.seqLines():
                    if seqLen + len(line) - 1 > maxSeqLen:
                        line = line[:maxSeqLen-seqLen] + '\n'
                    writer.write(line)
                    seqLen += len(line) - 1
                    if seqLen >= maxSeqLen:
                        break
                iRecOut += 1
            iRecIn += 1
            if iRecIn % 10000 == 0:
                print "Read %s fasta records, wrote %s." % (iRecIn,iRecOut)
        reader.close()
        writer.close()

    def makeSampleKmers(self):
        joinFastaByTaxid(self.knownSeqSampFasta,self.knownSeqSampKmerFasta)
        makeKmersTxt(self.knownSeqSampKmerFasta,self.kmerTxtPathSamp,self.kmerSize,self.kmerWindow)
        kmersTxtToBin(self.kmerTxtPathSamp)


        
    def extractNonViralSequence_bak(self):
        kmers = KmerBinReader(self.kmerPathNonVir)
        counts = kmers.readCounts()
        taxids = set(numpy.where(counts>0)[0])
        print "Extracting all sequence records for %s taxonomy ids." % (len(taxids),)
        reader = fastaReaderGzip(self.knownSeqFasta)
        writer = openGzip(self.knownSeqNonVirFasta,'w')
        iRecIn = 0
        iRecOut = 0
        for rec in reader.records():
            hdr = rec.header()
            taxid = taxidFromPhyFastaHeader(hdr)
            if taxid < len(counts) and counts[taxid] > 0:
                seqLen = 0
                maxSeqLen = self.kmerWindow * counts[taxid]
                lines = []
                for line in rec.seqLines():
                    if seqLen + len(line) - 1 > maxSeqLen:
                        line = line[:maxSeqLen-seqLen] + '\n'
                    lines.append(line)
                    seqLen += len(line) - 1
                    if seqLen >= maxSeqLen:
                        break
                kmersDone = int(seqLen / self.kmerWindow)
                if kmersDone > 0:
                    writer.write(hdr)
                    for line in lines:
                        writer.write(line)
                    counts[taxid] -= kmersDone
                    #print taxid, kmersDone, counts[taxid]
                    iRecOut += 1
            iRecIn += 1
            if iRecIn % 10000 == 0:
                print "Read %s fasta records, wrote %s." % (iRecIn,iRecOut)
        reader.close()
        writer.close()
        
    def assignNodeCounts(self):
        self.loadTaxaTree()
        counts = self.assignLeafCounts()
        count_aggregator = CountAggregateVisitor()
        self.taxaTree.visitDepthBottom(count_aggregator)
        counts[:] = 0
        def count_setter(node):
            counts[node.id] = node.data.kmerCnt
        self.taxaTree.visitDepthTop(count_setter)
        debug = False
        if debug:
            def node_printer(node):
                if node.data.rank == 'family':
                    print node.data.kmerCnt, node.lineageRanksStr()
            self.taxaTree.visitDepthTop(node_printer,self.viralRootTaxid)
        return counts
    
    def assignLeafCounts(self):
        kmerCounts = KmerCounts(self.kmerPath).load()
        nodes = self.taxaTree.getNodesDict()
        for node in nodes.itervalues():
            try:
                x = kmerCounts[node.data.taxid]
            except IndexError:
                x = 0
            node.data.kmerCnt = x
        return kmerCounts

    def subSampleKmersNonViral(self):
        self.loadTaxaTree()
        subCounts = self.subSampleCountsNonViral()
        subCountsSel = numpy.zeros(len(subCounts),dtype=int)
        kmers = KmerBinReader(self.kmerPathSamp)
        kmersSel = KmerBinWriter(self.kmerPathNonVir)
        for batch in kmers.readBatches(batchSize=10000):
            mask = numpy.zeros(len(batch['taxid']),dtype=bool)
            for (taxid,i_batch) in itertools.izip(batch['taxid'],itertools.count()):
                #pdb.set_trace()
                if subCountsSel[taxid] < subCounts[taxid]:
                    mask[i_batch] = 1
                    subCountsSel[taxid] += 1
            sel = batch[numpy.where(mask > 0)]
            kmersSel.writeBatch(sel)
            print "Read %s k-mer records out of %s, selected %s" % (kmers.posRec(),kmers.numRec(),kmersSel.numRec())
        kmersSel.close()
        kmers.close()
        kmersSel.saveCounts(subCountsSel)
        assert numpy.all( ( subCountsSel - subCounts ) == 0 )
        assert numpy.sum(subCountsSel) == kmersSel.numRec()
        
            
        
    def subSampleCountsNonViral(self):
        maskVir = self.taxaTree.buildSubtreeMask(ids=(self.viralRootTaxid,), values=(1,))
        kmerCounts = numpy.wrongResize(KmerCounts(self.kmerPathSamp).load(),len(maskVir))
        kmerCountsVir = numpy.wrongResize(KmerCounts(self.kmerPathVir).load(),len(maskVir))
        #assert len(maskVir) == len(kmerCounts)
        kmerSumVir = numpy.sum(kmerCountsVir.compress(maskVir > 0))
        kmerSumNonVir = kmerSumVir
        ind = numpy.arange(len(kmerCounts),dtype=int)
        indPerm = nrnd.permutation(ind)
        kmerCountsSel = numpy.zeros(len(kmerCounts),dtype=int)
        nSel = 0
        nSelPrev = 0
        done = False
        while not done:
            nSelPrev = nSel
            #TODO: should we randomize permutation at every pass?
            for taxid in indPerm: 
                if maskVir[taxid] != 1 and kmerCountsSel[taxid] < kmerCounts[taxid]:
                    kmerCountsSel[taxid] += 1
                    nSel += 1
                    if nSel >= kmerSumNonVir:
                        done = True 
                        break
            if nSel == nSelPrev:
                done = True
        #pdb.set_trace()
        assert numpy.sum(maskVir[kmerCountsSel.nonzero()] > 0) == 0, "POSTCONDITION: no viruses"
        return kmerCountsSel
        
    def writeViralVsNonViralTrainingSet(self):
        inpPaths = (self.kmerPathVir, self.kmerPathNonVir)
        labels = (1,-1)
        indTestVir = loadObj(self.testIndPathVirFamily)
        cntTaxidTestVir = indTestVir.cntTaxidTest
        kmersOut = KmerBinWriter(self.trainSetPathVirNonVir)
        for (inpPath,label) in zip(inpPaths,labels):
            kmersInp = KmerBinReader(inpPath)
            counts = kmersInp.readCounts()
            cntTaxidVir = numpy.wrongResize(max(len(cntTaxidVir),len(counts)))
            for batch in kmersInp.readBatches(batchSize=10000):
                batchSel = batch[cntTaxidTestVir[batch['taxid']] <= 0]
                batchSel['taxid'][:] = label
                kmersOut.writeBatch(batchSel)
                print "For label %s, read %s k-mer records out of %s, selected %s" % (label,kmersInp.posRec(),kmersInp.numRec(),kmersOut.numRec())
            kmersInp.close()
        kmersOut.close()

    def writeViralFamilyLabels(self):
        otherLabel = 0
        testRatio = 0.2
        minKmerRankTrain = 15
        maxKmerRankTrain = self.maxKmerRankTrain
        self.loadTaxaTree()
        kmersInp = KmerBinReader(self.kmerPathVir)        
        labelData = self.buildViralRankLabels(rank='family',
                                              indTaxid=kmersInp.readIndTaxid(),
                                              testRatio=testRatio,
                                              minKmerRankTrain=minKmerRankTrain,
                                              maxKmerRankTrain=maxKmerRankTrain)
        #dumpObj(labelData.indices,self.testIndPathVirFamily)
        #return
        kmersOutTrain = KmerToSvmWriter(self.trainSetPathVirFamily)
        kmersOutTest = KmerToSvmWriter(self.testSetPathVirFamily)
        iRec = 0
        for batch in kmersInp.readBatches(batchSize=10000):
            taxid = batch['taxid']
            iRecEnd = iRec + len(batch)
            indRank = labelData.indRank[iRec:iRecEnd]
            indTrain = labelData.indTrain[iRec:iRecEnd]
            taxid[:] = indRank
            batchSel = batch[indTrain>0]
            kmersOutTrain.writeBatch(batchSel)
            indTest = labelData.indTest[iRec:iRecEnd]
            batchSel = batch[indTest>0]
            kmersOutTest.writeBatch(batchSel)
            print "Read %s k-mer records out of %s, wrote %s training and %s testing samples." % \
            (kmersInp.posRec(),kmersInp.numRec(),kmersOutTrain.numRec(),kmersOutTest.numRec())
            iRec += len(batch)
        kmersInp.close()
        kmersOutTrain.close()
        kmersOutTest.close()
        dumpObj(labelData.indices,self.testIndPathVirFamily)

    def buildViralRankLabels(self,rank,indTaxid,testRatio,minKmerRankTrain,maxKmerRankTrain,permutMethod='sort',
                             trainWithNonSpecies=True,fixedTestCnt=False):
        """'permutMethod'='sort' will first select for testing species with smaller kmer counts, which will allow more
        fine grain control over the 'testRatio' actual selected value, because many ranks have very uneven distribution
        of kmer counts over species. The downside is a possible testing bias because smaller species counts are typical
        for partial sequence (e.g. protein coding sequences). 'permutMethod'='random' will shuffle the list of species
        and then select from left to right."""
        exclusionRank = "genus"
        #maskVir = self.taxaTree.buildSubtreeMask(ids=(self.viralRootTaxid,), values=(1,))
        #maskRank[0:maxTaxid+1] ( any taxid -> taxid of 'rank' (e.g. family taxid) )
        maskRank = self.taxaTree.buildRankSubtreeMask(rank=rank,other=0)
        #maskSpe[0:maxTaxid+1]  (any taxid -> taxid of species)
        maskSpe = self.taxaTree.buildRankSubtreeMask(rank=exclusionRank,other=0)
        #kmerCnt[0:maxTaxid+1] (any taxid -> number of kmer records with this taxid (itself, not its subtree))
        ##kmerCnt = numpy.wrongResize(self.assignNodeCounts(),len(maskRank))
        #kmerCnt calculated on actual kmer file, so it might not have high taxids if all of them have zero counts,
        #so we extend it with zeros to match size of mask arrays
        kmerCnt = numpy.bincount(indTaxid)
        kmerCnt = numpy.wrongResize(kmerCnt,len(maskRank))
        #[0:maxRankTaxid+1] (rank taxid -> kmer count in subtree of this taxid)
        rankKmerCnt = numpy.bincount(maskRank,kmerCnt).astype(int)
        #[0:maxSpeTaxid+1] (species taxid -> kmer count in subtree of this taxid)
        speKmerCnt = numpy.bincount(maskSpe,kmerCnt).astype(int)
        #Build index as a dict of sets (rank Taxid -> (rank's species Taxid's))
        indRankToSpe = {}
        #we ignore all 'rank' that are not represented by sufficient amount of current kmers
        for (taxidRank,taxidSpe) in ( x for x in numpy.column_stack((maskRank,maskSpe)) 
                                      if x[0] != 0 and x[1] != 0 and rankKmerCnt[x[0]] > minKmerRankTrain):
            try:
                indRankToSpe[taxidRank].add(taxidSpe)
            except KeyError:
                indRankToSpe[taxidRank] = set((taxidSpe,))
        rankSpeCnt = numpy.asarray([len(speSet) for speSet in indRankToSpe.itervalues()],numpy.int32)
        targetKmerRankTest = numpy.ceil(testRatio * numpy.median(rankKmerCnt[rankKmerCnt>0]))
        #from each family, randomly select 20% of its species for test set
        #but no less than 15 k-mers for training. If less than 15 k-mers at all, skip this
        #family
        #
        #[0:maxSpeTaxId] (species taxid -> 1 if selected for training set, 0 - otherwise
        speMaskTrain = numpy.zeros_like(speKmerCnt)
        #same as above but for testing set
        speMaskTest = numpy.zeros_like(speMaskTrain)
        ##[0:maxRankTaxid+1] ( rank taxid -> label (equal to rank taxid except than set to 0 for low kmer count )
        #rankToLabel = numpy.arange(len(rankKmerCnt))
        for (iRank,speSet) in indRankToSpe.iteritems():
            species = array(list(speSet),numpy.int32)
            if permutMethod == 'random':
                spePermut = nrnd.permutation(species)
                spePermutKmerCnt = speKmerCnt[spePermut]
            elif permutMethod == 'sort':
                spePermutKmerCnt = speKmerCnt[species]
                indOrder = spePermutKmerCnt.argsort()
                spePermut = species[indOrder]
                spePermutKmerCnt = spePermutKmerCnt[indOrder]
                
            spePermutKmerCntAccum = numpy.add.accumulate(spePermutKmerCnt)
            #we also consider kmers w/o species as valid 'rank' training sets,
            #and decrease the species defined minCount respectively
            minKmerRankTrainSpe = minKmerRankTrain - (rankKmerCnt[iRank] - spePermutKmerCntAccum[-1])
            maxKmerRankTestSpe = rankKmerCnt[iRank] - minKmerRankTrain
            targetKmerRankTestSpe = min(targetKmerRankTest,maxKmerRankTestSpe)
            if not fixedTestCnt:
                targetKmerRankTestSpe = min(targetKmerRankTestSpe,numpy.ceil(rankKmerCnt[iRank]*testRatio))
            speIndTargetCountTest = numpy.searchsorted(spePermutKmerCntAccum,targetKmerRankTestSpe,side='right')
            speTest = spePermut[:speIndTargetCountTest]
            speTrain = spePermut[speIndTargetCountTest:]
            speMaskTrain[speTrain] = 1
            speMaskTest[speTest] = 1
        
        #[0:Total kmer count] (kmer ind -> taxid of species)
        indSpe = maskSpe[indTaxid]
        #same as above, but for 'rank'. This will be classification labels.
        indRank=maskRank[indTaxid]
        #[0:Total kmer count] (kmer ind -> kmer count for 'rank')
        indRankCnt = rankKmerCnt[indRank]
        #[0:Total kmer count] (kmer ind -> 1 if kmer in testing set, 0 - otherwise)
        indTest = speMaskTest[indSpe]
        #same as above but for training set: all that is not testing and whose 'rank' has kmer count not too small
        #This will include into training set 'rank' members which do not have a species defined.
        #The alternative branch will omit 'unclassified' species from training set.
        #They are never included into testing set.
        if trainWithNonSpecies:
            indTrain = numpy.logical_and(numpy.logical_and(indRankCnt >= minKmerRankTrain,numpy.logical_not(indTest)),
                                         indRank>0)
        else:
            indTrain = speMaskTrain[indSpe]
        
        indTrain = self.clipMask(indRank, indTrain, maxKmerRankTrain)
        
        labelData = Struct(indTrain=indTrain,indTest=indTest,indRank=indRank)
        labelCounts = self.reportLabelCounts(labelData,indTaxid)
        labelData.indices = labelCounts
        print "Label counts:\n", labelCounts
        #pdb.set_trace()
        return labelData



    def reportLabelCounts(self,labelData,indTaxid):
        indRank = labelData.indRank
        indTrain = labelData.indTrain
        indTest = labelData.indTest
        #rank[taxid], non-zero only for taxid in test set
        rankTaxidTest = masksToInd(indRank*indTest, indTaxid)
        cntTaxidTest = numpy.bincount(indTaxid,indTest).astype(int)
        assert len(cntTaxidTest) <= len(rankTaxidTest)
        assert numpy.all(rankTaxidTest[len(cntTaxidTest):] == 0)
        cntTaxidTestWhere=whereItems(cntTaxidTest,cntTaxidTest>0)
        rankSelCnt = numpy.rec.fromarrays((numpy.bincount(indRank,indTrain).astype(int),numpy.bincount(indRank,indTest).astype(int)),names='train,test')
        rankSelCnt = numpy.sort(rankSelCnt[rankSelCnt['train']>0])
        assert numpy.sum(cntTaxidTestWhere['val']) == numpy.sum(rankSelCnt['test'])
        return Struct(rankSelCnt=rankSelCnt,
                      medianRankSelCnt=Struct(train=numpy.median(rankSelCnt['train']),test=numpy.median(rankSelCnt['test'])),
                      lenRankSelCnt=Struct(train=len(rankSelCnt['train']),test=len(numpy.where(rankSelCnt['test']>0)[0])),
                      rankTaxidTest=rankTaxidTest,
                      cntTaxidTest=cntTaxidTest,
                      cntTaxidTestWhere=cntTaxidTestWhere)

    def clipMask(self,indLabel,maskSelect,maxLabelCnt):
        #ATTENTION: I used var name 'counts' at first instead of the
        #'labelCounts' and 'counts' always came out all zero no matter what.
        #Must be a bug-like side effect due to some name clash.
        labelCounts = numpy.bincount(indLabel,maskSelect)
        #targetCount = numpy.median(labelCounts[labelCounts>0])
        labelCounts = labelCounts.clip(0,maxLabelCnt)
        #shuffle the labels, but remember the original positions
        labelsPerm = nrnd.permutation(numpy.rec.fromarrays((numpy.arange(len(indLabel)),
                                                            indLabel,
                                                            numpy.zeros_like(maskSelect)),
                                                            names='id,label,select'))
        countsSel = numpy.zeros_like(labelCounts)
        for item in labelsPerm:
            if countsSel[item['label']] < labelCounts[item['label']]:
                item['select'] = 1
                countsSel[item['label']] += 1
        newMaskSelect = numpy.zeros_like(maskSelect)
        newMaskSelect[labelsPerm['id']] = labelsPerm['select']
        return newMaskSelect


    
    def sampleTestByPredictLength(self,fastaFilePredict):
        reader = fastaReaderGzip(fastaFilePredict)
        bins =  numpy.arange(self.minTestSimSeqLen,
                             self.maxTestSimSeqLen,
                             (self.maxTestSimSeqLen-self.minTestSimSeqLen)/20.0,
                             dtype=float)
        genSeqLen = seqLengthDistribFromSample(reader.records(),bins=bins)
        reader.close()
        reader = fastaReaderGzip(self.knownSeqVirFasta)
        indices = loadObj(self.testIndPathVirFamily)
        writerChunks = openGzip(self.testSeqSimPathVirFamily,'w')
        taxidLast = -1
        seqTaxidStream = StringIO()
        for rec in reader.records():
            taxid = taxidFromPhyFastaHeader(rec.header())
            if taxid != taxidLast:
                seqTaxid = seqTaxidStream.getvalue()
                seqTaxidStream.close()
                seqTaxidStream = StringIO()
                if taxidLast > 0:
                    #now process the full taxid sequence
                    nChunks = 0
                    lenChunks = 0
                    lenSeqTaxid = len(seqTaxid)
                    nKmers = indices.cntTaxidTest[taxidLast]
                    while nChunks < nKmers:
                        lenChunk = int(numpy.round(genSeqLen()))
                        writerChunks.write(">%s\n" % (taxidLast,))
                        writeSeqByLines(writerChunks,seqTaxid[lenChunks:lenChunks+lenChunk])
                        nChunks += 1
                        lenChunks += lenChunk
                        if lenChunks > lenSeqTaxid:
                            #print "Warning: requested chunks length exceeded available sequence length. "+\
                            #"taxid:%s nChunks:%s nKmers:%s lenSeqTaxid:%s lenChunks:%s" % \
                            #(taxidLast,nChunks,nKmers,lenSeqTaxid,lenChunks)
                            break
                taxidLast = taxid
            
            for line in rec.seqLines():
                seqTaxidStream.write(line.rstrip("\n"))
        reader.close()
        writerChunks.close()
        lenHistPred = genSeqLen.histogram()
        lenHistPredNorm = lenHistPred[0].astype(float)/numpy.sum(lenHistPred[0])        
        reader = fastaReaderGzip(self.testSeqSimPathVirFamily)
        lenHistChunks = seqIterLengthsHistogram(reader.records(),bins=lenHistPred[1])
        reader.close()
        lenHistChunksNorm = lenHistChunks[0].astype(float)/numpy.sum(lenHistChunks[0])        
        print "Inp histogram: ", lenHistPredNorm, numpy.sum(lenHistPredNorm*lenHistPred[1])
        print "Out histogram: ", lenHistChunksNorm, numpy.sum(lenHistChunksNorm*lenHistChunks[1])
        print "Bins:          ", lenHistPred[1]
        inpHistArr = numpy.rec.fromarrays((lenHistPred[0],lenHistPred[1].astype(int)),names="count,bins")
        print "Inp counts:    ", inpHistArr, numpy.sum(inpHistArr['count'])

    def sampleTestByPredictLengthSvm(self,fastaFilePred):
        self.sampleTestByPredictLength(fastaFilePred)
        makeKmersTxt(self.testSeqSimPathVirFamily,self.testKmerTxtSimPathVirFamily,self.kmerSize,-1)
        kmersTxtToBin(self.testKmerTxtSimPathVirFamily)
        indices = loadObj(self.testIndPathVirFamily)
        kmerBinToSvmTxt(self.testSimPathVirFamily,self.testSvmSimPathVirFamily,indLabelTaxid=indices.rankTaxidTest)

    def makeViralSvm(self):
        #TMPEDIT: joinFastaByTaxid(self.knownSeqVirFasta,self.knownSeqVirKmerFasta)
        makeKmersTxt(self.knownSeqVirKmerFasta,self.kmerTxtPathVir,self.kmerSize,self.kmerWindow)
        kmersTxtToBin(self.kmerTxtPathVir)
        self.writeViralFamilyLabels()    
    
    def makeViralVsNonViralSvm(self):
        self.subSampleKmersNonViral()
        self.writeViralVsNonViralTrainingSet()
        kmerBinToSvmTxt(self.trainSetPathVirNonVir,self.trainSvmPathVirNonVir)
        kmerBinToSvmTxt(self.testSimPathVirFamily,self.testSvmPathVirNonVir,indLabelTaxid=numpy.ones(1000000,int))
    
    
    def makePredictSvm(self,rootPath):
        makeKmersTxt(rootPath+'.fas.gz',rootPath+'.6mers.gz',6,-1)
        kmerTxtToSvmTxt(rootPath+'.6mers.gz',rootPath+'.6mers.svm.inp')
        makeKmersTxt(rootPath+'.fas.gz',rootPath+'.7mers.gz',7,-1)
        kmerTxtToSvmTxt(rootPath+'.7mers.gz',rootPath+'.7mers.svm.inp')
    
    def predictViral(self,inpFasta,inpVnv,inpVfam):
        print inpFasta,inpVnv,inpVfam
        #self.loadTaxaTree()
        inp = open("taxonomy/taxdump/names.dmp",'r')
        taxaNames = {}
        for line in inp:
              flds = line.split('|')
              taxaNames[int(flds[0].strip())] = flds[1].strip()
        inp.close()
        iids = []
        lengths = []
        reader = fastaReaderGzip(inpFasta)
        for rec in reader.records():
            iid = int(rec.header().split(',')[1].split()[0])
            iids.append(iid)
            lengths.append(rec.seqLen())
        reader.close()
        seqIds = numpy.rec.fromarrays((numpy.array(iids),numpy.array(lengths)),names='iid,length')
        vnv = numpy.loadtxt(inpVnv,dtype=int)
        vfam = numpy.loadtxt(inpVfam,dtype=int)
        nSeq = len(seqIds[seqIds['length'] > 750])
        print "Input seq > 750:", nSeq
        vfamPred = vfam[numpy.logical_and(vnv == 1,seqIds['length'] > 750)]
        print "Predicted viral contigs:", len(vfamPred)
        vfamPredCnt = numpy.bincount(vfamPred)
        vfamPredCnt = whereItems(vfamPredCnt,vfamPredCnt>0)
        #self.taxaTree.getNode(item['ind']).lineageRanksStr()
        countNames = sorted([ (taxaNames[item['ind']],item['val']) for item in vfamPredCnt ],key=lambda item: item[1],reverse=True)
        out = open('pred.csv','w')
        for item in countNames:
            out.write("%s,%s,%d\n" % (item[0],item[1],float(item[1])/len(vfamPred)*100)) 
        out.write("Sequences classified,%s" % (nSeq,))
        out.close()

        
    
    def debug(self):
        pass
        #kmerTxt = KmerTxtReader(self.kmerTestFile)
        #kmerTxt.convertToBinary()
        #Pickling the taxa tree goes wrong apparently - memory consumption grows ~ 2x, and it takes longer than 
        #recreating from original flat file
        #self.taxaTree = objectDiskCacher(TaxaTree,os.path.basename(self.taxaNodesFile)+'.pkl')(ncbiDumpFile=self.taxaNodesFile)
        #kmerCounts = self.kmers.countValues(save=True,load=True)
        #print numpy.sum(kmerCounts)
        #return
        #rank_faker = RankFakerVisitor(topNode=viralRoot,lineage=viralRanksTemplate)
        #self.taxaTree.visitDepthTop(rank_faker,self.viralRootTaxid)
        #rank_reducer = RankReducerVisitor()
        #self.taxaTree.visitDepthTop(rank_reducer)
        

class SVMLibLinear:
    
    def __init__(self,workDir):
        self.workDir = workDir
        self.modelFileRel = "model.svm"
        makedir(workDir)
        self.binDir = os.environ['SVM_LIBLINEAR_BIN']
        self.binTrain = os.path.join(self.binDir,'train')
        self.binPredict = os.path.join(self.binDir,'predict')
        
    def train(self,trainFile):
        print "Starting training with %s samples of rank %s" % (iRec,len(rec.vals))        
        run(["svm_multiclass_learn","-c","1","-m","2000",trainFile,self.modelFileRel], cwd=self.workDir)
        print "Finished training"

        

class SVMLib:
    def __init__(self):
        import svm
        import cross_validation
        self.svm = svm
        self.cross_validation = cross_validation
        self.param = self.svm.svm_parameter(svm_type = self.svm.C_SVC, 
                                            kernel_type = self.svm.LINEAR, #self.svm.RBF,
                                            C = 2, gamma = 2, 
                                            cache_size=2000)
    
    def train(self,inpPath,testMode=False,testRecNum=10000):
        modelFile = inpPath + '.svm'
        kmersInp = KmerBinReader(inpPath)
        data = kmersInp.readAll()
        if testMode:
            data = numpy.concatenate((data[:testRecNum/2],data[-testRecNum/2:-1]))
        kmersInp.close()
        labels = data['taxid'].astype(numpy.float64)
        vals = data['vals'].astype(numpy.float64)
        print "labels: ", labels.shape, labels.dtype
        print "vectors: ", vals.shape, vals.dtype
        action = 'crossVal'
        if action == 'crossVal':
            print "Starting cross-validation with %s samples of rank %s" % (len(labels),len(vals[0]))
            self.cross_validation.do_cross_validation(vals,labels,self.param,5)
            print "Finished cross-validation"
        elif action == 'train':
            prob = self.svm.svm_problem(labels,vals)
            print "Starting training with %s samples of rank %s" % (len(labels),len(vals[0]))
            model = self.svm.svm_model(prob, self.param)
            print "Finished training"
            model.save(modelFile)
        else:
            raise ValueError("Unknown 'action' value: "+action)

class SVMMulticlass:
    
    def __init__(self,workDir='/export/tmp'):
        self.workDir = workDir
        self.sampleFileRel = "samples.svm"
        self.modelFileRel = "model.svm"
        makedir(workDir)
        
    def train(self,inpFile,testMode=False,testRecNum=10000):
        samples = open(os.path.join(self.workDir,self.sampleFileRel),'w', buffering=1024*1024)
        iRec = 0
        freqs = KmerFreqs(inpFile)    
        for rec in freqs.readValues():
            rec.write(samples)
            #samples.write(rec.asStr())
            samples.write("\n")
            iRec = iRec + 1
            if testMode and iRec >= testRecNum:
                break
        freqs.close()
        #important to flush, otherwise training program does not see the last records:
        samples.close()
        print "Starting training with %s samples of rank %s" % (iRec,len(rec.vals))        
        run(["svm_multiclass_learn","-c","1","-m","2000",self.sampleFileRel,self.modelFileRel], cwd=self.workDir)
        print "Finished training"



def buildViralRankLabels_bak(self,rank,indTaxid,testRatio,minKmerCountTrain):
    #maskVir = self.taxaTree.buildSubtreeMask(ids=(self.viralRootTaxid,), values=(1,))
    #maskRank[0:maxTaxid+1] ( any taxid -> taxid of 'rank' (e.g. family taxid) )
    maskRank = self.taxaTree.buildRankSubtreeMask(rank=rank,other=0)
    #maskSpe[0:maxTaxid+1]  (any taxid -> taxid of species)
    maskSpe = self.taxaTree.buildRankSubtreeMask(rank=exclusionRank,other=0)
    #kmerCnt[0:maxTaxid+1] (any taxid -> number of kmer records with this taxid (itself, not its subtree))
    ##kmerCnt = numpy.wrongResize(self.assignNodeCounts(),len(maskRank))
    #kmerCnt calculated on actual kmer file, so it might not have high taxids if all of them have zero counts,
    #so we extend it with zeros to match size of mask arrays
    kmerCnt = numpy.bincount(indTaxid)
    kmerCnt = numpy.wrongResize(kmerCnt,len(maskRank))
    #[0:maxRankTaxid+1] (rank taxid -> kmer count in subtree of this taxid)
    rankKmerCnt = numpy.bincount(maskRank,kmerCnt).astype(int)
    #[0:maxSpeTaxid+1] (species taxid -> kmer count in subtree of this taxid)
    speKmerCnt = numpy.bincount(maskSpe,kmerCnt).astype(int)
    #Build index as a dict of sets (rank Taxid -> (rank's species Taxid's))
    indRankToSpe = {}
    #we ignore all 'rank' that are not represented by any of the current kmers
    for (taxidRank,taxidSpe) in ( x for x in numpy.column_stack((maskRank,maskSpe)) 
                                  if x[0] != 0 and x[1] != 0 and rankKmerCnt[x[0]] > 0):
        try:
            indRankToSpe[taxidRank].add(taxidSpe)
        except KeyError:
            indRankToSpe[taxidRank] = set((taxidSpe,))
    rankSpeCnt = numpy.asarray([len(speSet) for speSet in indRankToSpe.itervalues()],numpy.int32)
    targetTestRankSpe = numpy.ceil(testRatio * numpy.median(rankSpeCnt))
    #from each family, randomly select 20% of its species for test set
    #but no less than 15 k-mers for training. If less than 15 k-mers at all, skip this
    #family
    #
    #[0:maxSpeTaxId] (species taxid -> 1 if selected for training set, 0 - otherwise
    speMaskTrain = numpy.zeros_like(speKmerCnt)
    #same as above but for testing set
    speMaskTest = numpy.zeros_like(speMaskTrain)
    ##[0:maxRankTaxid+1] ( rank taxid -> label (equal to rank taxid except than set to 0 for low kmer count )
    #rankToLabel = numpy.arange(len(rankKmerCnt))
    for (iRank,speSet) in indRankToSpe.iteritems():
        species = array(list(speSet),numpy.int32)
        spePermut = nrnd.permutation(species)
        spePermutKmerCnt = speKmerCnt[spePermut]
        spePermutKmerCntAccum = numpy.add.accumulate(spePermutKmerCnt)
        #we also consider kmers w/o species as valid 'rank' training sets,
        #and decrease the species defined minCount respectively
        minKmerRankSpe = minKmerCountTrain - (rankKmerCnt[iRank] - spePermutKmerCntAccum[-1])
        speIndMinCount = numpy.searchsorted(spePermutKmerCntAccum,minKmerRankSpe,side='right')
        nSpe = len(spePermut)
        nSpeTrain = nSpe - targetTestRankSpe
        nSpeTrain = max(nSpeTrain,speIndMinCount)
        nSpeTest = nSpe - nSpeTrain
        if spePermutKmerCntAccum[-1] < minKmerRankSpe:
            speTrain = spePermut[0:0]
        else:
            speTrain = spePermut[:nSpeTrain]
        speTest = spePermut[nSpeTrain:]
        speMaskTrain[speTrain] = 1
        speMaskTest[speTest] = 1
        if numpy.sum(spePermutKmerCnt) > 0:
            #pdb.set_trace()
            pass
        pass
    
    #[0:Total kmer count] (kmer ind -> taxid of species)
    indSpe = maskSpe[indTaxid]
    #same as above, but for 'rank'. This will be classification labels.
    indRank=maskRank[indTaxid]
    #[0:Total kmer count] (kmer ind -> kmer count for 'rank')
    indRankCnt = rankKmerCnt[indRank]
    #[0:Total kmer count] (kmer ind -> 1 if kmer in testing set, 0 - otherwise)
    indTest = speMaskTest[indSpe]
    #same as above but for training set: all that is not testing and whose 'rank' has kmer count not too small
    indTrain = numpy.logical_and(numpy.logical_and(indRankCnt >= minKmerCountTrain,numpy.logical_not(indTest)),
                                 indRank>0)        
    indTest = self.balanceMask(indRank, indTest)
    pdb.set_trace()
    labelData = Struct(indTrain=indTrain,indTest=indTest,indRank=indRank)
    return labelData