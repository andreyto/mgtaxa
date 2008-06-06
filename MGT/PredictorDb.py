"""Describe database of predictors for taxonomic tree."""

from MGT.Common import *
from MGT.Svm import *

class Predictor(MGTOptions):


    def __init__(self,predDb,node,id):
        MGTOptions.__init__(self)
        self.predDb = predDb
        if id is None:
            predDb.db.ddl("""
            insert into %s
            (id_node, pred_type, state)
            values
            (%s, '%s', '%s')
            """ % (predDb.predictorTable,
                   node.id,
                   predDb.predType(self),
                   predDb.predState.empty))
            id = predDb.db.selectScalar("SELECT LAST_INSERT_ID()")
        self.id = id
        self.workDir = os.path.join(predDb.predictorDir,"%07i" % id)
        makedir(self.workDir)
        self.trainFile = os.path.join(self.workDir,'train.svm')

    def trainingWriter(self):
        return SvmSparseFeatureWriterTxt(self.trainFile)


class MClassPredictor(Predictor):
    
    def __init__(self,**kw):
        Predictor.__init__(self,**kw)


class PredictorDb(MGTOptions):

    predState = Struct(empty=0,trainSubmitted=10)

    predTypeToClass = { 'ms' : MClassPredictor }
    predClassToType = dict( [ (value, key) for (key, value) in predTypeToClass.items() ] )

    def __init__(self,db,reset=False):
        MGTOptions.__init__(self)
        self.db = db
        self.workDir = self.predictorDir
        if reset:
            rmrf(self.workDir)
            dropList = ["table %s" % self.predictorTable]
            ignoreError = False
        else:
            dropList = []
            ignoreError = True
        makedir(self.workDir)
        db.ddl("""
        create table %s (id integer not null auto_increment,
                         id_node integer not null,
                         pred_type char(2), 
                         state smallint,
                         primary key (id),
                         key id_node(id_node))
        """ % self.predictorTable,
        dropList=dropList,
        ignoreError=ignoreError)

    def predType(self,predInst):
        return self.predClassToType[predInst.__class__]

    def predClass(self,predType):
        return self.predTypeToClass[predType]

    def newPredictor(self,node,predType,**kw):
        return self.predClass(predType)(predDb=self,node=node,id=None,**kw)
