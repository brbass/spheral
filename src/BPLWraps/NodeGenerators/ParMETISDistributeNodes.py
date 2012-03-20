import Spheral
import distributeNodesGeneric

#-------------------------------------------------------------------------------
# Domain decompose using ParMETIS (2d method).
#-------------------------------------------------------------------------------
def distributeNodes2d(*listOfNodeTuples):
    distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                  Spheral.DataBase2d,
                                                  Spheral.globalNodeIDs,
                                                  Spheral.ParmetisRedistributeNodes2d)

#-------------------------------------------------------------------------------
# Domain decompose using ParMETIS (3d method).
#-------------------------------------------------------------------------------
def distributeNodes3d(*listOfNodeTuples):
    distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                  Spheral.DataBase3d,
                                                  Spheral.globalNodeIDs,
                                                  Spheral.ParmetisRedistributeNodes3d)

