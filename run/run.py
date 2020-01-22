"""

"""

# Standard imports
import ROOT
import os
import sys
import math
import copy
import itertools
import json
import random

from RootTools.core.standard    import *

from nanoMET.core.Event         import Event
from nanoMET.core.Likelihood    import minLL
from nanoMET.core.JetResolution import JetResolution
from nanoMET.tools.progressbar  import update_progress

#
# Logger
#
import nanoMET.tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   "INFO", logFile = None)
logger_rt = logger_rt.get_logger("INFO", logFile = None)


def getNVtxWeight( hist, npv ):
    return hist.GetBinContent(hist.FindBin(npv))

class run:

    def __init__(self, samples, selection, jetResolution, outfile="results/tune", METPtVar="MET_pt", METPhiVar="MET_phi", JetCollection="Jet_pt", maxN=1e6, vetoEtaRegion=(10,10), jetThreshold=15., puWeight="puWeight", ttbarModifier=1, pTdepMetSig=True, nvtxHist=None):
        # Need fill a list in order to do the minimization, reading from the tree is too slow

        self.eventlist = []
        if samples[0].isData:
            self.variables = map( TreeVariable.fromString,  ["nJet/I", "fixedGridRhoFastjetAll/F", "%s/F"%METPtVar, "%s/F"%METPhiVar, "MET_sumPt/F"] )
        else:
            self.variables = map( TreeVariable.fromString,  ["weight/F", "puWeight/F", "puWeightUp/F", "puWeightDown/F", "nJet/I", "fixedGridRhoFastjetAll/F", "PV_npvsGood/I", "%s/F"%METPtVar, "%s/F"%METPhiVar, "MET_sumPt/F"] )
        self.variables += [VectorTreeVariable.fromString("Jet[pt/F,eta/F,phi/F,cleanmask/O,jetId/I,cleanmaskMETSig/I,cleanmaskMETSigRec/I]" ) ]
        if JetCollection=="Jet_pt_nom":
            self.variables += [VectorTreeVariable.fromString("Jet[pt_nom/F]")]
        elif JetCollection=="Jet_pt_jerUp":
            self.variables += [VectorTreeVariable.fromString("Jet[pt_jerUp/F]")]

        self.outfile = os.path.abspath( outfile )
        self.outdir  = os.path.dirname( self.outfile )

        if not os.path.isdir( self.outdir ):
            os.makedirs( self.outdir )        

        self.pTdepMetSig = pTdepMetSig
        self.paramRange  = 12 if self.pTdepMetSig else 7

        self.variableNames = ["a%i"%i for i in range(self.paramRange-2)] + ["u1","u2"]
        self.defaultStep   = [0.05]*self.paramRange
        self.defaultStart  = [1.5]*(self.paramRange-2) + [0., .5]

        for s in samples:
            logger.info( "Now working on sample: %s"%s.name )
            s.setSelectionString(selection)
            logger.info( "Getting number of events" )
            nEvents = s.getYieldFromDraw(selectionString="(1)", weightString="(1)")
            logger.info( "Running over %s events for sample %s."%(nEvents["val"], s.name) )
            nEvents = nEvents["val"]
            if nEvents > maxN:
                weightModifier = nEvents/float(maxN)
                fracToKeep = float(maxN)/nEvents
                logger.info( "Preliminary factor for modifying the event weight: %s"%weightModifier )
            else:
                weightModifier = 1
                fracToKeep = 1
            
            if "top" in s.name.lower():
                weightModifier *= ttbarModifier
            logger.info( "Filling the eventlist" )
            reader = s.treeReader(variables=self.variables)
            reader.start()
            i = 0
            tmp_eventlist = []
            while reader.run():
                i+=1
                nvtxWeight = 1
                if random.random() < fracToKeep:
                    if nvtxHist:
                        nvtxWeight = getNVtxWeight(nvtxHist, reader.event.PV_npvsGood)
                    tmp_eventlist += [Event(reader.event, jetResolution, weightModifier=weightModifier*nvtxWeight, isData=s.isData, METPtVar=METPtVar, METPhiVar=METPhiVar, JetCollection=JetCollection, vetoEtaRegion=vetoEtaRegion, jetThreshold=jetThreshold, puWeight=puWeight, pTdepMetSig=self.pTdepMetSig)]
                update_progress(i/nEvents)
            
            self.eventlist += tmp_eventlist


    def getLL(self, args):
        # recalculate MET Significance, Determinant and LL for given parameters
        # although it uses C methods it"s still rather slow
        LLs = map(lambda x: x.calcLL(args), self.eventlist)
        LL = sum(LLs)
        
        return LL

    def minimize(self, start=None, step=None, maxSig=25):
        gmin = ROOT.Math.Factory.CreateMinimizer("Minuit2")
        gmin.SetTolerance(10.0)
        gmin.SetStrategy(0)
        gmin.SetPrintLevel(3)
        
        step     = step  if step  else self.defaultStep
        start    = start if start else self.defaultStart

        LL = minLL(self, pTdepMetSig=self.pTdepMetSig)

        gmin.SetFunction(LL)

        logger.info( "Minimizing parameters %s"%self.variableNames )
        logger.info( "With stepsize of %s"%step )
        logger.info( "And starting values %s"%start )
    
        for i in range(self.paramRange):
            gmin.SetVariable( i, self.variableNames[i], start[i], step[i] )
    
        gmin.Minimize()

        #filter events with high significance
        self.eventlist = filter( lambda x: x.MET_sig < maxSig, self.eventlist )
        
        logger.info( "Now fitting after applying significance cut" )
        logger.info( "Total events: %i"%len(self.eventlist) )
        
        gmin.SetStrategy(1)
        gmin.Minimize()
        gmin.Hesse()
        
        pars = [ gmin.X()[i]      for i in range(self.paramRange) ]
        uncs = [ gmin.Errors()[i] for i in range(self.paramRange) ]
#        pars = map( float, gmin.X() )
#        uncs = map( float, gmin.Errors() )
        
        with open(self.outfile+".txt", "w") as of:
            json.dump(pars, of)
        with open(self.outfile+"_unc.txt", "w") as of:
            json.dump(uncs, of)
    



if __name__ == "__main__":

    from nanoMET.core.JetResolution   import JetResolution
    from nanoMET.tools.cutInterpreter import cutInterpreter
    from nanoMET.tools.metFilters     import getFilterCut

    # Logger
    import nanoMET.tools.logger as logger
    import RootTools.core.logger as logger_rt
    logger    = logger.get_logger(   "INFO", logFile = None)
    logger_rt = logger_rt.get_logger("INFO", logFile = None)

    postProcessing_directory = "2016_v22/dimuon/"
    from nanoMET.samples.nanoTuples_Summer16_postProcessed import *
    from nanoMET.samples.nanoTuples_Run2016_17Jul2018_postProcessed import *

    # define setting
    selection       = "diMuon-looseLeptonVeto-onZ"
    samples         = [DY_LO_16, Top_16, diboson_16, rare_16]
    for s in samples:
        s.reduceFiles( to=5 )
    trigger         = ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL", "HLT_IsoMu24", "HLT_IsoTkMu24"]
    jer             = "Summer16_25nsV1_MC"
    METPtVar        = "MET_pt_nom"
    METPhiVar       = "MET_phi_nom"
    JetCollection   = "Jet_pt_nom"
    vetoEtaRegion   = (10.,10.)
    minimize        = True
    ttbarModifier   = 1
    maxSig          = 25
    jetThreshold    = 15
    pTdepMetSig     = False

    # calculate setting
    preselection    = cutInterpreter.cutString(selection)
    triggerSel      = "(%s)"%"||".join(["Alt$(%s,0)"%trigg for trigg in trigger])
    eventfilter     = getFilterCut( 2016, isData=False )
    sel             = "&&".join([preselection, triggerSel, eventfilter])
    JR              = JetResolution(jer)
    version         = postProcessing_directory.split("/")[0]
    outfile         = "results/test_tune_%s_%s_puWeight_sumPt%i_max%i_ttbar%i_%s"%(jer,selection,jetThreshold,maxSig,ttbarModifier,version)

    r = run(samples, sel, JR, outfile=outfile, maxN=3e5, METPtVar=METPtVar, METPhiVar=METPhiVar, JetCollection=JetCollection, vetoEtaRegion=vetoEtaRegion, jetThreshold=jetThreshold, puWeight="puWeight", ttbarModifier=ttbarModifier, pTdepMetSig=pTdepMetSig)
    LL = r.getLL(r.defaultStart)
    if minimize:
        r.minimize(start=r.defaultStart, maxSig=maxSig)

