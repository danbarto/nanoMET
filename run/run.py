'''

'''

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
from nanoMET.tools.progressbar  import update_progress
from nanoMET.core.JetResolution import JetResolution


class run:

    def __init__(self, samples, selection, jetResolution, outfile="results/tune", METCollection='MET_pt', JetCollection="Jet_pt", maxN=1e6, vetoEtaRegion=(10,10)):
        # Need fill a list in order to do the minimization, reading from the tree is too slow

        self.eventlist = []
        if samples[0].isData:
            self.variables = map( TreeVariable.fromString,  ['nJet/I', 'fixedGridRhoFastjetAll/F', '%s/F'%METCollection, 'MET_phi/F', 'MET_sumPt/F'] )
        else:
            self.variables = map( TreeVariable.fromString,  ['weight/F', 'nJet/I', 'fixedGridRhoFastjetAll/F', '%s/F'%METCollection, 'MET_phi/F', 'MET_sumPt/F'] )
        self.variables += [VectorTreeVariable.fromString('Jet[pt/F,eta/F,phi/F,cleanmask/O,jetId/I,cleanmaskMETSig/I,cleanmaskMETSigRec/I]' ) ]
        if JetCollection=='Jet_pt_nom':
            self.variables += [VectorTreeVariable.fromString('Jet[pt_nom/F]')]
        self.outfile = outfile

        for s in samples:
            print
            print "Now working on sample: %s"%s.name
            s.setSelectionString(selection)
            print "Getting number of events"
            nEvents = s.getYieldFromDraw(selectionString='(1)', weightString='(1)')
            print "Running over %s events for sample %s."%(nEvents['val'], s.name)
            nEvents = nEvents['val']
            if nEvents > maxN:
                weightModifier = nEvents/float(maxN)
                fracToKeep = float(maxN)/nEvents
                print "Preliminary factor for modifying the event weight: %s"%weightModifier
            else:
                weightModifier = 1
                fracToKeep = 1
            
            print "Filling the eventlist"
            reader = s.treeReader(variables=self.variables)
            reader.start()
            i = 0
            tmp_eventlist = []
            while reader.run():
                i+=1
                if random.random() < fracToKeep:
                    tmp_eventlist += [Event(reader.event, jetResolution, weightModifier=weightModifier, isData=s.isData, METCollection=METCollection, JetCollection=JetCollection, vetoEtaRegion=(10,10))]
                update_progress(i/nEvents)
            
            print

            self.eventlist += tmp_eventlist


    def getLL(self, args):
        # recalculate MET Significance, Determinant and LL for given parameters
        # although it uses C methods it's still rather slow
        LLs = map(lambda x: x.calcLL(args), self.eventlist)
        LL = sum(LLs)
        
        return LL


    def minimize(self, start=[1.0, 1.0, 1.0, 1.0, 1.0, 0., .5], step=[0.05]*7, maxSig = 9):
        gmin = ROOT.Math.Factory.CreateMinimizer("Minuit2")
        gmin.SetTolerance(10.0)
        gmin.SetStrategy(0)
        gmin.SetPrintLevel(3)
        
        variable  = ['a1','a2','a3','a4','a5','u1','u2']
        step  = step
        start = start

        LL = minLL(self)

        gmin.SetFunction(LL)

        print 'Minimizing parameters',variable
        print 'With stepsize of', step
        print 'And starting values', start
    
        for i in range(0,7):
            gmin.SetVariable(i,variable[i],start[i], step[i])
    
        gmin.Minimize()

        #filter events with high significance
        self.eventlist = [x  for x in self.eventlist if x.MET_sig < maxSig ]
        
        print
        print 'Now fitting after applying significance cut'
        print 'Total events:',len(self.eventlist)
        
        gmin.SetStrategy(1)
        gmin.Minimize()
        gmin.Hesse()
        
        pars = [gmin.X()[i] for i in range(0,7)]
        uncs = [gmin.Errors()[i] for i in range(0,7)]
        
        with open(self.outfile+'.txt', 'w') as of:
            json.dump(pars, of)
        with open(self.outfile+'_unc.txt', 'w') as of:
            json.dump(uncs, of)
    



if __name__ == '__main__':

    postProcessing_directory = "2016_v6/dimuon/"
    from nanoMET.samples.nanoTuples_Summer16_postProcessed import *
    from nanoMET.samples.nanoTuples_Run2016_17Jul2018_postProcessed import *

    # define the selection
    preselection    = "jsonPassed && Sum$(Jet_pt>30&&Jet_jetId&&abs(Jet_eta)<2.4)>=0 && Sum$(Muon_pt>25&&Muon_isGoodMuon)==2 && Sum$(Electron_pt>10&&abs(Electron_eta)<2.5&&Electron_cutBased>0&&abs(Electron_pfRelIso03_all)<0.4)==0 && abs(dl_mass-91.2)<10"
    trigger         = " && ".join(['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL', 'HLT_IsoMu24', 'HLT_IsoTkMu24'])
    eventfilter     = " && ".join(['Flag_METFilters','Flag_goodVertices','Flag_HBHENoiseIsoFilter','Flag_HBHENoiseFilter','Flag_globalTightHalo2016Filter','Flag_EcalDeadCellTriggerPrimitiveFilter'])
    sel             = " && ".join([preselection,trigger,eventfilter])

    JR = JetResolution('Summer16_25nsV1_MC')

    #DoubleMuon_Run2016.reduceFiles(to=1)

    ## only run over max 1M event per sample, uncertainty is anyway low. Need to confirm that the parameters really converged then.
    r = run([DY_LO_16], sel, JR, outfile="results/tune_Summer16_small", maxN=1e4)

    LL = r.getLL( [1.0, 1.0, 1.0, 1.0, 1.0, 0., .5] )

    minimize = True
    if minimize:

        r.minimize()
