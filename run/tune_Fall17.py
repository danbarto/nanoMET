'''
Tuning of Run2016
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

from nanoMET.core.JetResolution import JetResolution
from nanoMET.tools.cutInterpreter import cutInterpreter

from Samples.Tools.metFilters   import getFilterCut

from run import run

postProcessing_directory = "2017_v11/dimuon/"
from nanoMET.samples.nanoTuples_Fall17_postProcessed import *

# define the selection
leptonSelection = "Sum$(Muon_pt>20&&Muon_isGoodMuon)==2&&Sum$(Muon_pt>35&&Muon_isGoodMuon)>0"
preselection    = cutInterpreter.cutString('looseLeptonVeto-offZ-mll50-nCleanJet1p')
#preselection    = cutInterpreter.cutString('looseLeptonVeto-offZ-mll50')
trigger         = "( %s )"%" || ".join(['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ', "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", 'HLT_IsoMu27'])
#EE_protection   = "Sum$((cos(Jet_phi-MET_phi)*Jet_pt*Jet_neEmEF)*(cos(Jet_phi-MET_phi)<cos(2*pi/3.)))/MET_pt"
#EE_jet_veto     = "Sum$(abs(Jet_eta)>2.5&&abs(Jet_eta)<3.2)==0"
EEveto          = "Sum$((2.6<abs(Jet_eta)&&abs(Jet_eta)<3&&Jet_pt>30))==0"
eventfilter     = getFilterCut( 2017, isData=False)
#genMet          = "GenMET_pt<10"

sel             = " && ".join([leptonSelection, preselection, trigger, eventfilter, EEveto])

JR = JetResolution('Fall17_V3_MC')

## only run over max 1M event per sample, uncertainty is anyway low. Need to confirm that the parameters really converged then.
#DoubleMuon_Run2016.reduceFiles(to=3)
#r = run([VVTo2L2Nu_17], sel, JR, outfile="results/tune_Fall17_test", maxN=1e5)
#samples = [DY_LO_17, Top_17, VVTo2L2Nu_17, WJets_17]
samples = [DY_LO_17, Top_17, VVTo2L2Nu_17]
#samples = [DY_LO_17]
r = run(samples, sel, JR, outfile="results/tune_Fall17_FixEE2017_1p5_puWeight_sumPt25_max25_EEveto_ttbar5_allMC_v11", METPtVar="METFixEE2017_pt", METPhiVar="METFixEE2017_phi", maxN=3e5, vetoEtaRegion=(10.,10.), jetThreshold=25., puWeight="puWeight", ttbarModifier=5.)

LL = r.getLL( [1.5, 1.5, 1.5, 1.5, 1.5, 0., .5] )

start=[1.5, 1.5, 1.5, 1.5, 1.5, 0., .5]

r.minimize(start=start, maxSig=25)

