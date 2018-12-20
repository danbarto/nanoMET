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

from run import run


from nanoMET.samples.nanoTuples_Run2016_17Jul2018_postProcessed import *

# define the selection
preselection    = "jsonPassed && Sum$(Jet_pt>30&&Jet_jetId&&abs(Jet_eta)<2.4)>=0 && Sum$(Muon_pt>25&&Muon_isGoodMuon)==2 && Sum$(Electron_pt>10&&abs(Electron_eta)<2.5&&Electron_cutBased>0&&abs(Electron_pfRelIso03_all)<0.4)==0 && abs(dl_mass-91.2)<10"
trigger         = "( %s )"%" || ".join(['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL', 'HLT_IsoMu24', 'HLT_IsoTkMu24'])
eventfilter     = " && ".join(['Flag_METFilters','Flag_goodVertices','Flag_HBHENoiseIsoFilter','Flag_HBHENoiseFilter','Flag_globalTightHalo2016Filter','Flag_EcalDeadCellTriggerPrimitiveFilter'])
sel             = " && ".join([preselection,trigger,eventfilter])

JR = JetResolution('Summer16_25nsV1_DATA')

## only run over max 1M event per sample, uncertainty is anyway low. Need to confirm that the parameters really converged then.
r = run([DoubleMuon_Run2016], sel, JR, outfile="results/tune_DoubleMuon_test2", maxN=1e5)

LL = r.getLL( [1.0, 1.0, 1.0, 1.0, 1.0, 0., .5] )

r.minimize()

