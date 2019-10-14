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

from Samples.Tools.metFilters   import getFilterCut

from run import run

postProcessing_directory = "2016_v6/dimuon/"
from nanoMET.samples.nanoTuples_Summer16_postProcessed import *

# define the selection
preselection    = "Sum$(Jet_pt>30&&Jet_jetId&&abs(Jet_eta)<2.4)>=0 && Sum$(Muon_pt>25&&Muon_isGoodMuon)==2 && Sum$(Electron_pt>10&&abs(Electron_eta)<2.5&&Electron_cutBased>0&&abs(Electron_pfRelIso03_all)<0.4)==0 && abs(dl_mass-91.2)<10"
trigger         = "( %s )"%" || ".join(['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL', 'HLT_IsoMu24', 'HLT_IsoTkMu24'])
eventfilter     = getFilterCut( 2016, isData=False)
sel             = " && ".join([preselection,trigger,eventfilter])

JR = JetResolution('Summer16_25nsV1_MC') # Spring16_25nsV6_MC, Summer16_25nsV1_MC

## only run over max 1M event per sample, uncertainty is anyway low. Need to confirm that the parameters really converged then.
samples = [DY_LO_16, Top_16, VVTo2L2Nu_16, WJets_16]
samples = [WJets_16]
r = run(samples, sel, JR, outfile="results/tune_Summer16_incl_puWeight_sumPt15_max25_v6", METPtVar="MET_pt_nom", METPhiVar="MET_phi_nom", JetCollection="Jet_pt_nom", maxN=3*1e5, puWeight='puWeight', jetThreshold=15.)

start=[1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 0., .5]

LL = r.getLL( start )

r.minimize(start=start, maxSig=25)

