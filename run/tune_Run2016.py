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
from nanoMET.samples.nanoTuples_Run2016_17Jul2018_postProcessed import *

# define the selection
preselection    = "Sum$(Jet_pt>30&&Jet_jetId&&abs(Jet_eta)<2.4)>=0 && Sum$(Muon_pt>25&&Muon_isGoodMuon)==2 && Sum$(Electron_pt>10&&abs(Electron_eta)<2.5&&Electron_cutBased>0&&abs(Electron_pfRelIso03_all)<0.4)==0 && abs(dl_mass-91.2)<10"
trigger         = "( %s )"%" || ".join(['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL', 'HLT_IsoMu24', 'HLT_IsoTkMu24'])
eventfilter     = getFilterCut( 2016, isData=True).replace('&&weight>0','')
sel             = " && ".join([preselection,trigger,eventfilter])

JR = JetResolution('Summer16_25nsV1_DATA')

## only run over max 1M event per sample, uncertainty is anyway low. Need to confirm that the parameters really converged then.
#DoubleMuon_Run2016.reduceFiles(to=3)
r = run([DoubleMuon_Run2016], sel, JR, outfile="results/tune_DoubleMuon_17Jul2018_incl_sumPt15_max25_v4", maxN=3*1e5, jetThreshold=15.)

LL = r.getLL( [1.5, 1.5, 1.5, 1.5, 1.5, 0., .5] )

start=[1.5, 1.5, 1.5, 1.5, 1.5, 0., .5]

r.minimize(start=start, maxSig=25)

