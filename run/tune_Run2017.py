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

postProcessing_directory = "2017_v8/dimuon/"
from nanoMET.samples.nanoTuples_Run2017_31Mar2018_postProcessed import *

# define the selection
leptonSelection = "Sum$(Muon_pt>20&&Muon_isGoodMuon)==2&&Sum$(Muon_pt>35&&Muon_isGoodMuon)>0"
preselection    = cutInterpreter.cutString('looseLeptonVeto-onZ')
trigger         = "( %s )"%" || ".join(['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ', "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", 'HLT_IsoMu27'])
EE_protection   = "Sum$((cos(Jet_phi-MET_phi)*Jet_pt*Jet_neEmEF)*(cos(Jet_phi-MET_phi)<cos(2*pi/3.)))/MET_pt"
eventfilter     = getFilterCut( 2017, isData=True)

sel             = " && ".join([leptonSelection, preselection, trigger, eventfilter])

JR = JetResolution('Fall17_V3_DATA') # similar to Fall17_25nsV1 and Fall17_V2

## only run over max 1M event per sample, uncertainty is anyway low. Need to confirm that the parameters really converged then.
#DoubleMuon_Run2016.reduceFiles(to=3)
r = run([DoubleMuon_Run2017], sel, JR, outfile="results/tune_DoubleMuon_Run2017_FixEE2017_noMax_sumPt25_v9", METPtVar="METFixEE2017_pt", METPhiVar="METFixEE2017_phi", maxN=5e5, vetoEtaRegion=(10,10), jetThreshold=25.)

LL = r.getLL( [1.0, 1.0, 1.0, 1.0, 1.0, 0., .5] )

r.minimize(maxSig=80)

