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
import pickle

from RootTools.core.standard    import *

from nanoMET.core.JetResolution import JetResolution
from nanoMET.tools.cutInterpreter import cutInterpreter

from Samples.Tools.metFilters   import getFilterCut

from run import run

postProcessing_directory = "2017_v11/dimuon/"
from nanoMET.samples.nanoTuples_Fall17_postProcessed import *

# define the selection
leptonSelection = "Sum$(Muon_pt>20&&Muon_isGoodMuon)==2&&Sum$(Muon_pt>35&&Muon_isGoodMuon)>0"
preselection    = cutInterpreter.cutString('looseLeptonVeto-onZ')
trigger         = "( %s )"%" || ".join(['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ', "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", 'HLT_IsoMu27'])
EEveto          = "Sum$((2.6<abs(Jet_eta)&&abs(Jet_eta)<3&&Jet_pt>30))==0"
eventfilter     = getFilterCut( 2017, isData=False)

sel             = " && ".join([leptonSelection, preselection, trigger, eventfilter, EEveto])

JR = JetResolution('Fall17_V3_MC')

## only run over max 1M event per sample, uncertainty is anyway low. Need to confirm that the parameters really converged then.
samples = [DY_LO_17, Top_17, VVTo2L2Nu_17]
nvtxHist = pickle.load(file('nvtx_hists.pkl', 'r'))['dataBCDE']
r = run(samples, sel, JR, outfile="results/tune_Fall17_FixEE2017_1p5_sumPt15_max9_EEveto_allMC_vetoEta_v11", METPtVar="METFixEE2017_pt", METPhiVar="METFixEE2017_phi", JetCollection="Jet_pt_nom", maxN=3e5, vetoEtaRegion=(2.65,3.14), jetThreshold=15., puWeight='puWeight', ttbarModifier=1., nvtxHist=None)

start=[1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 0., .5]

LL = r.getLL( start )

r.minimize(start=start, maxSig=9)

