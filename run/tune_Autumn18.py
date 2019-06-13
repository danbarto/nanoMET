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

postProcessing_directory = "2018_v12/dimuon/"
from nanoMET.samples.nanoTuples_Autumn18_postProcessed import *
#data_directory = "/afs/hephy.at/data/cms01/nanoTuples/"
#postProcessing_directory = "stops_2018_nano_v0p10/dilep/"
#from StopsDilepton.samples.nanoTuples_Autumn18_postProcessed import *

# define the selection
leptonSelection = "Sum$(Muon_pt>20&&Muon_isGoodMuon)==2&&Sum$(Muon_pt>35&&Muon_isGoodMuon)>0"
preselection    = cutInterpreter.cutString('looseLeptonVeto-onZ')
trigger         = "( %s )"%" || ".join(['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8', 'HLT_IsoMu24'])
eventfilter     = getFilterCut( 2018, isData=False)

sel             = " && ".join([leptonSelection, preselection, trigger, eventfilter])

JR = JetResolution('Autumn18_V1_MC')

## only run over max 1M event per sample, uncertainty is anyway low. Need to confirm that the parameters really converged then.
samples = [DY_LO_18, Top_18, diboson_18, rare_18]
#samples = [DY_LO_18, Top_pow_18, multiBoson_18, TTX_18]
r = run(samples, sel, JR, outfile="results/tune_Autumn18_incl_1p5_puWeight_sumPt25_max25_ttbar1_v12", maxN=3*1e5, METPtVar="MET_pt_nom", METPhiVar="MET_phi_nom", JetCollection="Jet_pt_nom", jetThreshold=25., puWeight='puWeight', ttbarModifier=1.)

LL = r.getLL( [1.5, 1.5, 1.5, 1.5, 1.5, 0., .5] )

start=[1.5, 1.5, 1.5, 1.5, 1.5, 0., .5]

r.minimize(start=start, maxSig=25)

