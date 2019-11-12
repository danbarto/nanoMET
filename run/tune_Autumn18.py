"""
Tuning of Autumn18
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

from RootTools.core.standard      import *

from nanoMET.core.JetResolution   import JetResolution
from nanoMET.tools.cutInterpreter import cutInterpreter
from nanoMET.tools.metFilters     import getFilterCut

from run import run

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',                      action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'])
argParser.add_argument('--pTdependent',                   action='store_true', help='run pT dependent METSig tuning', )
argParser.add_argument('--maxSig',                        action='store',      type=int, default=25 )
argParser.add_argument('--jetThreshold',                  action='store',      type=int, default=15 )
argParser.add_argument('--ttbarModifier',                 action='store',      type=int, default=1 )
args = argParser.parse_args()

# Logger
import nanoMET.tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

postProcessing_directory = "2018_v21/dimuon/"
from nanoMET.samples.nanoTuples_Autumn18_postProcessed import *

# define setting
selection       = "diMuon-looseLeptonVeto-onZ"
samples         = [DY_LO_18, Top_18, diboson_18, rare_18]
trigger         = ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", "HLT_IsoMu24"]
jer             = "Autumn18_V1_MC"
METPtVar        = "MET_pt_nom"
METPhiVar       = "MET_phi_nom"
JetCollection   = "Jet_pt_nom"
vetoEtaRegion   = (10.,10.)

# calculate setting
preselection    = cutInterpreter.cutString(selection)
triggerSel      = "(%s)"%"||".join(["Alt$(%s,0)"%trigg for trigg in trigger])
eventfilter     = getFilterCut( 2018, isData=False )
sel             = "&&".join([preselection, triggerSel, eventfilter])
JR              = JetResolution(jer)
version         = postProcessing_directory.split("/")[0]
outfile         = "results/tune_%s_%s_puWeight_sumPt%i_max%i_ttbar%i_%s"%(jer,selection,args.jetThreshold,args.maxSig,args.ttbarModifier,version)

r = run(samples, sel, JR, outfile=outfile, maxN=3e5, METPtVar=METPtVar, METPhiVar=METPhiVar, JetCollection=JetCollection, vetoEtaRegion=vetoEtaRegion, jetThreshold=args.jetThreshold, puWeight="puWeight", ttbarModifier=args.ttbarModifier, pTdepMetSig=args.pTdependent)
LL = r.getLL(r.defaultStart)
r.minimize(start=r.defaultStart, maxSig=args.maxSig)

