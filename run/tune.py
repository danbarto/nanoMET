"""
Tuning of MetSig
"""

# Standard imports
import ROOT
ROOT.gROOT.SetBatch(True)
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
argParser.add_argument('--selection',                     action='store',      type=str, default="diMuon-looseLeptonVeto-onZ" )
argParser.add_argument('--year',                          action='store',      type=int, default=2016, choices=[2016,2017,2018] )
argParser.add_argument('--runData',                       action='store_true', help='run data tuning', )
args = argParser.parse_args()

# Logger
import nanoMET.tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

# define setting
if args.year == 2016:
    postProcessing_directory = "2016_v21/dimuon/"
    trigger                  = ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL", "HLT_IsoMu24", "HLT_IsoTkMu24"]
    METPtVar                 = "MET_pt_nom"
    METPhiVar                = "MET_phi_nom"
    JetCollection            = "Jet_pt_nom"
    vetoEtaRegion            = (10.,10.)

    if args.runData:
        from nanoMET.samples.nanoTuples_Run2016_17Jul2018_postProcessed import *
        samples              = [DoubleMuon_Run2016]
        jer                  = "Summer16_25nsV1_DATA"
    else:
        from nanoMET.samples.nanoTuples_Summer16_postProcessed import *
        samples              = [DY_LO_16, Top_16, diboson_16, rare_16]
        jer                  = "Summer16_25nsV1_MC"

elif args.year == 2017:
    postProcessing_directory = "2017_v21/dimuon/"
    args.selection          += "-BadEEJetVeto"
    trigger                  = ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", "HLT_IsoMu27"]
    METPtVar                 = "METFixEE2017_pt"
    METPhiVar                = "METFixEE2017_phi"
    JetCollection            = "Jet_pt_nom"
    vetoEtaRegion            = (2.65,3.14)

    if args.runData:
        from nanoMET.samples.nanoTuples_Run2017_31Mar2018_postProcessed import *
        samples              = [DoubleMuon_Run2017BCDE]
        jer                  = "Fall17_V3_DATA"
    else:
        from nanoMET.samples.nanoTuples_Fall17_postProcessed import *
        samples              = [DY_LO_17, Top_17, diboson_17, rare_17]
        jer                  = "Fall17_V3_MC"

elif args.year == 2018:
    postProcessing_directory = "2018_v21/dimuon/"
    trigger                  = ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", "HLT_IsoMu24"]
    METPtVar                 = "MET_pt_nom"
    METPhiVar                = "MET_phi_nom"
    JetCollection            = "Jet_pt_nom"
    vetoEtaRegion            = (10.,10.)

    if args.runData:
        from nanoMET.samples.nanoTuples_Run2018_17Sep2018_postProcessed import *
        samples              = [DoubleMuon_Run2018]
        jer                  = "Autumn18_V1_DATA"
    else:
        from nanoMET.samples.nanoTuples_Autumn18_postProcessed import *
        samples              = [DY_LO_18, Top_18, diboson_18, rare_18]
        jer                  = "Autumn18_V1_MC"


# calculate setting
preselection    = cutInterpreter.cutString(args.selection)
triggerSel      = "(%s)"%"||".join(["Alt$(%s,0)"%trigg for trigg in trigger])
eventfilter     = getFilterCut( args.year, isData=args.runData )
sel             = "&&".join([preselection, triggerSel, eventfilter])
JR              = JetResolution(jer)
version         = postProcessing_directory.split("/")[0]
outfile         = "results/tune_%s_%i_%s_%s_sumPt%i_max%i_%s"%("Data" if args.runData else "MC", args.year, jer, args.selection, args.jetThreshold, args.maxSig, version)
if args.pTdependent:
    outfile    += "_pTdep"

# run
r = run(samples, sel, JR, outfile=outfile, maxN=1e5, METPtVar=METPtVar, METPhiVar=METPhiVar, JetCollection=JetCollection, vetoEtaRegion=vetoEtaRegion, jetThreshold=args.jetThreshold, puWeight="puWeight", ttbarModifier=args.ttbarModifier, pTdepMetSig=args.pTdependent)
LL = r.getLL(r.defaultStart)
r.minimize(start=r.defaultStart, maxSig=args.maxSig)

