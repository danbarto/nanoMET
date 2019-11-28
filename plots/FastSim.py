''' Analysis script for 1D 2l plots (RootTools)
'''

# Standard imports
import ROOT
import os
import copy
import itertools
from math import sqrt, cos, sin, pi, acos

# RootTools
from RootTools.core.standard import *

import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',  action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',     action='store_true', help='Small?')
args = argParser.parse_args()

data_directory              = "/afs/hephy.at/data/dspitzbart03/nanoAOD/"
postProcessing_directory    = "dimuon/"
plot_directory              = "/afs/hephy.at/user/d/dspitzbart/www/nanoMET/FastSim/"

# Logger
import nanoMET.tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

dirs = {}
dirs['DYJetsToLL_M50']      = ["DYJetsToLL_M50_MLM_S16_80X_priv"]
dirs['DYJetsToLL_M50_FS']   = ['DYJetsToLL_M50_MLM_FS_S16_80X_priv']
dirs['TTJets']              = ['TTJets_MLM_S16_80X_priv']
dirs['TTJets_FS']           = ['TTJets_MLM_FS_S16_80X_priv']
directories = { key : [ os.path.join( data_directory, postProcessing_directory, dir) for dir in dirs[key]] for key in dirs.keys()}

DYJetsToLL_M50      = Sample.fromDirectory(name="DY", treeName="Events", isData=False, color=1, texName="DY FullSim", directory=directories['DYJetsToLL_M50'])
DYJetsToLL_M50_FS   = Sample.fromDirectory(name="DY_FS", treeName="Events", isData=False, color=1, texName="DY FastSim", directory=directories['DYJetsToLL_M50_FS'])

TTJets      = Sample.fromDirectory(name="TTJets", treeName="Events", isData=False, color=1, texName="tt FullSim", directory=directories['TTJets'])
TTJets_FS   = Sample.fromDirectory(name="TTJets_FS", treeName="Events", isData=False, color=1, texName="tt FastSim", directory=directories['TTJets_FS'])

sample1 = TTJets
sample2 = TTJets_FS

# Sequence
read_variables =    ["weight/F","puWeight/F","PV_npvsGood/I",
                    #"jet[pt/F,eta/F,phi/F,btagCSV/F,DFb/F,DFbb/F,id/I,btagDeepCSV/F]", "njet/I","nJetSelected/I",
                    "MET_pt/F", "MET_phi/F", "MET_sumPt/F", "MET_significance/F", "MET_SignificanceRec/F","dl_mass/F", "MET_sumPt/F", "MET_sumEt/F",
                    #"Z_l1_index/I", "Z_l2_index/I", "nonZ_l1_index/I", "nonZ_l2_index/I",
                    #"Z_phi/F","Z_pt/F", "Z_mass/F", "Z_eta/F","Z_lldPhi/F", "Z_lldR/F"
]

sequence = []

# Plotting
lumi_scale = 35.9
noData = True

def drawObjects( plotData, dataMCScale, lumi_scale ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.15, 0.95, 'CMS Preliminary' if plotData else 'CMS Simulation'),
      (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) Scale %3.2f'% ( lumi_scale, dataMCScale ) ) if plotData else (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % lumi_scale)
    ]
    return [tex.DrawLatex(*l) for l in lines]

def drawPlots(plots, dataMCScale):
  for log in [False, True]:
    ext = "_small" if args.small else ""
    ext += "_log" if log else ""
    plot_directory_ = os.path.join(plot_directory, sample1.name, 'v1', ext)
    for plot in plots:
      if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot
      extensions_ = ["pdf", "png", "root"]

      plotting.draw(plot,
        plot_directory = plot_directory_,
        extensions = extensions_,
        ratio = {'yRange':(0.1,1.9)},
        logX = False, logY = log, sorting = False,
        scaling = {0:1},
        yRange = (0.03, "auto") if log else (0.001, "auto"),
        legend = [ (0.15,0.9-0.03*sum(map(len, plot.histos)),0.9,0.9), 2],
        drawObjects = drawObjects( not noData, dataMCScale , lumi_scale ),
        copyIndexPHP = True,
      )

# Samples
samples = [sample1, sample2]

sample1.color = ROOT.kRed+1
sample2.color = ROOT.kGreen+2
sample1.style = styles.lineStyle(color=ROOT.kRed+1)
sample2.style = styles.lineStyle(color=ROOT.kGreen+1)

for s in samples:
    #s.setSelectionString([getFilterCut(isData=False), tr.getSelection("MC")])
    if args.small: s.reduceFiles(to=3)
    s.read_variables = read_variables
    #s.weight         = lambda event, s: event.puWeight
    #s.style = styles.fillStyle(s.color)
    s.scale = lumi_scale

stack = Stack([sample1], [sample2])

weight_ = lambda event, sample: event.weight

selection = "dl_mass>0&&(Sum$(Muon_isGoodMuon&&Muon_pt>25&&abs(Muon_eta)<2.4)+Sum$(Electron_pt>25&&abs(Electron_eta)<2.4))==2"#&&Sum$(Electron_pt>15&&abs(Electron_eta)<2.4)==0"

Plot.setDefaults(stack = stack, weight = staticmethod(weight_), selectionString = selection, addOverFlowBin='upper')

plots = []

plots.append(Plot(
  name = 'dl_mass', texX = 'M(ll) (GeV)', texY = 'Number of Events',
  attribute = TreeVariable.fromString( "dl_mass/F" ),
  binning=[100,50.,150.],
))

plots.append(Plot(
  name = 'MET_pt', texX = 'MET (GeV)', texY = 'Number of Events',
  attribute = TreeVariable.fromString( "MET_pt/F" ),
  binning=[40,0.,400.],
))

plots.append(Plot(
  name = 'MET_sumEt', texX = 'MET Sum(E_{T}) (GeV)', texY = 'Number of Events',
  attribute = TreeVariable.fromString( "MET_sumEt/F" ),
  binning=[60,0.,3000.],
))

plots.append(Plot(
  name = 'MET_sumPt', texX = 'MET Sum(p_{T}) (GeV)', texY = 'Number of Events',
  attribute = TreeVariable.fromString( "MET_sumPt/F" ),
  binning=[50,0.,2000.],
))

plots.append(Plot(
  name = 'MET_significance', texX = 'MET sig (GeV)', texY = 'Number of Events',
  attribute = TreeVariable.fromString( "MET_significance/F" ),
  binning=[50,0.,50.],
))

plots.append(Plot(
  name = 'MET_SignificanceRec', texX = 'MET sig rec (GeV)', texY = 'Number of Events',
  attribute = TreeVariable.fromString( "MET_SignificanceRec/F" ),
  binning=[50,0.,50.],
))

plots.append(Plot(
  name = 'PV_npvsGood', texX = 'N_{PV good}', texY = 'Number of Events',
  attribute = TreeVariable.fromString( "PV_npvsGood/I" ),
  binning=[50,0.,50.],
))

plotting.fill(plots, read_variables = read_variables, sequence = sequence)

dataMCScale = 1
drawPlots(plots, dataMCScale)

