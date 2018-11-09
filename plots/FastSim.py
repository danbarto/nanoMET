''' Analysis script for 1D 2l plots (RootTools)
'''

#Standard imports
import ROOT
from math import sqrt, cos, sin, pi, acos
import itertools,os
import copy

import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
args = argParser.parse_args()


#RootTools
from RootTools.core.standard import *

data_directory              = "/afs/hephy.at/data/dspitzbart03/nanoAOD/"
postProcessing_directory    = "dimuon/"
plot_directory              = "/afs/hephy.at/user/d/dspitzbart/www/nanoMET/FastSim/"

import nanoMET.tools.logger as logger

import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)


dirs = {}
dirs['DYJetsToLL_M50']          = ["DYJetsToLL_M50_MLM_S16_80X_priv"]
dirs['DYJetsToLL_M50_FS'] = ['DYJetsToLL_M50_MLM_FS_S16_80X_priv']
directories = { key : [ os.path.join( data_directory, postProcessing_directory, dir) for dir in dirs[key]] for key in dirs.keys()}

sample1 = Sample.fromDirectory(name="DY", treeName="Events", isData=False, color=1, texName="DY FullSim", directory=directories['DYJetsToLL_M50'])
sample2 = Sample.fromDirectory(name="DY_FS", treeName="Events", isData=False, color=1, texName="DY FastSim", directory=directories['DYJetsToLL_M50_FS'])

## Sequence
read_variables =    ["weight/F","puWeight/F",
                    #"jet[pt/F,eta/F,phi/F,btagCSV/F,DFb/F,DFbb/F,id/I,btagDeepCSV/F]", "njet/I","nJetSelected/I",
                    "MET_pt/F", "MET_phi/F", "MET_sumPt/F", "MET_significance/F", "MET_SignificanceRec/F","dl_mass/F",
                    #"Z_l1_index/I", "Z_l2_index/I", "nonZ_l1_index/I", "nonZ_l2_index/I",
                    #"Z_phi/F","Z_pt/F", "Z_mass/F", "Z_eta/F","Z_lldPhi/F", "Z_lldR/F"
]

sequence = []

## Plotting
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
    ext = "_small" if small else ""
    ext += "_log" if log else ""
    plot_directory_ = os.path.join(plot_directory, 'v1')
    for plot in plots:
      if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot
      extensions_ = ["pdf", "png", "root"]

      plotting.draw(plot,
        plot_directory = plot_directory_,
        extensions = extensions_,
        ratio = {'yRange':(0.1,1.9)} if not noData else None,
        logX = False, logY = log, sorting = False,
        yRange = (0.03, "auto") if log else (0.001, "auto"),
        legend = [ (0.15,0.9-0.03*sum(map(len, plot.histos)),0.9,0.9), 2],
        drawObjects = drawObjects( not noData, dataMCScale , lumi_scale ),
        copyIndexPHP = True,
      )

# Samples
samples = [sample1, sample2]

sample1.color = ROOT.kRed+1
sample2.color = ROOT.kGreen+2

for s in samples:
    #s.setSelectionString([getFilterCut(isData=False), tr.getSelection("MC")])
    s.read_variables = read_variables
    s.weight         = lambda event, s: event.puWeight
    s.style = styles.fillStyle(s.color)
    s.scale = lumi_scale

stack = Stack([sample1], [sample2])

weight_ = lambda event, sample: event.weight

selection = "dl_mass>0"

Plot.setDefaults(stack = stack, weight = staticmethod(weight_), selectionString = selection, addOverFlowBin='upper')

plots = []

plots.append(Plot(
  name = 'dl_mass', texX = 'M(ll) (GeV)', texY = 'Number of Events',
  attribute = TreeVariable.fromString( "dl_mass/F" ),
  binning=[200,0.,200.],
))


plotting.fill(plots, read_variables = read_variables, sequence = sequence)

dataMCScale = 1
drawPlots(plots, dataMCScale)

#for sample in [sample1, sample2]:
#    sample.read_variables = 


