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
argParser.add_argument('--small', action='store_true', help='Small?')
args = argParser.parse_args()


#RootTools
from RootTools.core.standard import *

plot_directory              = "/afs/hephy.at/user/d/dspitzbart/www/nanoMET/validation/"

import nanoMET.tools.logger as logger

import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

from nanoMET.samples.color import *

data_directory              = "/afs/hephy.at/data/dspitzbart03/nanoSamples/"
postProcessing_directory    = "2016_v6/Inclusive/"
dirs = {}
dirs['WJetsToLNu']      = ["WJetsToLNu_ext"]
dirs['TTZToQQ']         = ['TTZToQQ'] 
directories = { key : [ os.path.join( data_directory, postProcessing_directory, dir) for dir in dirs[key]] for key in dirs.keys()}

WJetsToLNu_16  = Sample.fromDirectory(name="WJets", treeName="Events", isData=False, color=color.WJets, texName="W+Jets", directory=directories['WJetsToLNu'])
TTZToQQ_16     = Sample.fromDirectory(name="TTZToQQ", treeName="Events", isData=False, color=color.TTZ, texName="t#bar{t}Z, had", directory=directories['TTZToQQ'])

data_directory              = "/afs/hephy.at/data/dspitzbart03/nanoSamples/"
postProcessing_directory    = "2016_v6/dimuon/"
dirs = {}
dirs['DYJetsToLL']          = ["DYJetsToLL_M50_LO_ext1"]
directories = { key : [ os.path.join( data_directory, postProcessing_directory, dir) for dir in dirs[key]] for key in dirs.keys()}

DYJetsToLL_16  = Sample.fromDirectory(name="DYJets_16", treeName="Events", isData=False, color=color.DY, texName="DY", directory=directories['DYJetsToLL'])


data_directory              = "/afs/hephy.at/data/dspitzbart03/nanoSamples/"
postProcessing_directory    = "2017_v6/dimuon/"
dirs = {}
dirs['DYJetsToLL']          = ["DYJetsToLL_M50_LO"]
directories = { key : [ os.path.join( data_directory, postProcessing_directory, dir) for dir in dirs[key]] for key in dirs.keys()}

DYJetsToLL_17  = Sample.fromDirectory(name="DYJets_17", treeName="Events", isData=False, color=color.DY, texName="DY", directory=directories['DYJetsToLL'])

sample1 = DYJetsToLL_17

selection = "Sum$(Muon_pt>25&&Muon_isGoodMuon)==2 && Sum$(Electron_pt>10&&abs(Electron_eta)<2.5&&Electron_cutBased>0&&abs(Electron_pfRelIso03_all)<0.4)==0 && Sum$(Jet_pt>30&&Jet_jetId&&abs(Jet_eta)<2.4)>=0 && abs(dl_mass-91.2)<10 && GenMET_pt<10"
#selection = "nElectron==0 && nMuon==0 && GenMET_pt<10"


## Sequence
read_variables = ["weight/F", "MET_pt/F", "MET_phi/F", "MET_sumPt/F", "fixedGridRhoFastjetAll/F", "Muon[pt/F,eta/F,phi/F]", "Jet[pt/F,eta/F,phi/F,cleanmask/I,cleanmaskMETSig/I]", "nJet/I"]

sequence = []

from nanoMET.core.JetResolution import *
from nanoMET.core.Event         import Event
#JERData = JetResolution('Spring16_25nsV6_DATA')
#JERMC   = JetResolution('Spring16_25nsV6_MC')
#paramsData  = [1.38, 1.27, 1.22, 1.16, 1.10, 0.0, 0.58]
#paramsMC    = [1.39, 1.26, 1.21, 1.23, 1.28, -0.26, 0.62]

# 2016 tuning v3
JERData     = JetResolution('Summer16_25nsV1_DATA')
JERMC       = JetResolution('Summer16_25nsV1_MC')
paramsData  = [1.843242937068234, 1.64107911184195, 1.567040591823117, 1.5077143780804294, 1.614014783345394, -0.0005986196920895609, 0.6071479349467596]
paramsMC    = [1.617529475909303, 1.4505983036866312, 1.411498565372343, 1.4087559908291813, 1.3633674107893856, 0.0019861227075085516, 0.6539410816436597]

# 2017 tuning v3
JERData     = JetResolution('Fall17_25nsV1_DATA')
JERMC       = JetResolution('Fall17_25nsV1_MC')
paramsData  = []
paramsMC    = [0.7908154690397596, 0.8274420527567241, 0.8625204829478312, 0.9116933716967324, 1.1863207810108252, -0.0021905431583211926, 0.6620237657886061]


def calcMETSig( event, sample ):
    if sample.isData:
        JER = JERData
        params = paramsData
    else:
        JER = JERMC
        params = paramsMC
    ev = Event(event, JER, isData=sample.isData)
    ev.calcMETSig(params)
    event.MET_significance_rec = ev.MET_sig

sequence += [ calcMETSig ]


## Plotting
lumi_scale = 35.9
noData = True

small = args.small


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

def drawPlots(plots, dataMCScale, objects=None):
  dO = drawObjects( not noData, dataMCScale , lumi_scale )
  if objects:
    dO += copy.deepcopy(objects)

  for log in [True, False]:
    ext = "_small" if small else ""
    ext += "_log" if log else ""
    plot_directory_ = os.path.join(plot_directory, sample1.name + ext)
    for plot in plots:
      if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot
      extensions_ = ["pdf", "png", "root"]

      plotting.draw(plot,
        plot_directory = plot_directory_,
        extensions = extensions_,
        #ratio = {'yRange':(0.1,1.9)},
        logX = False, logY = log, sorting = False,
        #scaling = {0:1},
        yRange = (0.03, "auto") if log else (0.001, "auto"),
        legend = (0.50,0.88-0.05*sum(map(len, plot.histos)),0.9,0.88),
        drawObjects = dO,
        copyIndexPHP = True,
      )

# Samples
samples = [sample1]

sample1.style = styles.fillStyle(color=sample1.color)



for s in samples:
    #s.setSelectionString([getFilterCut(isData=False), tr.getSelection("MC")])
    if small: s.reduceFiles(to=1)
    s.read_variables = read_variables
    #s.weight         = lambda event, s: event.puWeight
    #s.style = styles.fillStyle(s.color)
    s.scale = lumi_scale

stack = Stack([sample1])

weight_ = lambda event, sample: event.weight

Plot.setDefaults(stack = stack, weight = staticmethod(weight_), selectionString = selection, addOverFlowBin='upper')
dataMCScale = 1

#plots = []
#
##plots.append(Plot(
##  name = 'dl_mass', texX = 'M(ll) (GeV)', texY = 'Number of Events',
##  attribute = TreeVariable.fromString( "dl_mass/F" ),
##  binning=[100,50.,150.],
##))
#
#plots.append(Plot(
#  name = 'MET_pt', texX = 'MET (GeV)', texY = 'Number of Events',
#  attribute = TreeVariable.fromString( "MET_pt/F" ),
#  binning=[40,0.,400.],
#))
#
#plots.append(Plot(
#  name = 'MET_sumEt', texX = 'MET Sum(E_{T}) (GeV)', texY = 'Number of Events',
#  attribute = TreeVariable.fromString( "MET_sumEt/F" ),
#  binning=[60,0.,3000.],
#))
#
#plots.append(Plot(
#  name = 'MET_sumPt', texX = 'MET Sum(p_{T}) (GeV)', texY = 'Number of Events',
#  attribute = TreeVariable.fromString( "MET_sumPt/F" ),
#  binning=[50,0.,2000.],
#))
#
##plots.append(Plot(
##  name = 'MET_significance', texX = 'MET sig (GeV)', texY = 'Number of Events',
##  attribute = TreeVariable.fromString( "MET_significance/F" ),
##  binning=[50,0.,50.],
##))
#
#plots.append(Plot(
#  name = 'PV_npvsGood', texX = 'N_{PV good}', texY = 'Number of Events',
#  attribute = TreeVariable.fromString( "PV_npvsGood/I" ),
#  binning=[50,0.,50.],
#))
#
#plotting.fill_with_draw(plots)
#
#drawPlots(plots, dataMCScale)

metSigPlots = []
maxSig = 100.
nBins = 50

metSigPlots.append(Plot(
  name = 'MET_SignificanceRec', texX = 'MET sig rec (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.MET_significance_rec,
  binning=[nBins,0.,maxSig],
))

plotting.fill(metSigPlots, read_variables = read_variables, sequence = sequence)

const = metSigPlots[0].histos[0][0].GetBinContent(1)
const = const/ROOT.TMath.Exp(-(float(maxSig)/nBins)/4.)
chi2_2 = ROOT.TF1("chi22","TMath::Exp(-x/2)*"+str(const),0,maxSig)
chi2_2.SetLineColor(ROOT.kRed)
chi2_2.SetLineWidth(2)

drawPlots(metSigPlots, dataMCScale, [chi2_2])



