''' Analysis script for 1D 2l plots (RootTools)
'''

#Standard imports
import ROOT
from math import sqrt, cos, sin, pi, acos
import itertools,os
import copy

import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',            action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--year',                action='store',      default=2016)
argParser.add_argument('--tuneEra',             action='store',      default=False)
argParser.add_argument('--small', action='store_true', help='Small?')
args = argParser.parse_args()


#RootTools
from RootTools.core.standard import *

plot_directory              = "/afs/hephy.at/user/d/dspitzbart/www/nanoMET/validation/"

if args.tuneEra:                      plot_directory += "tune%s/"%args.tuneEra

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
postProcessing_directory    = "2017_v8/dimuon/"
dirs = {}
dirs['DYJetsToLL']          = ["DYJetsToLL_M50_LO"]
dirs['TTLep_pow']          = ["TTLep_pow"]
directories = { key : [ os.path.join( data_directory, postProcessing_directory, dir) for dir in dirs[key]] for key in dirs.keys()}

from nanoMET.samples.nanoTuples_Run2017_31Mar2018_postProcessed import *
DYJetsToLL_17  = Sample.fromDirectory(name="DYJets_17", treeName="Events", isData=False, color=color.DY, texName="DY", directory=directories['DYJetsToLL'])
Top_17              = Sample.fromDirectory(name="Top_17",              treeName="Events", isData=False, color=color.TTJets,        texName="t(#bar{t})",               directory=directories['TTLep_pow'])


year = int(args.year)

from nanoMET.core.JetResolution import *
from nanoMET.core.Event         import Event

if year == 2016:
    JERData     = JetResolution('Summer16_25nsV1_DATA')
    JERMC       = JetResolution('Summer16_25nsV1_MC')
    sample1 = DYJetsToLL_16
elif year == 2017:
    sample1 = DYJetsToLL_17
    #sample1 = DoubleMuon_Run2017
    #sample1 = Top_17
    JERData     = JetResolution('Fall17_25nsV1_DATA')
    JERMC       = JetResolution('Fall17_25nsV1_MC')


#selection = "Sum$(Muon_pt>25&&Muon_isGoodMuon)==2 && Sum$(Electron_pt>10&&abs(Electron_eta)<2.5&&Electron_cutBased>0&&abs(Electron_pfRelIso03_all)<0.4)==0 && Sum$(Jet_pt>30&&Jet_jetId&&abs(Jet_eta)<2.4)>=0 && abs(dl_mass-91.2)<10"
selection = "Sum$(Muon_pt>25&&Muon_isGoodMuon)==2 && Sum$(Electron_pt>10&&abs(Electron_eta)<2.5&&Electron_cutBased>0&&abs(Electron_pfRelIso03_all)<0.4)==0 && Sum$(Jet_pt>30&&Jet_jetId&&abs(Jet_eta)<2.4)>=0"
#if not sample1.isData: selection += "&& GenMET_pt<10"
#selection = "nElectron==0 && nMuon==0 && GenMET_pt<10"


## Sequence
read_variables = ["weight/F", "RawMET_pt/F", "RawMET_phi/F", "MET_pt/F", "MET_phi/F", "MET_sumPt/F", "fixedGridRhoFastjetAll/F", "Muon[pt/F,eta/F,phi/F,pfRelIso03_all/F,isGoodMuon/I]", "Jet[pt/F,eta/F,phi/F,cleanmask/O,cleanmaskMETSig/I,neEmEF/F,jetId/I,neHEF/F]", "nJet/I", "nPhoton/I","nMuon/I","nElectron/I"]
if year == 2017:
    read_variables += ["METFixEE2017_pt/F", "METFixEE2017_phi/F"]
if year == 2018:
    read_variables += ["Jet[pt_nom/F]", "MET_pt_nom/F"]

sequence = []


tuneParams = {
    2016: {
            'data': [1.843242937068234, 1.64107911184195, 1.567040591823117, 1.5077143780804294, 1.614014783345394, -0.0005986196920895609, 0.6071479349467596],
            'mc':   [1.617529475909303, 1.4505983036866312, 1.411498565372343, 1.4087559908291813, 1.3633674107893856, 0.0019861227075085516, 0.6539410816436597]
            },
    20167: {
            'data': [1.843242937068234, 1.64107911184195, 1.567040591823117, 1.5077143780804294, 1.614014783345394, -0.0005986196920895609, 0.7502485988002288],
            'mc':   [1.617529475909303, 1.4505983036866312, 1.411498565372343, 1.4087559908291813, 1.3633674107893856, 0.0019861227075085516, 0.7397929625473758]
            },
    2017: {
            'data': [1.5622583144490318, 1.540194388842639, 1.566197393467264, 1.5132500586113067, 1.9398948489956538, -0.00028476818675941427, 0.7502485988002288],
            'mc':   [1.1106319092645605, 1.1016751869920842, 1.0725643000703053, 1.0913641155398053, 1.8499497840145123, -0.0015646588275905911, 0.7397929625473758]
            },
    2018: {
            'data': [1.4416589250258958, 1.592549070456071, 1.400707548171599, 1.4213958262324593, 2.0868635081187348, -0.0007745499117034968, 0.7261267509272097],
            'mc':   [1.0117455874431338, 1.2986232760320007, 1.1414800394963855, 0.9209396460085367, 1.312067503147174, 0.0012299929784571964, 0.681951027334837]
            },
    0: {
            'data': [0.7,0.7,0.7,0.7,0.7,0.,0.75],
            'mc':   [0.7,0.7,0.7,0.7,0.7,0.,0.75]
        },
    1: {
            'data': [1.,1.,1.,1.,1.,0.,0.75],
            'mc':   [1.,1.,1.,1.,1.,0.,0.75]
        },
    2: {
            'data': [1.5, 1.5, 1.5, 1.5, 1.5,0.,0.75],
            'mc':   [1.5, 1.5, 1.5, 1.5, 1.5,0.,0.75]
        },
    3: {
            'data': [2.,2.,2.,2.,2.,0.,0.75],
            'mc':   [2.,2.,2.,2.,2.,0.,0.75]
        },
    4: {
            'data': [4.,4.,4.,4.,4.,0.,0.75],
            'mc':   [4.,4.,4.,4.,4.,0.,0.75]
        }
    }

tuneEra = year if not args.tuneEra else int(args.tuneEra)

paramsData  = tuneParams[tuneEra]['data']
paramsMC    = tuneParams[tuneEra]['mc']

vetoEtaRegion = (2.65, 3.14) if year == False else (10,10)
jetThreshold = 25 if (year == 2017 or year == 2018) else 15


def calcMETSig( event, sample ):
    if sample.isData:
        JER = JERData
        params = paramsData
    else:
        JER = JERMC
        params = paramsMC
    if year == 2017:
        METPtVar = "METFixEE2017_pt"
        METPhiVar = "METFixEE2017_phi"
    elif year == 2018:
        METPtVar = "MET_pt_nom"
        METPhiVar = "MET_phi"
    else:
        METPtVar = "MET_pt"
        METPhiVar = "MET_phi"
    if year == 2018: JetCollection = "Jet_pt_nom"
    else: JetCollection = "Jet_pt"
    ev = Event(event, JER, isData=sample.isData, METPtVar=METPtVar, METPhiVar=METPhiVar, JetCollection=JetCollection, vetoEtaRegion=vetoEtaRegion, jetThreshold=jetThreshold)
    ev.calcMETSig(params)
    event.MET_significance_rec = ev.MET_sig
    event.Jet_dpt = [ x for x in ev.Jet_dpt ]

#    if sample.isData:
#        JER = JERData
#        params = paramsData
#    else:
#        JER = JERMC
#        params = paramsMC
#    ev = Event(event, JER, isData=sample.isData)
#    ev.calcMETSig(params)
#    event.MET_significance_rec = ev.MET_sig

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
    plot_directory_ = os.path.join(plot_directory, args.year, sample1.name + ext)
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



