#!/usr/bin/env python
''' Analysis script for standard plots
'''
# Standard imports and batch mode
import ROOT, os
ROOT.gROOT.SetBatch(True)
import itertools
import numpy

from math                               import sqrt, cos, sin, pi

from RootTools.core.standard            import *

from nanoMET.tools.user                 import plot_directory
from nanoMET.tools.puReweighting        import getReweightingFunction
from nanoMET.tools.cutInterpreter       import cutInterpreter
from nanoMET.tools.metFilters           import getFilterCut
from nanoMET.tools.helpers              import deltaPhi

from nanoMET.samples.color              import color

from nanoMET.core.JetResolution         import *
from nanoMET.core.Event                 import Event


# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',            action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--noData',              action='store_true', default=False,           help='also plot data?')
argParser.add_argument('--small',               action='store_true', help='Run only on a small subset of the data?', )
argParser.add_argument('--verySmall',           action='store_true', help='Run only on a small subset of the data?', )
argParser.add_argument('--plot_directory',      action='store',      default='v1')
argParser.add_argument('--year',                action='store',      default=2016, type=int)
argParser.add_argument('--selection',           action='store',      default='looseLeptonVeto-onZ-nCleanJet1p')
args = argParser.parse_args()

# Logger
import StopsDilepton.tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.small:     args.plot_directory += "_small"
if args.verySmall: args.plot_directory += "_verySmall"
if args.noData:    args.plot_directory += "_noData"

# Make samples, will be searched for in the postProcessing directory
if year == 2016:
    postProcessing_directory = "2016_v22/dimuon/"
    from nanoMET.samples.nanoTuples_Summer16_postProcessed import *
    postProcessing_directory = "2016_v22/dimuon/"
    from nanoMET.samples.nanoTuples_Run2016_17Jul2018_postProcessed import *
    data_sample = DoubleMuon_Run2016
    mc          = [DY_LO_16, Top_16, diboson_16, rare_16]
    dy          = DY_LO_16
    top         = Top_16
    triggers    = ['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL', 'HLT_IsoMu24', 'HLT_IsoTkMu24']

    JERData     = JetResolution('Summer16_25nsV1_DATA')
    JERMC       = JetResolution('Summer16_25nsV1_MC')
    paramsData  = [1.843242937068234, 1.64107911184195, 1.567040591823117, 1.5077143780804294, 1.614014783345394, -0.0005986196920895609, 0.6071479349467596]
    paramsMC    = [1.617529475909303, 1.4505983036866312, 1.411498565372343, 1.4087559908291813, 1.3633674107893856, 0.0019861227075085516, 0.6539410816436597]

elif year == 2017:
    postProcessing_directory = "2017_v22/dimuon/"
    from nanoMET.samples.nanoTuples_Fall17_postProcessed import *
    postProcessing_directory = "2017_v22/dimuon/"
    from nanoMET.samples.nanoTuples_Run2017_31Mar2018_postProcessed import *
    data_sample = DoubleMuon_Run2017
    mc          = [DY_LO_17, Top_17, diboson_17, rare_17]
    dy          = DY_LO_17
    top         = Top_17
    triggers    = ['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ', 'HLT_IsoMu27']

    JERData     = JetResolution('Fall17_V3_DATA')
    JERMC       = JetResolution('Fall17_V3_MC')
    paramsData  = [1.743319492995906, 1.6882972548344242, 1.6551185757422577, 1.4185872885319166, 1.5923201986159454, -0.0002185734915505621, 0.6558819144933438]
    paramsMC    = [0.7908154690397596, 0.8274420527567241, 0.8625204829478312, 0.9116933716967324, 1.1863207810108252, -0.0021905431583211926, 0.6620237657886061]
    
elif year == 2018:
    postProcessing_directory = "2018_v22/dimuon/"
    from nanoMET.samples.nanoTuples_Autumn18_postProcessed import *
    postProcessing_directory = "2018_v22/dimuon/"
    from nanoMET.samples.nanoTuples_Run2018_17Sep2018_postProcessed import *
    data_sample = DoubleMuon_Run2018
    mc          = [DY_LO_18, Top_18, diboson_18, rare_18]
    dy          = DY_LO_18
    top         = Top_18
    triggers    = ['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8', 'HLT_IsoMu24']

    JERData     = JetResolution('Autumn18_V1_DATA')
    JERMC       = JetResolution('Autumn18_V1_MC')
    paramsData  = [1.8901832149541773, 2.026001195551111, 1.7805585857080317, 1.5987158841135176, 1.4509516794588302, 0.0003365079273751142, 0.6697617770737838]
    paramsMC    = [1.3889924894064565, 1.4100950862040742, 1.388614360360041, 1.2352876826748016, 1.0377595808114612, 0.004479319982990152, 0.6269386702181299]

# Text on the plots
def drawObjects( plotData, dataMCScale, lumi_scale ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.15, 0.95, 'CMS Preliminary' if plotData else 'CMS Simulation'), 
      (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) Scale %3.2f'% ( lumi_scale, dataMCScale ) ) if plotData else (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % lumi_scale)
    ]
    if "mt2ll100" in args.selection and args.noData: lines += [(0.55, 0.5, 'M_{T2}(ll) > 100 GeV')] # Manually put the mt2ll > 100 GeV label
    return [tex.DrawLatex(*l) for l in lines] 

def drawPlots(plots, mode, dataMCScale):
  for log in [False, True]:
    plot_directory_ = os.path.join(plot_directory, 'analysis_plots', str(year), args.plot_directory, mode + ("_log" if log else ""), args.selection)
    for plot in plots:
      if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot
      if not args.noData: 
        if mode == "all": plot.histos[1][0].legendText = "Data"
        if mode == "SF":  plot.histos[1][0].legendText = "Data (SF)"

      plotting.draw(plot,
	    plot_directory = plot_directory_,
	    ratio = {'yRange':(0.1,1.9)} if not args.noData else None,
	    logX = False, logY = log, sorting = True,
	    yRange = (0.03, "auto") if log else (0.001, "auto"),
	    scaling = {},
	    legend = (0.50,0.88-0.04*sum(map(len, plot.histos)),0.9,0.88) if not args.noData else (0.50,0.9-0.047*sum(map(len, plot.histos)),0.85,0.9),
	    drawObjects = drawObjects( not args.noData, dataMCScale , lumi_scale ),
        copyIndexPHP = True,
      )

# Read variables and sequences
read_variables = ["weight/F", "MET_pt/F", "MET_phi/F", "MET_sumPt/F",
                  "fixedGridRhoFastjetAll/F", "Muon[pt/F,eta/F,phi/F]",
                  "Jet[pt/F,eta/F,phi/F,cleanmask/I,cleanmaskMETSig/I,neEmEF/F,jetId/I]",
                  "nJet/I", "GenMET_pt/F", "PV_npvsGood/I"
                 ]

sequence = []

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
    event.Jet_dpt = [ x for x in ev.Jet_dpt ]

def calcHT( event, sample ):
    HT = 0
    for i in range(event.nJet):
        if event.Jet_pt[i]>30 and abs(event.Jet_eta[i])<2.4 and event.Jet_jetId[i]>0:
            HT += event.Jet_pt[i]
    event.HT = HT
    if HT>0:
        event.MET_significance_eventBased = event.MET_pt/math.sqrt(HT)
    else:
        event.MET_significance_eventBased = 100 # what to assign? significance should be high

sequence += [ calcMETSig, calcHT ]

def getLeptonSelection( mode ):
    if mode=="mumu":
        if year == 2016:
            return cutInterpreter.cutString("diMuon16")
        else:
            # slower trigger turn-on in 2017&2018
            return cutInterpreter.cutString("diMuon1718")

nTrueInt36fb_puRW        = getReweightingFunction(data="PU_2016_36000_XSecCentral", mc="Summer16")
nTrueInt36fb_puRWUp      = getReweightingFunction(data="PU_2016_36000_XSecUp",      mc="Summer16")

# Loop over channels
yields     = {}
allPlots   = {}
allModes   = ['mumu']

for index, mode in enumerate(allModes):
    
    lumi_scale = 35.9
    weight_ = lambda event, sample: event.weight

    for sample in mc: sample.style = styles.fillStyle(sample.color)

    for sample in mc:
      sample.scale          = lumi_scale
      sample.read_variables = ['puWeight/F','Pileup_nTrueInt/F']
      sample.weight         = lambda event, sample: event.puWeight
      sample.setSelectionString([getFilterCut(isData=False, year=year), getLeptonSelection(mode), "(%s)"%"||".join(triggers)])

    stack = Stack(mc)

    if args.small:
      for sample in stack.samples:
          sample.normalization=1.
          sample.reduceFiles( factor=40 )
          sample.scale /= sample.normalization
          #sample.reduceFiles( to = 3 )
    if args.verySmall:
      for sample in stack.samples:
          sample.normalization=1.
          sample.reduceFiles( to = 1 )
          sample.scale /= sample.normalization

    # Use some defaults
    dummyWeight_ = lambda event, sample: float(0) 
    Plot.setDefaults(stack=stack, weight = staticmethod(dummyWeight_), selectionString = cutInterpreter.cutString(args.selection), histo_class=ROOT.TH1D)
    
    # 1D dummy plot
    dummyPlot = Plot(
        name = 'ROC', texX = 'signal efficiency', texY = 'background rejection',
        attribute = lambda event, sample: 0.5 + index, # yields are somehow acting weird
        binning=[1, 0, 1],
    )

    # Use some defaults
    Plot2D.setDefaults(weight = weight_, selectionString = cutInterpreter.cutString(args.selection), histo_class=ROOT.TH2D)

    ### 2D DY plots ###
    plots2D = []

    plots2D.append(Plot2D(
      stack = stack,
      texX = 'MET Sig', texY = 'gen MET',
      name = "METSig_vs_genMET",
      attribute = [ lambda event, sample: event.MET_significance_rec, lambda event, sample: event.GenMET_pt>1 ],
      binning = [200,0,100, 2, 0, 2]
    ))

    plots2D.append(Plot2D(
      stack = stack,
      texX = 'MET pt', texY = 'gen MET',
      name = "MET_pt_vs_genMET",
      attribute = [ lambda event, sample: event.MET_pt, lambda event, sample: event.GenMET_pt>1 ],
      binning = [200,0,400, 2, 0, 2]
    ))

    plots2D.append(Plot2D(
      stack = stack,
      texX = 'MET/sqrt(HT)', texY = 'gen MET',
      name = "METSig_eventBased_vs_genMET",
      attribute = [ lambda event, sample: event.MET_significance_eventBased, lambda event, sample: event.GenMET_pt>1 ],
      binning = [100,0,30, 2, 0, 2]
    ))

    plots2D.append(Plot2D(
      stack = stack,
      texX = 'MET Sig', texY = 'N_{PV}',
      name = "METSig_vs_PV_npvsGood",
      attribute = [ lambda event, sample: event.MET_significance_rec, TreeVariable.fromString( "PV_npvsGood/I" ) ],
      binning = [50,0,50, 100, 0, 100]
    ))

    plots2D.append(Plot2D(
      stack = stack,
      texX = 'MET pt', texY = 'N_{PV}',
      name = "MET_pt_vs_PV_npvsGood",
      attribute = [ lambda event, sample: event.MET_pt, TreeVariable.fromString( "PV_npvsGood/I" ) ],
      binning = [50,0,250, 100, 0, 100]
    ))

    plots2D.append(Plot2D(
      stack = stack,
      texX = 'MET/sqrt(HT)', texY = 'N_{PV}',
      name = "METSig_eventBased_vs_PV_npvsGood",
      attribute = [ lambda event, sample: event.MET_significance_eventBased, TreeVariable.fromString( "PV_npvsGood/I" ) ],
      binning = [50,0,25, 100, 0, 100]
    ))

    plotting.fill( plots2D + [dummyPlot], read_variables = read_variables, sequence = sequence)

    for plot in plots2D:
      for log in [True, False]:
          plotting.draw2D(
              plot = plot,
              plot_directory = os.path.join(plot_directory, 'analysis_plots', str(year), args.plot_directory, mode + ("_log" if log else ""), args.selection),
              logX = False, logY = False, logZ = log,
              drawObjects = drawObjects( not args.noData, 1.00 , lumi_scale ),
              copyIndexPHP = True
          )

    colors = [ROOT.kGreen+2, ROOT.kRed+1, ROOT.kBlue+2]
    metSigRoc = []
    for j,p in enumerate(plots2D[:3]):
        signal  = [0]*p.binning[0]
        bkg     = [0]*p.binning[0]
        for h in p.histos[0]:

            for i in range(p.binning[0]):
                signal[i]   += h.GetBinContent(i+1, 2)
                bkg[i]      += h.GetBinContent(i+1, 1)

        # get cumulative
        cumulativeSignal    = numpy.cumsum(signal)
        cumulativeBkg       = numpy.cumsum(bkg)
        totalSignal         = cumulativeSignal[-1]
        totalBkg            = cumulativeBkg[-1]

        revCumSignal        = [totalSignal]*p.binning[0]
        revCumBkg           = [totalBkg]*p.binning[0]


        for i in range(p.binning[0]):
            revCumSignal[i] += -cumulativeSignal[i]
            revCumBkg[i]    += -cumulativeBkg[i]

        tmp = ROOT.TGraph()
        for i in range(p.binning[0]):
            tmp.SetPoint(i, revCumSignal[i]/totalSignal, 1-revCumBkg[i]/totalBkg)
        tmp.SetPoint(p.binning[0],0,1)
        #tmp.style          = styles.lineStyle(colors[j])
        tmp.SetLineColor(colors[j])
        tmp.SetLineWidth(2)
        tmp.SetMarkerSize(0)
        metSigRoc.append(copy.deepcopy(tmp))
        
    plotting.draw(
        plot = dummyPlot,
        plot_directory = os.path.join(plot_directory, 'analysis_plots', str(year), args.plot_directory, mode, args.selection),
        logX = False, logY = False,
        legend = None,
        yRange = (0.0, 1.0),
        drawObjects = drawObjects( not args.noData, 1.00 , lumi_scale ) + metSigRoc,
    )


logger.info( "Done." )

