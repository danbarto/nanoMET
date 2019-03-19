#!/usr/bin/env python
''' Analysis script for standard plots
'''
#
# Standard imports and batch mode
#
import ROOT, os
ROOT.gROOT.SetBatch(True)
import itertools

from math                               import sqrt, cos, sin, pi
from RootTools.core.standard            import *
from nanoMET.tools.user                 import plot_directory
from nanoMET.samples.color              import color
from nanoMET.tools.puReweighting        import getReweightingFunction
from nanoMET.tools.cutInterpreter       import cutInterpreter
from Samples.Tools.metFilters           import getFilterCut
from nanoMET.tools.helpers              import deltaPhi
#from nanoMET.tools.cutInterpreter  import cutInterpreter

#
# Arguments
# 
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',            action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--noData',              action='store_true', default=False,           help='also plot data?')
argParser.add_argument('--small',               action='store_true',     help='Run only on a small subset of the data?', )
argParser.add_argument('--verySmall',               action='store_true',     help='Run only on a small subset of the data?', )
argParser.add_argument('--plot_directory',      action='store',      default='v8_noSigMax')
argParser.add_argument('--year',                action='store',      default=2016)
argParser.add_argument('--selection',           action='store',      default='looseLeptonVeto-onZ')
args = argParser.parse_args()



year = int(args.year)
#
# Logger
#
import StopsDilepton.tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

from nanoMET.core.JetResolution import *
from nanoMET.core.Event         import Event

if args.small:                        args.plot_directory += "_small"
if args.verySmall:                        args.plot_directory += "_verySmall"
if args.noData:                       args.plot_directory += "_noData"
#
# Make samples, will be searched for in the postProcessing directory
#
if year == 2016:
    postProcessing_directory = "2016_v6/dimuon/"
    from nanoMET.samples.nanoTuples_Summer16_postProcessed import *
    postProcessing_directory = "2016_v6/dimuon/"
    from nanoMET.samples.nanoTuples_Run2016_17Jul2018_postProcessed import *
    data_sample = DoubleMuon_Run2016
    mc          = [DY_LO_16, Top_16, VVTo2L2Nu_16, WJets_16]
    dy          = DY_LO_16
    top         = Top_16
    vv          = VVTo2L2Nu_16
    triggers    = ['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL', 'HLT_IsoMu24', 'HLT_IsoTkMu24']

    JERData     = JetResolution('Summer16_25nsV1_DATA')
    JERMC       = JetResolution('Summer16_25nsV1_MC')
    paramsData  = [1.843242937068234, 1.64107911184195, 1.567040591823117, 1.5077143780804294, 1.614014783345394, -0.0005986196920895609, 0.6071479349467596]
    paramsMC    = [1.617529475909303, 1.4505983036866312, 1.411498565372343, 1.4087559908291813, 1.3633674107893856, 0.0019861227075085516, 0.6539410816436597]

elif year == 2017:
    postProcessing_directory = "2017_v8/dimuon/"
    from nanoMET.samples.nanoTuples_Fall17_postProcessed import *
    postProcessing_directory = "2017_v8/dimuon/"
    from nanoMET.samples.nanoTuples_Run2017_31Mar2018_postProcessed import *
    data_sample = DoubleMuon_Run2017
    mc          = [DY_LO_17, Top_17, VVTo2L2Nu_17, WJets_17]
    dy          = DY_LO_17
    top         = Top_17
    vv          = VVTo2L2Nu_17
    triggers    = ['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ', 'HLT_IsoMu27']

    JERData     = JetResolution('Fall17_25nsV1_DATA')
    JERMC       = JetResolution('Fall17_25nsV1_MC')

## no EE jets in METSig (and tuning)
    paramsData  = [1.957222895585727, 1.9729235168459986, 1.90948681548405, 1.6360930500389212, 1.4258381193442973, 0.0004736982371569354, 0.6723294208092309]
    paramsMC    = [1.2241980571984066, 1.1869013615009936, 1.2102170428976644, 1.2081321118729322, 1.3823907717329618, -0.39694728791067135, 0.6841860339541063]

### no sig max in tuning
#    paramsData  = [1.966462143930697, 1.937183271447419, 1.8657691786964459, 1.6690843411711644, 1.4288816064105585, -0.00040472779927774936, 0.6722445090422673]
#    paramsMC    = [1.215607616120827, 1.2063625156082602, 1.2306271643996056, 1.2462203247952433, 1.393142812137987, 0.009617412335929804, 0.6826006402909576]

## EE Fix tuning
#    paramsData  = [1.8203205027930778, 1.831631046506774, 1.7677237626986682, 1.5918292402822833, 1.367835925196143, -0.00019188634750270158, 0.6536545242198241]
#    paramsMC    = [0.9407084182022997, 0.9523526365256799, 1.0575730226561166, 1.172357726178777, 1.322739950549703, -0.006175237894865501, 0.6658458528486568]

## old tuning
#    paramsData  = [1.743319492995906, 1.6882972548344242, 1.6551185757422577, 1.4185872885319166, 1.5923201986159454, -0.0002185734915505621, 0.6558819144933438]
#    paramsMC    = [0.7908154690397596, 0.8274420527567241, 0.8625204829478312, 0.9116933716967324, 1.1863207810108252, -0.0021905431583211926, 0.6620237657886061]
    
elif year == 2018:
    postProcessing_directory = "2018_v8/dimuon/"
    from nanoMET.samples.nanoTuples_Autumn18_postProcessed import *
    postProcessing_directory = "2018_v8/dimuon/"
    from nanoMET.samples.nanoTuples_Run2018_17Sep2018_postProcessed import *
    data_sample = DoubleMuon_Run2018
    mc          = [DY_LO_18, Top_18, diboson_18, rare_18]#, WJets_18]
    dy          = DY_LO_18
    top         = Top_18
    vv          = diboson_18
    triggers    = ['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8', 'HLT_IsoMu24']

    JERData     = JetResolution('Fall17_25nsV1_DATA')
    JERMC       = JetResolution('Fall17_25nsV1_MC')
    paramsData  = [1.8901832149541773, 2.026001195551111, 1.7805585857080317, 1.5987158841135176, 1.4509516794588302, 0.0003365079273751142, 0.6697617770737838]
    paramsMC    = [1.3889924894064565, 1.4100950862040742, 1.388614360360041, 1.2352876826748016, 1.0377595808114612, 0.004479319982990152, 0.6269386702181299]
#
# Text on the plots
#
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

#
# Read variables and sequences
#
read_variables = ["weight/F", "RawMET_pt/F", "RawMET_phi/F", "MET_pt/F", "MET_phi/F", "MET_sumPt/F", "fixedGridRhoFastjetAll/F", "Muon[pt/F,eta/F,phi/F,pfRelIso03_all/F,isGoodMuon/I]", "Jet[pt/F,eta/F,phi/F,cleanmask/O,cleanmaskMETSig/I,neEmEF/F,jetId/I,neHEF/F]", "nJet/I", "nPhoton/I","nMuon/I","nElectron/I"]
if year == 2018:
    read_variables += ["Jet[pt_nom/F]", "MET_pt_nom/F"]
#read_variables = ["weight/F", "MET_pt/F", "MET_phi/F", "MET_sumPt/F", "fixedGridRhoFastjetAll/F", "Muon[pt/F,eta/F,phi/F]", "Jet[pt/F,eta/F,phi/F,cleanmask/I,cleanmaskMETSig/I,neEmEF/F]", "nJet/I"]

#read_variables += ["Electron[pt/F,eta/F,phi/F,pfRelIso03_all/F]", "dl_mass/F", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8/O", "HLT_IsoMu24/O"]
#read_variables += ["Flag_goodVertices/O","Flag_globalSuperTightHalo2016Filter/O", "Flag_HBHENoiseFilter/O", "Flag_HBHENoiseIsoFilter/O", "Flag_EcalDeadCellTriggerPrimitiveFilter/O", "Flag_BadPFMuonFilter/O", "Flag_BadChargedCandidateFilter/O", "Flag_ecalBadCalibFilter/O", "Flag_eeBadScFilter/O"]


sequence = []

vetoEtaRegion = (2.65, 3.14) if year == False else (10,10)
jetThreshold = 25 if year == 2017 else 15

def calcMETSig( event, sample ):
    if sample.isData:
        JER = JERData
        params = paramsData
    else:
        JER = JERMC
        params = paramsMC
    if year == 2017: METCollection = "METFixEE2017_pt"
    else: METCollection = "MET_pt"
    if year == 2018: JetCollection = "Jet_pt_nom"
    else: JetCollection = "Jet_pt"
    ev = Event(event, JER, isData=sample.isData, METCollection=METCollection, JetCollection=JetCollection, vetoEtaRegion=vetoEtaRegion, jetThreshold=jetThreshold)
    ev.calcMETSig(params)
    event.MET_significance_rec = ev.MET_sig
    event.Jet_dpt = [ x for x in ev.Jet_dpt ]

def getMET_neEmEBalace( event, sample ):
    Jet_dPhiJetMET = []
    MET_min = []
    MET_minX = []
    MET_min_pseudoJet = []
    sumEt = 0
    badJet_Energy = 0

    MET_x = event.MET_pt * math.cos(event.MET_phi)
    MET_y = event.MET_pt * math.sin(event.MET_phi)
    
    pseudoJet_x = 0.
    pseudoJet_y = 0.
    pseudoJet_neEmEF = 0.
    pseudoJet_sumPt  = 0.

    EE_jets = []

    #print "Next event"
    for i in range(event.nJet):
        dPhiJetMET = deltaPhi(event.Jet_phi[i], event.MET_phi)
        #print event.Jet_phi[i], event.MET_phi, event.Jet_pt[i], event.Jet_neEmEF[i], dPhiJetMET
        if 2.5<abs(event.Jet_eta[i])<3.0 and (abs(dPhiJetMET - pi) < (pi/3.)):
            sumEt += -cos(dPhiJetMET)*event.Jet_pt[i]*event.Jet_neEmEF[i]
        if 2.5<abs(event.Jet_eta[i])<3.0 and event.Jet_pt[i] < 50:
            badJet_Energy += event.Jet_neEmEF[i]*event.Jet_pt[i]*math.cosh(event.Jet_eta[i])
        #print sumEt

        # do the cool stuff
        if 2.65<abs(event.Jet_eta[i])<3.14:# and event.Jet_neEmEF[i]*event.Jet_pt[i]*math.cosh(event.Jet_eta[i])>15:
            EE_jets.append({'pt':event.Jet_pt[i], 'eta':event.Jet_eta[i], 'phi':event.Jet_phi[i], 'neEmEF':event.Jet_neEmEF[i]})
            pseudoJet_sumPt += event.Jet_pt[i]
            pseudoJet_neEmEF += event.Jet_neEmEF[i]*event.Jet_pt[i]

            Jet_x = event.Jet_pt[i] * math.cos(event.Jet_phi[i])
            Jet_y = event.Jet_pt[i] * math.sin(event.Jet_phi[i])

            pseudoJet_x += Jet_x
            pseudoJet_y += Jet_y

            if event.Jet_neEmEF[i] > 0:
                alpha_j = 1 + (MET_x*Jet_x + MET_y*Jet_y)/(event.Jet_neEmEF[i]*(Jet_x**2 + Jet_y**2))
                alphaX_j = 1 + (MET_x*Jet_x + MET_y*Jet_y)/((Jet_x**2 + Jet_y**2))
            else:
                alpha_j = 1
                alphaX_j = 1
            #print alpha_j
            alpha_j = max(min(alpha_j,1),0)
            alphaX_j = max(min(alphaX_j,1),0)
            
            Jet_pt_min = math.sqrt((Jet_x*event.Jet_neEmEF[i]*alpha_j)**2 + (Jet_y*event.Jet_neEmEF[i]*alpha_j)**2)
            MET_min_j = math.sqrt((MET_x - (alpha_j-1) * event.Jet_neEmEF[i] * Jet_x)**2 + (MET_y - (alpha_j-1) * event.Jet_neEmEF[i] * Jet_y)**2)
            #print MET_min_j
            MET_minX_j = math.sqrt((MET_x - (alphaX_j-1) * Jet_x)**2 + (MET_y - (alphaX_j-1) * Jet_y)**2)
            MET_min.append((alpha_j, MET_min_j, event.Jet_pt[i], event.Jet_phi[i], event.Jet_eta[i], Jet_pt_min))
            MET_minX.append((alphaX_j, MET_minX_j, event.Jet_pt[i], event.Jet_phi[i], event.Jet_eta[i], Jet_pt_min))
    
    # get all combinations of pseudo-jets
    nEE_jets = len(EE_jets)
    MET_min_pseudoJet = copy.deepcopy(MET_min)
    nEE_jets = min(nEE_jets, 6)
    for i in range(2,nEE_jets+1):
        #print "Now pseudo-jets"
        combinations = itertools.combinations(EE_jets, i)
        for comb in combinations:
            pseudoJet_x = 0.
            pseudoJet_y = 0.
            pseudoJet_neEmEF = 0.
            pseudoJet_sumPt  = 0.

            for jet in comb:
                Jet_x = jet['pt'] * math.cos(jet['phi'])
                Jet_y = jet['pt'] * math.sin(jet['phi'])

                pseudoJet_x += Jet_x
                pseudoJet_y += Jet_y

                pseudoJet_sumPt     += jet['pt']
                pseudoJet_neEmEF    += jet['neEmEF'] * jet['pt']

            
            if pseudoJet_neEmEF>0:
                pseudoJet_neEmEF = pseudoJet_neEmEF/pseudoJet_sumPt
                alpha_j = 1 + (MET_x*pseudoJet_x + MET_y*pseudoJet_y)/(pseudoJet_neEmEF*(pseudoJet_x**2 + pseudoJet_y**2))
            else:
                alpha_j = 1
            alpha_j = max(min(alpha_j,1),0)
            MET_min_j = math.sqrt((MET_x - (alpha_j-1) * pseudoJet_neEmEF * pseudoJet_x)**2 + (MET_y - (alpha_j-1) * pseudoJet_neEmEF * pseudoJet_y)**2)
            MET_min_pseudoJet.append((alpha_j, MET_min_j, -99, -99, -99, -99))
            #print MET_min_j
    
    event.nJet_EE = len(MET_min)
    event.nJet_pseudoJets = len(MET_min_pseudoJet)
                   
    event.Jet_pt_EE     = -99
    event.Jet_eta_EE    = -99
    event.Jet_phi_EE    = -99
    event.Jet_pt_EE_corr = -99

    #print len(MET_min)

    if len(MET_min)>0:
        MET_min_t = min(MET_min, key = lambda t: t[1])
        MET_minX_t = min(MET_minX, key = lambda t: t[1])
        MET_min_pseudoJet_t = min(MET_min_pseudoJet, key = lambda t: t[1])

    if len(MET_min) == 0:
        # fill with defaults
        event.MET_min_MET40     = event.MET_pt
        event.alpha_MET40       = 1
        event.MET_min           = event.MET_pt
        event.MET_minX          = event.MET_pt
        event.MET_min_pseudoJet = event.MET_pt
        event.alpha             = 1
    elif event.MET_pt < 40:
        event.MET_min_MET40     = event.MET_pt
        event.alpha_MET40       = 1
        event.MET_min           = MET_min_t[1]
        event.MET_minX          = MET_minX_t[1]
        event.MET_min_pseudoJet = MET_min_pseudoJet_t[1]
        event.alpha             = MET_min_pseudoJet_t[0]
        event.Jet_pt_EE         = MET_min_t[2]
        event.Jet_phi_EE        = MET_min_t[3]
        event.Jet_eta_EE        = MET_min_t[4]
        event.Jet_pt_EE_corr    = MET_min_t[5]
    else:
        event.MET_min_MET40     = MET_min_t[1]
        event.alpha_MET40       = MET_min_t[0]
        event.MET_min           = MET_min_t[1]
        event.MET_minX          = MET_minX_t[1]
        event.MET_min_pseudoJet = MET_min_pseudoJet_t[1]
        event.alpha             = MET_min_pseudoJet_t[0]
        event.Jet_pt_EE         = MET_min_t[2]
        event.Jet_phi_EE        = MET_min_t[3]
        event.Jet_eta_EE        = MET_min_t[4]
        event.Jet_pt_EE_corr    = MET_min_t[5]
        #print MET_min, event.MET_pt

    event.MET_min_ratio = abs(event.MET_min/event.MET_pt - 1)
    event.MET_min_ratio_MET40 = abs(event.MET_min_MET40/event.MET_pt - 1)


    event.badJet_Energy = badJet_Energy
    if event.MET_pt > 0:
        event.MET_neEmEBalace = sumEt/event.MET_pt
    else:
        event.MET_neEmEBalace = 0
    #print event.MET_neEmEBalace

def getNJet( event, sample ):
    allJets = []
    for i in range(event.nJet):
        allJets.append({'pt':event.Jet_pt[i], 'eta':event.Jet_eta[i], 'phi':event.Jet_phi[i], 'neEmEF':event.Jet_neEmEF[i], 'jetId':event.Jet_jetId[i], 'cleanmask':event.Jet_cleanmask[i], 'neHEF':event.Jet_neHEF[i]})

    highPtJets = filter( lambda j: (j['pt']>30 and abs(j['eta'])<2.4 and j['jetId']>5 and j['cleanmask']>0), allJets )
    highPtJetsNoCleanmask = filter( lambda j: (j['pt']>30 and abs(j['eta'])<2.4 and j['jetId']>5), allJets )
    event.nJet_good = len(highPtJets)
    event.nJet_good_noCleanmask = len(highPtJetsNoCleanmask)
    event.nJet_noLeptons = len( filter( lambda j:j['cleanmask']>0, allJets) )
    lowPtJets = filter( lambda j: j['pt']<30., allJets)
    event.nJet_lowPt                = len( lowPtJets )
    event.nJet_lowPt_eta0to0p8      = len( filter( lambda j:abs(j['eta'])<0.8, lowPtJets) )
    event.nJet_lowPt_eta0p8to1p3    = len( filter( lambda j:0.8 <= abs(j['eta']) < 1.3, lowPtJets) )
    event.nJet_lowPt_eta1p3to1p9    = len( filter( lambda j:1.3 <= abs(j['eta']) < 1.9, lowPtJets) )
    event.nJet_lowPt_eta1p9to2p5    = len( filter( lambda j:1.9 <= abs(j['eta']) < 2.5, lowPtJets) )
    event.nJet_lowPt_eta2p5to3p1    = len( filter( lambda j:2.5 <= abs(j['eta']) < 3.1, lowPtJets) )
    event.nJet_lowPt_eta3p1toInf    = len( filter( lambda j:3.1 <= abs(j['eta']), lowPtJets) )

    sumPtneEmEF = 0
    sumPtneHEF  = 0
    for jet in filter( lambda j:2.5 <= abs(j['eta']) < 3.1, lowPtJets):
        sumPtneEmEF += jet['pt']*jet['neEmEF']
        sumPtneHEF  += jet['pt']*jet['neHEF']
        
    event.sumPtneEmEF_lowPt_eta2p5to3p1 = sumPtneEmEF
    event.sumPtneHEF_lowPt_eta2p5to3p1 = sumPtneHEF

sequence += [ calcMETSig, getMET_neEmEBalace, getNJet ]

def getLeptonSelection( mode ):
  if   mode=="mumu":
    if year == 2016:
        return "Sum$(Muon_pt>20&&Muon_isGoodMuon)==2&&Sum$(Muon_pt>25&&Muon_isGoodMuon)>0"
    else:
        # slower trigger turn-on in 2017&2018
        return "Sum$(Muon_pt>20&&Muon_isGoodMuon)==2&&Sum$(Muon_pt>35&&Muon_isGoodMuon)>0"

nTrueInt36fb_puRW        = getReweightingFunction(data="PU_2016_36000_XSecCentral", mc="Summer16")
nTrueInt36fb_puRWUp      = getReweightingFunction(data="PU_2016_36000_XSecUp",      mc="Summer16")

#
# Loop over channels
#
yields     = {}
allPlots   = {}
allModes   = ['mumu']

for index, mode in enumerate(allModes):
  yields[mode] = {}
  if   mode=="mumu":
    data_sample.texName = "data (2 #mu)"
    data_selectionString = '&&'.join([getFilterCut(isData=True, year=year).replace('&&weight>0',''), getLeptonSelection(mode), "( %s )"%" || ".join(triggers)])
    # maybe add singleMu backup
    print data_selectionString
    data_sample.setSelectionString([data_selectionString])
    data_sample.scale = 1
  if   mode=="mumu": data_sample.texName = "data (2 #mu)"

  data_sample.name           = "data"
  data_sample.read_variables = ["event/I","run/I","jsonPassed/I"]
  data_sample.style          = styles.errorStyle(ROOT.kBlack)
  lumi_scale                 = data_sample.lumi/1000

  if args.noData: lumi_scale = 35.9
  weight_ = lambda event, sample: event.weight

  #weightJetEn = lambda event, sample: event.Jet_pt*

  for sample in mc: sample.style = styles.fillStyle(sample.color)

  for sample in mc:
    sample.scale          = lumi_scale
    sample.read_variables = ['puWeight/F','puWeightUp/F', 'Pileup_nTrueInt/F']
    #sample.weight         = lambda event, sample: event.puWeight*nTrueInt36fb_puRW(event.Pileup_nTrueInt)
    sample.weight         = lambda event, sample: event.puWeight
    sample.setSelectionString([getFilterCut(isData=False, year=year), getLeptonSelection(mode), "( %s )"%" || ".join(triggers)])

  if not args.noData:
    stack = Stack(mc, data_sample)
  else:
    stack = Stack(mc)

  # seperate stacks for the 2D plots
  stackData = Stack(data_sample)
  stackDY   = Stack(dy)
  stackTT   = Stack(top)

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
  Plot.setDefaults(stack = stack, weight = staticmethod(weight_), selectionString = cutInterpreter.cutString(args.selection), addOverFlowBin='upper', histo_class=ROOT.TH1D)
  Plot2D.setDefaults(weight = weight_, selectionString = cutInterpreter.cutString(args.selection), histo_class=ROOT.TH2D)
  
  plots = []

  plots.append(Plot(
    name = 'yield', texX = 'yield', texY = 'Number of Events',
    attribute = lambda event, sample: 0.5 + index, # yields are somehow acting weird
    binning=[1, 0, 1],
  ))

  plots.append(Plot(
    name = 'nVtxs', texX = 'vertex multiplicity', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "PV_npvsGood/I" ),
    binning=[80,0,80],
  ))

  plots.append(Plot(
      texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
      attribute = TreeVariable.fromString( "MET_pt/F" ),
      binning=[400/20,0,400],
  ))

  plots.append(Plot(
      texX = 'E_{T}^{miss} (raw) (GeV)', texY = 'Number of Events / 20 GeV',
      attribute = TreeVariable.fromString( "RawMET_pt/F" ),
      binning=[400/20,0,400],
  ))

  plots.append(Plot(
      texX = 'E_{T}^{miss} neEmEBalance', texY = 'Number of Events',
      name = "MET_neEmEBalace",
      attribute = lambda event, sample: event.MET_neEmEBalace,
      binning=[80,0,2.0],
  ))

  plots.append(Plot(
      texX = 'E_{T}^{miss} min.', texY = 'Number of Events',
      name = "MET_min",
      attribute = lambda event, sample: event.MET_min,
      binning=[400/20,0,400],
  ))

  plots.append(Plot(
      texX = 'E_{T}^{miss} min. pseudo jet', texY = 'Number of Events',
      name = "MET_min_pseudoJet",
      attribute = lambda event, sample: event.MET_min_pseudoJet,
      binning=[400/20,0,400],
  ))

  plots.append(Plot(
      texX = 'E_{T}^{miss} min. all jet', texY = 'Number of Events',
      name = "MET_minX",
      attribute = lambda event, sample: event.MET_minX,
      binning=[400/20,0,400],
  ))

  plots.append(Plot(
      texX = 'E_{T}^{miss} min. (MET>40)', texY = 'Number of Events',
      name = "MET_min_MET40",
      attribute = lambda event, sample: event.MET_min_MET40,
      binning=[400/20,0,400],
  ))

  plots.append(Plot(
      texX = 'E_{T}^{miss} min. ratio', texY = 'Number of Events',
      name = "MET_min_ratio",
      attribute = lambda event, sample: event.MET_min_ratio,
      binning=[80,0,2.0],
  ))

  plots.append(Plot(
      texX = 'E_{T}^{miss} min. ratio (MET>40)', texY = 'Number of Events',
      name = "MET_min_ratio_MET40",
      attribute = lambda event, sample: event.MET_min_ratio_MET40,
      binning=[80,0,2.0],
  ))

  plots.append(Plot(
      texX = '#alpha', texY = 'Number of Events',
      name = "alpha",
      attribute = lambda event, sample: event.alpha,
      binning=[50,0,1.0],
  ))

  plots.append(Plot(
      texX = '#alpha (MET>40)', texY = 'Number of Events',
      name = "alpha_MET40",
      attribute = lambda event, sample: event.alpha_MET40,
      binning=[50,0,1.0],
  ))

  plots.append(Plot(
      texX = 'p_{T}(mis. Jet) (GeV)', texY = 'Number of Events',
      name = "Jet_pt_EE",
      attribute = lambda event, sample: event.Jet_pt_EE,
      binning=[50,0,200],
  ))

  plots.append(Plot(
      texX = 'p_{T}(mis. Jet, corr.) (GeV)', texY = 'Number of Events',
      name = "Jet_pt_EE_corr",
      attribute = lambda event, sample: event.Jet_pt_EE_corr,
      binning=[50,0,200],
  ))

  plots.append(Plot(
      texX = '#phi (mis. Jet)', texY = 'Number of Events',
      name = "Jet_phi_EE",
      attribute = lambda event, sample: event.Jet_phi_EE,
      binning=[64,-3.2,3.2],
  ))

  plots.append(Plot(
      texX = '#eta (mis. Jet)', texY = 'Number of Events',
      name = "Jet_eta_EE",
      attribute = lambda event, sample: event.Jet_eta_EE,
      binning=[64,-3.2,3.2],
  ))

  plots.append(Plot(
      texX = '#phi(E_{T}^{miss})', texY = 'Number of Events / 20 GeV',
      attribute = TreeVariable.fromString( "MET_phi/F" ),
      binning=[10,-pi,pi],
  ))

  plots.append(Plot(
      texX = '#phi(E_{T}^{miss}) (raw)', texY = 'Number of Events / 20 GeV',
      attribute = TreeVariable.fromString( "RawMET_phi/F" ),
      binning=[10,-pi,pi],
  ))

  if year == 2017:
    #METFixEE2017_pt
    
    plots.append(Plot(
        texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
        attribute = TreeVariable.fromString( "METFixEE2017_pt/F" ),
        binning=[400/20,0,400],
    ))

    plots.append(Plot(
        texX = '#phi(E_{T}^{miss})', texY = 'Number of Events / 20 GeV',
        attribute = TreeVariable.fromString( "METFixEE2017_phi/F" ),
        binning=[10,-pi,pi],
    ))
    
    plots.append(Plot(
        texX = 'badJet Energy', texY = 'Number of Events / 20 GeV',
        name = 'badJet_Energy',
        attribute = lambda event, sample: event.badJet_Energy,
        binning=[40,0,400],
    ))

  # removed from nanoAOD
  if year == 2018:
    plots.append(Plot(
      texX = 'E_{T}^{miss} Significance', texY = 'Number of Events',
      attribute = TreeVariable.fromString('MET_significance/F'),
      binning= [50,0,100],
    ))
    
    plots.append(Plot(
      texX = 'E_{T}^{miss} recorr', texY = 'Number of Events',
      attribute = TreeVariable.fromString('MET_pt_nom/F'),
      binning= [400/20,0,400],
    ))
    
  plots.append(Plot(
    texX = 'E_{T}^{miss} Significance rec.', texY = 'Number of Events',
    name = "MET_significance_rec",
    attribute = lambda event, sample: event.MET_significance_rec,
    binning= [50,0,100],
  ))

  plots.append(Plot(
    texX = '#Sigma p_{T} (GeV)', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "MET_sumPt/F" ),
    binning=[50, 0,2500],
  ))

  plots.append(Plot(
    texX = '#Sigma E_{T} (GeV)', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "MET_sumEt/F" ),
    binning=[40, 500,4500],
  ))

  plots.append(Plot(
    texX = 'm(ll) of leading dilepton (GeV)', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "dl_mass/F" ),
    binning=[40,81.2,101.2],
  ))

  plots.append(Plot(
    texX = 'm(ll) of leading dilepton (GeV)', texY = 'Number of Events',
    name = "dl_mass_wide",
    attribute = TreeVariable.fromString( "dl_mass/F" ),
    binning=[50,20,220],
  ))

  plots.append(Plot(
    texX = 'p_{T}(l_{1}) (GeV)', texY = 'Number of Events / 15 GeV',
    name = 'l1_pt',
    attribute = lambda event, sample: event.Muon_pt[0],
    binning=[20,0,300],
  ))

  plots.append(Plot(
    texX = '#phi(l_{1}) (GeV)', texY = 'Number of Events / 15 GeV',
    name = 'l1_phi',
    attribute = lambda event, sample: event.Muon_phi[0],
    binning=[16,-3.2,3.2],
  ))


  plots.append(Plot(
    texX = '#eta(l_{1}) (GeV)', texY = 'Number of Events / 15 GeV',
    name = 'l1_eta',
    attribute = lambda event, sample: event.Muon_eta[0],
    binning=[15,-3.0,3.0],
  ))

  plots.append(Plot(
    texX = 'p_{T}(l_{1}) (GeV)', texY = 'Number of Events / 15 GeV',
    name = 'l1_pt',
    attribute = lambda event, sample: event.Muon_pt[0],
    binning=[30,25,325],
  ))

  plots.append(Plot(
    texX = '#phi(l_{2}) (GeV)', texY = 'Number of Events / 15 GeV',
    name = 'l2_phi',
    attribute = lambda event, sample: event.Muon_phi[1],
    binning=[16,-3.2,3.2],
  ))

  plots.append(Plot(
    texX = '#eta(l_{2}) (GeV)', texY = 'Number of Events / 15 GeV',
    name = 'l2_eta',
    attribute = lambda event, sample: event.Muon_eta[1],
    binning=[15,-3.0,3.0],
  ))

  plots.append(Plot(
    texX = 'p_{T}(l_{2}) (GeV)', texY = 'Number of Events / 15 GeV',
    name = 'l2_pt',
    attribute = lambda event, sample: event.Muon_pt[1],
    binning=[28,20,300],
  ))

  plots.append(Plot(
    texX = 'N_{j} (EE)', texY = 'Number of Events',
    name = 'nJet_EE',
    attribute = lambda event, sample: event.nJet_EE,
    binning=[15,-0.5,14.5],
  ))

  plots.append(Plot(
    texX = 'N_{j} (EE, pseud-jets)', texY = 'Number of Events',
    name = 'nJet_pseudoJets',
    attribute = lambda event, sample: event.nJet_pseudoJets,
    binning=[20,0,500],
  ))

  plots.append(Plot(
    texX = 'N_{j} (all jets)', texY = 'Number of Events',
    name = 'nJet',
    attribute = lambda event, sample: event.nJet,
    binning=[40,0,40],
  ))

  plots.append(Plot(
    texX = 'N_{j} (p_{T}<30 GeV)', texY = 'Number of Events',
    name = 'nJet_lowPt',
    attribute = lambda event, sample: event.nJet_lowPt,
    binning=[40,0,40],
  ))

  plots.append(Plot(
    texX = 'N_{j} (p_{T}<30 GeV, |#eta|<0.8)', texY = 'Number of Events',
    name = 'nJet_lowPt_eta0to0p8',
    attribute = lambda event, sample: event.nJet_lowPt_eta0to0p8,
    binning=[20,0,20],
  ))

  plots.append(Plot(
    texX = 'N_{j} (p_{T}<30 GeV, 0.8<|#eta|<1.3)', texY = 'Number of Events',
    name = 'nJet_lowPt_eta0p8to1p3',
    attribute = lambda event, sample: event.nJet_lowPt_eta0p8to1p3,
    binning=[20,0,20],
  ))

  plots.append(Plot(
    texX = 'N_{j} (p_{T}<30 GeV, 1.3<|#eta|<1.9)', texY = 'Number of Events',
    name = 'nJet_lowPt_eta1p3to1p9',
    attribute = lambda event, sample: event.nJet_lowPt_eta1p3to1p9,
    binning=[20,0,20],
  ))

  plots.append(Plot(
    texX = 'N_{j} (p_{T}<30 GeV, 1.9<|#eta|<2.5)', texY = 'Number of Events',
    name = 'nJet_lowPt_eta1p9to2p5',
    attribute = lambda event, sample: event.nJet_lowPt_eta1p9to2p5,
    binning=[20,0,20],
  ))

  plots.append(Plot(
    texX = 'N_{j} (p_{T}<30 GeV, 2.5<|#eta|<3.1)', texY = 'Number of Events',
    name = 'nJet_lowPt_eta2p5to3p1',
    attribute = lambda event, sample: event.nJet_lowPt_eta2p5to3p1,
    binning=[20,0,20],
  ))

  plots.append(Plot(
    texX = 'N_{j} (p_{T}<30 GeV, 3.1<|#eta|)', texY = 'Number of Events',
    name = 'nJet_lowPt_eta3p1toInf',
    attribute = lambda event, sample: event.nJet_lowPt_eta3p1toInf,
    binning=[20,0,20],
  ))

  plots.append(Plot(
    texX = 'N_{j} (p_{T}>30 GeV, |#eta|<2.4)', texY = 'Number of Events',
    name = 'nJet_good',
    attribute = lambda event, sample: event.nJet_good,
    binning=[10,0,10],
  ))

  plots.append(Plot(
    texX = 'N_{j} (p_{T}>30 GeV, |#eta|<2.4, no cleanmask)', texY = 'Number of Events',
    name = 'nJet_good_noCleanmask',
    attribute = lambda event, sample: event.nJet_good_noCleanmask,
    binning=[10,0,10],
  ))

  plots.append(Plot(
    texX = 'N_{j} (all jets, no leptons)', texY = 'Number of Events',
    name = 'nJet_noLeptons',
    attribute = lambda event, sample: event.nJet_noLeptons,
    binning=[10,0,10],
  ))

  plots.append(Plot(
    texX = 'N_{e}', texY = 'Number of Events',
    name = 'nElectron',
    attribute = lambda event, sample: event.nElectron,
    binning=[10,0,10],
  ))

  plots.append(Plot(
    texX = 'N_{#mu}', texY = 'Number of Events',
    name = 'nMuon',
    attribute = lambda event, sample: event.nMuon,
    binning=[10,0,10],
  ))

  plots.append(Plot(
    texX = 'N_{#gamma}', texY = 'Number of Events',
    name = 'nPhoton',
    attribute = lambda event, sample: event.nPhoton,
    binning=[10,0,10],
  ))

  plots.append(Plot(
    texX = '#Sigma p_T*neEmEF (jets, p_{T}<30 GeV, 2.5<|#eta|<3.1)', texY = 'Number of Events',
    name = 'sumPtneEmEF_lowPt_eta2p5to3p1',
    attribute = lambda event, sample: event.sumPtneEmEF_lowPt_eta2p5to3p1,
    binning=[20,0,100],
  ))

  plots.append(Plot(
    texX = '#Sigma p_T*neHEF (jets, p_{T}<30 GeV, 2.5<|#eta|<3.1)', texY = 'Number of Events',
    name = 'sumPtneHEF_lowPt_eta2p5to3p1',
    attribute = lambda event, sample: event.sumPtneHEF_lowPt_eta2p5to3p1,
    binning=[20,0,100],
  ))

  ### 2D Data plots ###
  plots2D = []
  plots2D.append(Plot2D(
    stack = stackData,
    texX = 'E_{T}^{miss} (GeV)', texY = 'E_{T}^{miss, min} (GeV)',
    name = "Data_METvsMETmin",
    attribute = [ lambda event, sample: event.MET_pt, lambda event, sample: event.MET_min_pseudoJet ],
    binning = [40,0,200, 40, 0, 200]
  ))

  plots2D.append(Plot2D(
    stack = stackData,
    texX = '#eta(jet, EE)', texY = '#phi(jet, EE)',
    name = "Data_Jet_eta_EE_vs_Jet_phi_EE_plus",
    attribute = [ lambda event, sample: event.Jet_eta_EE, lambda event, sample: event.Jet_phi_EE ],
    binning = [20,2.6,3.2, 32, -3.2, 3.2]
  ))

  plots2D.append(Plot2D(
    stack = stackData,
    texX = '#eta(jet, EE)', texY = '#phi(jet, EE)',
    name = "Data_Jet_eta_EE_vs_Jet_phi_EE_minus",
    attribute = [ lambda event, sample: event.Jet_eta_EE, lambda event, sample: event.Jet_phi_EE ],
    binning = [20,-3.2,-2.6, 32, -3.2, 3.2]
  ))
  
  ### 2D DY plots ###
  plots2D.append(Plot2D(
    stack = stackDY,
    texX = 'E_{T}^{miss} (GeV)', texY = 'E_{T}^{miss, min} (GeV)',
    name = "DY_METvsMETmin",
    attribute = [ lambda event, sample: event.MET_pt, lambda event, sample: event.MET_min_pseudoJet ],
    binning = [40,0,200, 40, 0, 200]
  ))

  plots2D.append(Plot2D(
    stack = stackDY,
    texX = '#eta(jet, EE)', texY = '#phi(jet, EE)',
    name = "DY_Jet_eta_EE_vs_Jet_phi_EE_plus",
    attribute = [ lambda event, sample: event.Jet_eta_EE, lambda event, sample: event.Jet_phi_EE ],
    binning = [20,2.6,3.2, 32, -3.2, 3.2]
  ))

  plots2D.append(Plot2D(
    stack = stackDY,
    texX = '#eta(jet, EE)', texY = '#phi(jet, EE)',
    name = "DY_Jet_eta_EE_vs_Jet_phi_EE_minus",
    attribute = [ lambda event, sample: event.Jet_eta_EE, lambda event, sample: event.Jet_phi_EE ],
    binning = [20,-3.2,-2.6, 32, -3.2, 3.2]
  ))

  ### 2D DY plots ###
  plots2D.append(Plot2D(
    stack = stackTT,
    texX = 'E_{T}^{miss} (GeV)', texY = 'E_{T}^{miss, min} (GeV)',
    name = "Top_METvsMETmin",
    attribute = [ lambda event, sample: event.MET_pt, lambda event, sample: event.MET_min_pseudoJet ],
    binning = [40,0,200, 40, 0, 200]
  ))

  plots2D.append(Plot2D(
    stack = stackTT,
    texX = '#eta(jet, EE)', texY = '#phi(jet, EE)',
    name = "Top_Jet_eta_EE_vs_Jet_phi_EE_plus",
    attribute = [ lambda event, sample: event.Jet_eta_EE, lambda event, sample: event.Jet_phi_EE ],
    binning = [20,-3.2,2.6, 32, -3.2, 3.2]
  ))

  plots2D.append(Plot2D(
    stack = stackTT,
    texX = '#eta(jet, EE)', texY = '#phi(jet, EE)',
    name = "Top_Jet_eta_EE_vs_Jet_phi_EE_minus",
    attribute = [ lambda event, sample: event.Jet_eta_EE, lambda event, sample: event.Jet_phi_EE ],
    binning = [20,-3.2,-2.6, 32, -3.2, 3.2]
  ))

  plotting.fill( plots+plots2D, read_variables = read_variables, sequence = sequence)

  ## other 1D plots
  #print data_sample.selectionString
  #print dy.selectionString
  #print lumi_scale
  #for sample in [data_sample, dy, top, vv]:
  #  sample.chain.SetBranchStatus("*",1)
  #data_jetEta   = data_sample.get1DHistoFromDraw("Jet_eta", [50,-5.,5.], selectionString = cutInterpreter.cutString(args.selection), weightString="(1)")
  #dy_jetEta     = dy.get1DHistoFromDraw("Jet_eta", [50,-5.,5.], selectionString = cutInterpreter.cutString(args.selection), weightString="weight*puWeight*%s"%lumi_scale)
  #top_jetEta    = top.get1DHistoFromDraw("Jet_eta", [50,-5.,5.], selectionString = cutInterpreter.cutString(args.selection), weightString="weight*puWeight*%s"%lumi_scale)
  #vv_jetEta     = vv.get1DHistoFromDraw("Jet_eta", [50,-5.,5.], selectionString = cutInterpreter.cutString(args.selection), weightString="weight*puWeight*%s"%lumi_scale)
  
  plots1D = []
  #plots1D.append(Plot.fromHisto(
  #  name = "Jet_eta", texX="#eta(all jets)", texY="Events",
  #  histos = [[dy_jetEta,top_jetEta,vv_jetEta], [data_jetEta]]
  #))

  # Get normalization yields from yield histogram
  for plot in plots:
    if plot.name == "yield":
      for i, l in enumerate(plot.histos):
        for j, h in enumerate(l):
          yields[mode][plot.stack[i][j].name] = h.GetBinContent(h.FindBin(0.5+index))
          h.GetXaxis().SetBinLabel(1, "#mu#mu")
  if args.noData: yields[mode]["data"] = 0

  yields[mode]["MC"] = sum(yields[mode][s.name] for s in mc)
  dataMCScale        = yields[mode]["data"]/yields[mode]["MC"] if yields[mode]["MC"] != 0 else float('nan')

  drawPlots(plots, mode, dataMCScale)
  allPlots[mode] = plots

  for plot in plots2D:
    for log in [True, False]:
        plotting.draw2D(
            plot = plot,
            plot_directory = os.path.join(plot_directory, 'analysis_plots', str(year), args.plot_directory, mode + ("_log" if log else ""), args.selection),
            logX = False, logY = False, logZ = log,
            drawObjects = drawObjects( not args.noData, dataMCScale , lumi_scale ),
            copyIndexPHP = True
        )

  #for plot in plots1D:
  #  for log in [True, False]:
  #      plotting.draw(
  #          plot = plot,
  #          plot_directory = os.path.join(plot_directory, 'analysis_plots', str(year), args.plot_directory, mode + ("_log" if log else ""), args.selection),
  #          logX = False, logY = False,
  #          ratio={},
  #          drawObjects = drawObjects( not args.noData, dataMCScale , lumi_scale ),
  #          copyIndexPHP = True
  #      )

# Add the different channels into SF and all
if len(allModes) == 3: #just added for slimmedPlots since we only use doubleMu channel
    for mode in ["SF","all"]:
      yields[mode] = {}
      for y in yields[allModes[0]]:
        try:    yields[mode][y] = sum(yields[c][y] for c in (['ee','mumu'] if mode=="SF" else ['ee','mumu','mue']))
        except: yields[mode][y] = 0
      dataMCScale = yields[mode]["data"]/yields[mode]["MC"] if yields[mode]["MC"] != 0 else float('nan')
    
      for plot in allPlots['mumu']:
        for plot2 in (p for p in (allPlots['ee'] if mode=="SF" else allPlots["mue"]) if p.name == plot.name):  #For SF add EE, second round add EMu for all
          for i, j in enumerate(list(itertools.chain.from_iterable(plot.histos))):
    	    for k, l in enumerate(list(itertools.chain.from_iterable(plot2.histos))):
    	        if i==k:
    	            j.Add(l)
    
    drawPlots(allPlots['mumu'], mode, dataMCScale)

logger.info( "Done." )

