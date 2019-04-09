#!/usr/bin/env python
''' Analysis script for standard plots
'''
#
# Standard imports and batch mode
#
import ROOT
ROOT.gROOT.SetBatch(True)

from math                           import sqrt, cos, sin, pi
from RootTools.core.standard        import *
from RootTools.plot.helpers         import copyIndexPHP
from nanoMET.tools.user             import plot_directory
from Samples.Tools.metFilters       import getFilterCut
from nanoMET.tools.cutInterpreter   import cutInterpreter
from nanoMET.tools.lock             import waitForLock, removeLock

import pickle, os, time
import errno
#
# Arguments
# 
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',          action='store',      default='INFO',     nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--signal',            action='store',      default='None',        nargs='?', choices=['None', 'ewkDM'], help="Add signal to plot")
argParser.add_argument('--noData',            action='store_true', default=False,       help='also plot data?')
argParser.add_argument('--plot_directory',    action='store',      default='v1')
argParser.add_argument('--selection',         action='store',            default='looseLeptonVeto-onZ')
argParser.add_argument('--selectSys',         action='store',      default='all')
#argParser.add_argument('--noMultiThreading',  action='store_true', default='False', help="noMultiThreading?") # Need no multithreading when doing batch-to-natch
argParser.add_argument('--showOnly',          action='store',      default=None)
argParser.add_argument('--small',             action='store_true',     help='Run only on a small subset of the data?', )
argParser.add_argument('--runLocal',             action='store_true',     help='Run local or submit?', )
argParser.add_argument('--isChild',           action='store_true', default=False)
argParser.add_argument('--normalizeBinWidth', action='store_true', default=False,       help='normalize wider bins?')
argParser.add_argument('--dryRun',            action='store_true', default=False,       help='do not launch subjobs')
argParser.add_argument("--year",              action='store', type=int,      default=2016, choices = [ 2016, 2017, 2018 ], help='Which year?')
args = argParser.parse_args()


#
# Logger
#
import nanoMET.tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

jetSelection    = "nJetSelected"
bJetSelectionM  = "nBTag"

logger.info("Starting")

#
# Systematics to run over
#
jet_systematics    = ['jesTotalUp','jesTotalDown']# 'JERDown','JECVUp','JECVDown']
met_systematics    = ['unclustEnUp', 'unclustEnDown']
jme_systematics    = jet_systematics + met_systematics
weight_systematics = ['puWeightUp', 'puWeightDown']

if args.selectSys != "all" and args.selectSys != "combine": all_systematics = [args.selectSys if args.selectSys != 'None' else None]
else:                                                       all_systematics = [None] + weight_systematics + jme_systematics
#else:                                                       all_systematics = [None] + jet_systematics


sys_pairs = [\
    ('JES',         'jesTotalUp', 'jesTotalDown'),
    ('unclustEn',   'unclustEnUp', 'unclustEnDown'),
    ('puWeight',    'puWeightUp', 'puWeightDown')
]

#
# If this is the mother process, launch the childs and exit (I know, this could potententially be dangereous if the --isChild and --selection commands are not given...)
#

def wrapper(com):
  import os
  os.system(com)

#if not args.isChild and args.selection is None and (args.selectSys == "all" or args.selectSys == "combine"):
if not args.isChild and (args.selectSys == "all" or args.selectSys == "combine"):
  jobs = []
  for sys in (all_systematics if args.selectSys == "all" else ["combine"]):
    command = "python systematicsPlots.py --selection=" + args.selection + (" --noData" if args.noData else "")\
               + (" --isChild")\
               + (" --small" if args.small else "")\
               + (" --plot_directory=" + args.plot_directory)\
               + (" --logLevel=" + args.logLevel)\
               + (" --selectSys=" + str(sys))\
               + (" --signal=" + args.signal)\
               + (" --year=" + str(args.year))\
               + (" --normalizeBinWidth" if args.normalizeBinWidth else "")
    if args.selectSys == 'combine':
        jobs.append(command)
    elif args.selectSys == 'all':
        if args.runLocal:
            jobs.append(command)
        else:
            jobs.append( "submitBatch.py --title='sys' '%s'"%command )

#  if args.noMultiThreading: 
  logger.info("Running/submitting all systematics.")
  results = map(wrapper, jobs)
  logger.info("Done with running/submitting systematics.")
  exit(0)

if args.noData:                   args.plot_directory += "_noData"
if args.signal == "DM":           args.plot_directory += "_DM"
if args.signal == "T2tt":         args.plot_directory += "_T2tt"
if args.small:                    args.plot_directory += "_small"

year = int(args.year)
try: os.makedirs(os.path.join(plot_directory, 'systematicsPlots', str(args.year), args.plot_directory))
except: pass

logger.info("Logger still working? 1")
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
    triggers    = ['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL', 'HLT_IsoMu24', 'HLT_IsoTkMu24']

elif year == 2017:
    postProcessing_directory = "2017_v11/dimuon/"
    from nanoMET.samples.nanoTuples_Fall17_postProcessed import *
    postProcessing_directory = "2017_v11/dimuon/"
    from nanoMET.samples.nanoTuples_Run2017_31Mar2018_postProcessed import *
    data_sample = DoubleMuon_Run2017
    mc          = [DY_LO_17, Top_17, VVTo2L2Nu_17, WJets_17]
    dy          = DY_LO_17
    top         = Top_17
    triggers    = ['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ', 'HLT_IsoMu27']

elif year == 2018:
    postProcessing_directory = "2018_v11/dimuon/"
    from nanoMET.samples.nanoTuples_Autumn18_postProcessed import *
    postProcessing_directory = "2018_v11/dimuon/"
    from nanoMET.samples.nanoTuples_Run2018_17Sep2018_postProcessed import *
    data_sample = DoubleMuon_Run2018
    mc          = [DY_LO_18, Top_18]#, VVTo2L2Nu_18, WJets_18]
    dy          = DY_LO_18
    top         = Top_18
    triggers    = ['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8', 'HLT_IsoMu24']


logger.info("Logger still working? 2")
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
      (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) Scale %3.2f'% ( lumi_scale, dataMCScale ) ) if False else (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % lumi_scale)
    ]
    return [tex.DrawLatex(*l) for l in lines] 


def addSys( selectionString , sys = None ):
    if   sys in jet_systematics: return selectionString.replace('Jet_pt', 'Jet_pt_' + sys).replace('MET_pt', 'MET_pt_' + sys)
    elif sys in met_systematics: return selectionString.replace('MET_pt', 'MET_pt_' + sys)
    else:                        return selectionString


def weightMC( sys = None ):
    # weights used right now: PU, BTag, trigger, tracking
    if sys is None:
        return (lambda event, sample:event.weight*event.puWeight, "weight * puWeight")
        #weight_ = lambda event, sample: event.puWeight
        #return (staticmethod(weight_), "puWeight")#event.weight*event.puWeight,"weight * puWeight")
        #return weight_
    elif 'puWeight' in sys:
        return (lambda event, sample:event.weight*getattr(event, sys), "weight * "+sys)
    elif sys in jme_systematics :
        return weightMC( sys = None )
    else:
        raise ValueError( "Systematic %s not known"%sys )
    
#
# Read variables and sequences
#
read_variables = ["weight/F", "Jet[pt/F,eta/F,phi/F]",
                  "MET_pt/F", "MET_phi/F",
                  ]

sequence = []

logger.info("Logger still working? 3")


def getLeptonSelection( mode ):
  if   mode=="mumu":
    if year == 2016:
        return "Sum$(Muon_pt>20&&Muon_isGoodMuon)==2&&Sum$(Muon_pt>25&&Muon_isGoodMuon)>0"
    else:
        # slower trigger turn-on in 2017&2018
        return "Sum$(Muon_pt>20&&Muon_isGoodMuon)==2&&Sum$(Muon_pt>35&&Muon_isGoodMuon)>0"


#
# Loop over channels
#
allPlots   = {}
allModes   = ['mumu']
logger.info("Logger still working?")
for index, mode in enumerate(allModes):
    logger.info("Modes")
    logger.info('Working on mode ' + str(mode))


    data_selectionString = '&&'.join([getFilterCut(isData=True, year=year).replace('&&weight>0',''), getLeptonSelection(mode), "( %s )"%" || ".join(triggers)])
    data_sample.setSelectionString([data_selectionString])
    data_sample.scale = 1

    lumi_scale = data_sample.lumi/1000
    data_sample.texName = "data"
    data_sample.read_variables = ['weight/F']
    data_sample.style   = styles.errorStyle( ROOT.kBlack )

    data_weight = lambda event, sample: event.weight
    data_weight_string = "weight"

    print lumi_scale
    logger.info('Lumi scale is ' + str(lumi_scale))


    if args.small:
        for sample in mc + [data_sample]:# + ([data_sample] if type(data_sample)!=type([]) else data_sample):
            sample.reduceFiles( to = 1 )#factor = 40 )

    #def weightMC2( ):
    #    def func( event, sample ):
    #        return event.weight
    #    return func
    #weight_ = weightMC2()

    weight_ = lambda event,sample: event.weight

    for sample in mc:
        sample.scale           = lumi_scale
        sample.style           = styles.fillStyle(sample.color, lineColor = sample.color)
        sample.read_variables  = ["puWeight/F"]
        sample.read_variables += ["%s/F"%s for s in weight_systematics]
        #sample.read_variables += ["Jet_pt_%s/F"%s   for s in jet_systematics]# + ["MET_pt_%s/F"%s   for s in jet_systematics]# + ["MET_pt_min_%s/F"%s   for s in jet_systematics]
        if year == 2017:
            sample.read_variables += ["METFixEE2017_pt_nom/F"]
            sample.read_variables += ["METFixEE2017_pt_%s/F"%s   for s in jet_systematics]
            sample.read_variables += ["METFixEE2017_pt_%s/F"%s   for s in met_systematics]
        else:
            sample.read_variables += ["MET_pt_nom/F"]
            sample.read_variables += ["MET_pt_%s/F"%s   for s in jet_systematics]
            sample.read_variables += ["MET_pt_%s/F"%s   for s in met_systematics]
        #sample.read_variables += ["MET_pt_%s/F"%s   for s in jet_systematics] + ["MET_pt_min_%s/F"%s    for s in jet_systematics] + ["MET_significance_%s/F"%s    for s in jet_systematics]
        sample.read_variables += ["MET_significance_%s/F"%s    for s in jet_systematics]
        sample.read_variables += ["MET_significance_%s/F"%s    for s in met_systematics]

        #sample.setWeightString("weight*puWeight")
        #sample.weight = lambda event, sample: event.puWeight
        sample.setSelectionString([getFilterCut(isData=False, year=args.year), getLeptonSelection(mode), "( %s )"%" || ".join(triggers)])

    # Use some defaults
    #Plot.setDefaults( selectionString = cutInterpreter.cutString(args.selection) )
  
    stack_mc        = Stack( mc )

    stack_data = Stack( data_sample )
    sys_stacks = {sys:copy.deepcopy(stack_mc) for sys in [None] + weight_systematics + jme_systematics }
    plots = []
  

    if year == 2016 or year == 2018:
        ## MET_pt
        if not args.noData and (args.selectSys == 'None' or args.selectSys == 'combine'):
            met_data  = Plot( 
                name            = "MET_pt_data",
                texX            = 'E_{T}^{miss} (GeV)', 
                texY            = 'Number of Events' if args.normalizeBinWidth else "Number of Event / 20 GeV",
                stack           = stack_data, 
                attribute       = TreeVariable.fromString( "MET_pt/F" ),
                binning         = [20,0,400],
                selectionString = cutInterpreter.cutString(args.selection),
                weight          = data_weight,
                )
            plots.append( met_data )
        else:
            met_data = None

        met_mc = {}
        for sys in all_systematics:
            met_mc[sys] = Plot(
                name            = "MET_pt" if sys is None else "MET_pt_mc_%s" % sys,
                texX            = 'E_{T}^{miss} (GeV)',
                texY            = 'Number of Events' if args.normalizeBinWidth else "Number of Event / 20 GeV",
                stack           = sys_stacks[sys],
                attribute       = TreeVariable.fromString('MET_pt/F') if sys not in jme_systematics else TreeVariable.fromString( "MET_pt_%s/F" % sys ),
                binning         = [20,0,400],
                selectionString = addSys(cutInterpreter.cutString(args.selection), sys),
                weight          = weightMC(sys)[0]#weight_,#lambda event, sample: int(1),
                )
        plots.extend( met_mc.values() )


    if year == 2017:
        ## METFixEE2017_pt
        if not args.noData and (args.selectSys == 'None' or args.selectSys == 'combine'):
            met_data  = Plot( 
                name            = "MET_pt_data",
                texX            = 'E_{T}^{miss} (GeV)', 
                texY            = 'Number of Events' if args.normalizeBinWidth else "Number of Event / 20 GeV",
                stack           = stack_data, 
                attribute       = TreeVariable.fromString( "METFixEE2017_pt_nom/F" ),
                binning         = [20,0,400],
                selectionString = cutInterpreter.cutString(args.selection),
                weight          = data_weight,
                )
            plots.append( met_data )
        else:
            met_data = None

        met_mc = {}
        for sys in all_systematics:
            met_mc[sys] = Plot(
                name            = "MET_pt" if sys is None else "MET_pt_mc_%s" % sys,
                texX            = 'E_{T}^{miss} (GeV)',
                texY            = 'Number of Events' if args.normalizeBinWidth else "Number of Event / 20 GeV",
                stack           = sys_stacks[sys],
                attribute       = TreeVariable.fromString('METFixEE2017_pt_nom/F') if sys not in jme_systematics else TreeVariable.fromString( "METFixEE2017_pt_%s/F" % sys ),
                binning         = [20,0,400],
                selectionString = addSys(cutInterpreter.cutString(args.selection), sys),
                weight          = weightMC(sys)[0]#weight_,#lambda event, sample: int(1),
                )
        plots.extend( met_mc.values() )


    ## MET_significance
    if not args.noData and (args.selectSys == 'None' or args.selectSys == 'combine'):
        met_sig_data  = Plot( 
            name            = "MET_significance_data",
            texX            = 'S(E_{T}^{miss})', 
            texY            = 'Number of Events',
            stack           = stack_data, 
            attribute       = TreeVariable.fromString( "MET_significance/F" ),
            binning         = [20,0,100],
            selectionString = cutInterpreter.cutString(args.selection),
            weight          = data_weight,
            )
        plots.append( met_sig_data )
    else:
        met_sig_data = None

    met_sig_mc = {}
    for sys in all_systematics:
        met_sig_mc[sys] = Plot(
            name            = "MET_significance" if sys is None else "MET_significance_mc_%s" % sys,
            texX            = 'S(E_{T}^{miss})',
            texY            = 'Number of Events',
            stack           = sys_stacks[sys],
            attribute       = TreeVariable.fromString('MET_significance/F') if sys not in jme_systematics else TreeVariable.fromString( "MET_significance_%s/F" % sys ),
            binning         = [20,0,100],
            selectionString = addSys(cutInterpreter.cutString(args.selection), sys),
            weight          = weightMC(sys)[0]#weight_,#lambda event, sample: int(1),
            )
    plots.extend( met_sig_mc.values() )


    ## MET_significance fine binning
    if not args.noData and (args.selectSys == 'None' or args.selectSys == 'combine'):
        met_sig_data_fine  = Plot( 
            name            = "MET_significance_data",
            texX            = 'S(E_{T}^{miss})', 
            texY            = 'Number of Events',
            stack           = stack_data, 
            attribute       = TreeVariable.fromString( "MET_significance/F" ),
            binning         = [50,0,100],
            selectionString = cutInterpreter.cutString(args.selection),
            weight          = data_weight,
            )
        plots.append( met_sig_data_fine )
    else:
        met_sig_data_fine = None

    met_sig_mc_fine = {}
    for sys in all_systematics:
        met_sig_mc_fine[sys] = Plot(
            name            = "MET_significance" if sys is None else "MET_significance_mc_%s" % sys,
            texX            = 'S(E_{T}^{miss})',
            texY            = 'Number of Events',
            stack           = sys_stacks[sys],
            attribute       = TreeVariable.fromString('MET_significance/F') if sys not in jme_systematics else TreeVariable.fromString( "MET_significance_%s/F" % sys ),
            binning         = [50,0,100],
            selectionString = addSys(cutInterpreter.cutString(args.selection), sys),
            weight          = weightMC(sys)[0]#weight_,#lambda event, sample: int(1),
            )
    plots.extend( met_sig_mc_fine.values() )

    ### MET_pt_min
    #if not args.noData and (args.selectSys == 'None' or args.selectSys == 'combine'):
    #    met_min_data  = Plot(
    #        name            = "MET_pt_min_data",
    #        texX            = 'E_{T}^{miss,min} (GeV)',
    #        texY            = 'Number of Events' if args.normalizeBinWidth else "Number of Event / 20 GeV",
    #        stack           = stack_data,
    #        attribute       = TreeVariable.fromString( "MET_pt_min/F" ),
    #        binning         = [20,0,400],
    #        selectionString = cutInterpreter.cutString(args.selection),
    #        weight          = data_weight,
    #        )
    #    plots.append( met_min_data )
    #else:
    #    met_min_data = None

    #met_min_mc = {}
    #for sys in all_systematics:
    #    met_min_mc[sys] = Plot(
    #        name            = "MET_pt_min" if sys is None else "MET_pt_min_mc_%s" % sys,
    #        texX            = 'E_{T}^{miss,min} (GeV)',
    #        texY            = 'Number of Events' if args.normalizeBinWidth else "Number of Event / 20 GeV",
    #        stack           = sys_stacks[sys],
    #        attribute       = TreeVariable.fromString('MET_pt_min/F') if sys not in jme_systematics else TreeVariable.fromString( "MET_pt_min_%s/F" % sys ),
    #        binning         = [20,0,400],
    #        selectionString = addSys(cutInterpreter.cutString(args.selection), sys),
    #        weight          = weightMC(sys)[0]#weight_,#lambda event, sample: int(1),
    #        )
    #plots.extend( met_min_mc.values() )


    plotConfigs = [\
            [ met_mc, met_data, 1],
            [ met_sig_mc, met_sig_data, 1],
            [ met_sig_mc_fine, met_sig_data_fine, 1],
            #[ met_min_mc, met_min_data, 1],

    ]

    print args.year
    result_dir  = os.path.join(plot_directory, "systematicsPlots", str(args.year), args.plot_directory, mode, args.selection)
    result_file = os.path.join(result_dir, 'results.pkl')
    try: os.makedirs(result_dir)
    except: pass

    ## get the norm etc - not needed for ttZ!
    if args.selectSys != "combine": 
      
        plotting.fill(plots, read_variables = read_variables, sequence = sequence)

        waitForLock( result_file ) 
        if os.path.exists(result_file):
            allPlots = pickle.load(file( result_file ))
            allPlots.update({p.name : p.histos for p in plots})
        else:                           
            allPlots = {p.name : p.histos for p in plots}
        print result_file
        pickle.dump( allPlots, file( result_file, 'w' ) )
        removeLock( result_file ) 
        logger.info( "Done for sys " + args.selectSys )

    else:
        print result_file
        allPlots = pickle.load(file( result_file ))

        from RootTools.plot.Plot import addOverFlowBin1D
        for p in plots:
            p.histos = allPlots[p.name]
            for s in p.histos:
                for h in s:
                    addOverFlowBin1D(h, "upper")
                    if h.Integral()==0: logger.warning( "Found empty histogram %s in results file %s", h.GetName(), result_file )

        for plot_mc, plot_data, bin_width in plotConfigs:
            if args.normalizeBinWidth and bin_width>0:
                for p in plot_mc.values() + [plot_data]:
                    for histo in sum(p.histos, []):
                        histo.Scale(1,"width")

            
            # For now, keep Toms way of adding shape uncertainties (the following part is only used for that)
            if None in plot_mc.keys():
                shapeHists = {comp:plot_mc[None].histos[0][mc.index(comp)] for comp in mc}
            else:
                print "Couldn't find central histogram! Taking %s insted."%plot_mc.keys()[0]
                shapeHists = {comp:plot_mc[plot_mc.keys()[0]].histos[0][comp] for comp in mc}

            print shapeHists

            #Calculating systematics
            h_summed = {k: plot_mc[k].histos_added[0][0].Clone() for k in plot_mc.keys()}

            h_rel_err = h_summed[None].Clone()
            h_rel_err.Reset()

            #MC statistical error
            for ib in range( 1 + h_rel_err.GetNbinsX() ):
                h_rel_err.SetBinContent(ib, h_summed[None].GetBinError(ib)**2 )

            h_sys = {}
            goOn = False
            for k, s1, s2 in ([s for s in sys_pairs if s[0] == args.showOnly] if args.showOnly else sys_pairs):
              goOn = True
              h_sys[k] = h_summed[s1].Clone()
              h_sys[k].Scale(-1)
              h_sys[k].Add(h_summed[s2])
            if not goOn: continue

            # Adding in quadrature
            for k in h_sys.keys():
                print k
                for ib in range( 1 + h_rel_err.GetNbinsX() ):
                  print h_sys[k].GetBinContent(ib)
                  h_rel_err.SetBinContent(ib, h_rel_err.GetBinContent(ib) + h_sys[k].GetBinContent(ib)**2 )

            ## In case one wants to add uncertainties to specific backgrounds (like x-sec), that can be done here
            #if True:
            #    for ib in range(1 + h_rel_err.GetNbinsX() ):
            #        counts = [ shapeHists[x].GetBinContent(ib) for x in mc ]
            #        totalCount = sum(counts)
            #        print "Count in bin %s: %.2f"%(ib, totalCount)
            #        shapeUnc = [ 0 ]
            #        print rare_mc
            #        print rare_mc.name
            #        shapeUnc.append( (0.50*shapeHists[rare_mc].GetBinContent(ib))**2 )
            #        shapeUnc.append( (0.025*(totalCount))**2 )
            #        h_rel_err.SetBinContent(ib, h_rel_err.GetBinContent(ib) + sum( shapeUnc ) )

            # take sqrt
            for ib in range( 1 + h_rel_err.GetNbinsX() ):
                h_rel_err.SetBinContent(ib, sqrt( h_rel_err.GetBinContent(ib) ) )

            # Divide
            h_rel_err.Divide(h_summed[None])

            plot = plot_mc[None]
            if args.normalizeBinWidth: plot.name += "_normalizeBinWidth"
            if not args.noData:
                data_histo    = plot_data.histos[0][0]
                for h in plot_data.histos[0][1:]:
                    data_histo.Add(h)

                data_histo.style = styles.errorStyle( ROOT.kBlack )
                plot.histos += [[ data_histo ]]
                plot.stack += [[ plot_data.stack[0][0] ]]

            boxes = []
            ratio_boxes = []
            if not args.noData:
                for ib in range(1, 1 + h_rel_err.GetNbinsX() ):
                    val = h_summed[None].GetBinContent(ib)
                    if val<0: continue
                    sys = h_rel_err.GetBinContent(ib)
                    box = ROOT.TBox( h_rel_err.GetXaxis().GetBinLowEdge(ib),  max([0.03, (1-sys)*val]), h_rel_err.GetXaxis().GetBinUpEdge(ib), max([0.03, (1+sys)*val]) )
                    box.SetLineColor(ROOT.kBlack)
                    box.SetFillStyle(3444)
                    box.SetFillColor(ROOT.kBlack)
                    r_box = ROOT.TBox( h_rel_err.GetXaxis().GetBinLowEdge(ib),  max(0.1, 1-sys), h_rel_err.GetXaxis().GetBinUpEdge(ib), min(1.9, 1+sys) )
                    r_box.SetLineColor(ROOT.kBlack)
                    r_box.SetFillStyle(3444)
                    r_box.SetFillColor(ROOT.kBlack)

                    boxes.append( box )
                    ratio_boxes.append( r_box )

                    ratio = {'yRange':(0.1,1.9), 'drawObjects':ratio_boxes}
            else:
                ratio = None
            #print "plot.histos[0][pos_top].Integral()", pos_top,plot.histos 
            #print "plot.histos[0][pos_top].Integral()", plot.histos[0][pos_top].Integral()    
            for log in [False, True]:
                plotDir = os.path.join(plot_directory, 'systematicsPlots', str(args.year), args.plot_directory, mode + ("_log" if log else "") + "_scaled", args.selection)
                if args.showOnly: plotDir = os.path.join(plotDir, "only_" + args.showOnly)
                plotting.draw(plot,
                    plot_directory = plotDir,
                    ratio = ratio,
                    legend = (0.50,0.88-0.04*sum(map(len, plot.histos)),0.95,0.88),
                    logX = False, logY = log, sorting = True,
                    yRange = (0.03, "auto"),
                    #drawObjects = drawObjects( True, top_sf[None], lumi_scale ) + boxes,
                    drawObjects = drawObjects( True, 1, lumi_scale ) + boxes,
                    copyIndexPHP = True
                )

