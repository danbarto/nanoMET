''' 



'''
# Standrd imports
from math import sqrt, atan2, cos, sin, atan, exp

# RootTools 
from RootTools.core.standard import *

# Logging
import JetMET.tools.logger as logger
logger    = logger.get_logger('INFO', logFile = None)
import RootTools.core.logger as logger_rt
logger_rt = logger_rt.get_logger('INFO', logFile = None)

# JECconfig
from JetMET.JetCorrector.JetCorrector import JetCorrector, correction_levels_data, correction_levels_mc

# PANDAS
import pandas as pd

# 
def getDictFromObj(obj):
    return {'pt':obj.pt(), 'eta':obj.eta(), 'phi':obj.phi(), 'px': obj.px(), 'py': obj.py(), 'pdgId':obj.pdgId(), 'charge':obj.charge(), 'energy':obj.energy()}

def getTupleFromObj(obj):
    try:
        #print obj.pt(), obj.eta(), obj.phi()
        return (obj.pt(), obj.eta(), obj.phi(), obj.px(), obj.py(), obj.pdgId(), obj.charge(), obj.energy())
    except:
        #print "No footprint found"
        return (-99, -99, -99, -99, -99, -99, -99, -99)

def getDictFromTuple(obj):
    return {'pt':obj[0], 'eta':obj[1], 'phi':obj[2], 'px': obj[3], 'py': obj[4], 'pdgId':obj[5], 'charge':obj[6], 'energy':obj[7]}

## variations
def getUnclusteredDelta(cand, correlated=True, debug=False):
    ## 
    #res_x, res_y = 0, 0
    res_x, res_y = {}, {}
    for source in ['HF', 'charged', 'neutral', 'photon']:
        res_x[source], res_y[source] = 0,0
    for c in cand:
        # variations similar to 
        source = None
        # HF
        if abs(c['pdgId'])==1 or abs(c['pdgId'])==2:
            source = 'HF'
            res = sqrt( 1./c['energy'] + 0.0025 )
        # charged
        elif abs(c['pdgId'])==211: # this was a BS test or abs(c['pdgId'])==13 or abs(c['pdgId'])==11:
            source = 'charged'
            res = sqrt( (0.00009*c['energy'])**2 + ( 0.0085/(sin(2*atan(exp(-c['eta'])))) )**2 )
        # neutral
        elif abs(c['pdgId'])==130:
            source = 'neutral'
            res = min([0.25,sqrt(0.64/c['energy']+0.0025)]) if abs(c['eta'])<1.3 else min([0.30,sqrt(1.0/c['energy']+0.0016)])
        # photon
        elif abs(c['pdgId'])==22:
            source = 'photon'
            res = sqrt(0.0009/c['energy']+0.000001)
        
        if source:
            res_x[source] += (-1)*res*cos(c['phi'])*c['pt']
            res_y[source] += (-1)*res*sin(c['phi'])*c['pt']
    
    if debug:
        print "x variations", res_x['HF'], res_x['charged'], res_x['neutral'], res_x['photon']
        print "y variations", res_y['HF'], res_y['charged'], res_y['neutral'], res_y['photon']

    if correlated:
        return res_x['HF'] + res_x['charged'] + res_x['neutral'] + res_x['photon'], res_y['HF'] + res_y['charged'] + res_y['neutral'] + res_y['photon']
    else:
        return sqrt(res_x['HF']**2 + res_x['charged']**2 + res_x['neutral']**2 + res_x['photon']**2), sqrt(res_y['HF']**2 + res_y['charged']**2 + res_y['neutral']**2 + res_y['photon']**2)

## Run2018D
#corrector_old     = JetCorrector.fromTarBalls( [(1, 'Fall17_17Nov2017F_V32_DATA') ], correctionLevels = correction_levels_data ) # 
# from DAS
corrector     = JetCorrector.fromTarBalls( [(1, 'Autumn18_V19_MC') ], correctionLevels = correction_levels_mc )
corrector_L1  = JetCorrector.fromTarBalls( [(1, 'Autumn18_V19_MC') ], correctionLevels = ['L1FastJet'] )
# new JEC
#corrector     = JetCorrector.fromTarBalls( [(1, 'Autumn18_RunD_V8_DATA') ], correctionLevels = correction_levels_data )
#corrector_L1  = JetCorrector.fromTarBalls( [(1, 'Autumn18_RunD_V8_DATA') ], correctionLevels = ['L1FastJet'] )

#corrector = JetCorrector.fromTarBalls( [(1, 'Fall17_17Nov2017B_V6_DATA') ], correctionLevels = correction_levels_data ) #this is what's used in prompt reco
#corrector_new     = JetCorrector.fromTarBalls( [(1, 'Fall17_17Nov2017B_V32_DATA') ], correctionLevels = correction_levels_data )

# products  
products = { 
    'jets':{'type':'vector<pat::Jet>', 'label': ( "slimmedJets" )},
    'mets':{'type':'vector<pat::MET>',  'label':( "slimmedMETs" )},
    'rho': {'type':'double', 'label': ("fixedGridRhoFastjetAll")},
    'pf':  {'type':'vector<pat::PackedCandidate>', 'label': ("packedPFCandidates")}, 
    'ele': {'type':'vector<pat::Electron>', 'label':("slimmedElectrons")},
    'mu':  {'type':'vector<pat::Muon>', 'label':("slimmedMuons")},
    'pho': {'type':'vector<pat::Photon>', 'label':("slimmedPhotons")},
    'tau': {'type':'vector<pat::Tau>', 'label':("slimmedTaus")},
   }

# sample  
#sample = FWLiteSample.fromFiles("test", files = ['file:/afs/hephy.at/data/dspitzbart01/event_0.root'])
sample = FWLiteSample.fromFiles("test", files = ['file:./corMETMiniAOD.root'])
#sample = FWLiteSample.fromFiles("test", files = ['file:/afs/hephy.at/work/r/rschoefbeck/CMS/tmp/CMSSW_10_2_9/src/StopsDilepton/plots/plotsRobert/tails/2018_mumu_lepSel-POGMetSig12-njet1-btag1p-relIso0.12-looseLeptonVeto-mll20-dPhiJet0-dPhiJet1/tmp_DoubleMuon_Run2018D-PromptReco-v2_MINIAOD.root'])

#2018D
#sample = FWLiteSample.fromFiles("test", files = ['root://cms-xrd-global.cern.ch//store/data/Run2018D/DoubleMuon/MINIAOD/PromptReco-v2/000/322/106/00000/DC6DE7DE-ACB2-E811-815C-FA163E0F815B.root'])
#2017 test
#sample = FWLiteSample.fromFiles("test", files = ['file:dm_2017B_31Mar.root'] )
# Event loop
r = sample.fwliteReader( products = products )
r.start()
while r.run():

    #try:
    #    if not r.run(): break
    #except:
    #    print "Smth wrong in this event!", r.evt
    #    r.position+=1
    #    continue

    footprints = []
    footprints_unique = set()
    candidates = []
    candidates_unique = set()
    jets = []

    # footprints of clustered particles, or leptons/photons
    # add footprints from jets
    for j in r.event.jets:

        # only add footprints from jets with pt>15 GeV
        if j.pt()>15:
            print "Found a jet for MET"
            #jets.append(getDictFromObj(j))
            for i in range(j.numberOfSourceCandidatePtrs()):
                #footprints.append(getDictFromObj(j.sourceCandidatePtr(i)))
                footprints_unique.add(getTupleFromObj(j.sourceCandidatePtr(i)))
        else:
            print "Jet below 15 GeV"

    nJetFootPrints = len(footprints_unique)
    print "Added %s footprints from jets"%nJetFootPrints

    # add footprints from electrons
    for e in r.event.ele:
        for i in range(e.numberOfSourceCandidatePtrs()):
            footprints_unique.add(getTupleFromObj(e.sourceCandidatePtr(i)))

    nEleFootPrints = len(footprints_unique) - nJetFootPrints
    print "Added %s footprints from electrons"%nEleFootPrints

    # add footprints from muons
    for m in r.event.mu:
        for i in range(m.numberOfSourceCandidatePtrs()):
            footprints_unique.add(getTupleFromObj(m.sourceCandidatePtr(i)))

    nMuFootPrints = len(footprints_unique) - nJetFootPrints - nEleFootPrints
    print "Added %s footprints from muons"%nMuFootPrints
    
    # add footprints from photons
    for p in r.event.pho:
        for i in range(p.numberOfSourceCandidatePtrs()):
            footprints_unique.add(getTupleFromObj(p.sourceCandidatePtr(i)))

    nPhoFootPrints = len(footprints_unique) - nJetFootPrints - nEleFootPrints - nMuFootPrints
    print "Added %s footprints from photons"%nPhoFootPrints

    # add footprints from taus
    for t in r.event.tau:
        for i in range(t.numberOfSourceCandidatePtrs()):
            footprints_unique.add(getTupleFromObj(t.sourceCandidatePtr(i)))

    nTauFootPrints = len(footprints_unique) - nJetFootPrints - nEleFootPrints - nMuFootPrints -nPhoFootPrints
    print "Added %s footprints from taus"%nTauFootPrints

    print "Total footprints: %s"%len(footprints_unique)
    # that's it

    for c in r.event.pf:
        candidates.append(getDictFromObj(c))
        candidates_unique.add(getTupleFromObj(c))
    
    # all PF candidates - used to calculate RAW MET
    can_all = pd.DataFrame( candidates )
    
    #print "Raw MET:", sqrt(can_all['py'].sum()**2 + can_all['px'].sum()**2), atan2(-can_all['py'].sum(), -can_all['px'].sum())

    #print "Raw MET, sample:", r.event.mets[0].pt(), r.event.mets[0].phi()

    print "Total candidates: %s"%len(candidates_unique)
    candidates_unclustered = candidates_unique - footprints_unique
    
    print "Total unclustered candidates: %s"%len(candidates_unclustered)

    cand_uncl = [ getDictFromTuple(c) for c in candidates_unclustered ]

    print "raw MET_x:", can_all['px'].sum()
    print "raw MET_y:", can_all['py'].sum()
    print "raw MET_pt:", sqrt(can_all['py'].sum()**2+can_all['px'].sum()**2)

    can = pd.DataFrame( cand_uncl )
    print "unclustered MET_x:", can['px'].sum()
    print "unclustered MET_y:", can['py'].sum()
    print "unclustered MET_pt:", sqrt(can['py'].sum()**2+can['px'].sum()**2)

    
    print "Delta x:", getUnclusteredDelta(cand_uncl)[0]
    print "Delta y:", getUnclusteredDelta(cand_uncl)[1]

    print "Delta x, uncorr:", getUnclusteredDelta(cand_uncl, correlated=False)[0]
    print "Delta y, uncorr:", getUnclusteredDelta(cand_uncl, correlated=False)[1]

    #jet = pd.DataFrame( jets )

    

    #foo = pd.DataFrame(footprints)

    #intersection  = pd.merge(can, foo)

    shift_x, shift_y = 0., 0.
    shift_x_noL1, shift_y_noL1 = 0., 0.
    for j in r.event.jets:
#        print "pt (mAOD) % 7.5f pt(Uncorr) % 7.5f pT recorr(new) % 7.5f pT recorr(old,DAS) % 7.5f" %( 
#                j.pt(), 
#                j.correctedJet("Uncorrected").pt(),
#                j.correctedJet("Uncorrected").pt()*corrector_new.correction( j.correctedJet("Uncorrected").pt(), j.eta(), j.jetArea(), r.event.rho[0], r.event.run ), 
##                j.correctedJet("Uncorrected").pt()*corrector_old.correction( j.correctedJet("Uncorrected").pt(), j.eta(), j.jetArea(), r.event.rho[0], r.event.run ), 
#                j.correctedJet("Uncorrected").pt()*corrector.correction( j.correctedJet("Uncorrected").pt(), j.eta(), j.jetArea(), r.event.rho[0], r.event.run ), 
#            )
#        print j.pt(), j.correctedJet("Uncorrected").pt(), j.electronEnergyFraction(), j.photonEnergyFraction()

        if j.electronEnergyFraction()+j.photonEnergyFraction()<0.9:

            pt_raw   = j.correctedJet("Uncorrected").pt() 
            mu_ef    = j.muonEnergyFraction()
            phi      = j.phi()
            #jec      = corrector_new.correction( (1-mu_ef)*pt_raw, j.eta(), j.jetArea(), r.event.rho[0], r.event.run )
            jec      = corrector.correction( pt_raw, j.eta(), j.jetArea(), r.event.rho[0], r.event.run )
            jec_L1   = corrector_L1.correction( pt_raw, j.eta(), j.jetArea(), r.event.rho[0], r.event.run )
            #shift_x += j.correctedJet("Uncorrected").p)t()*cos(j.correctedJet("Uncorrected").phi()) - j.pt()*cos(j.phi())
            #shift_y += j.correctedJet("Uncorrected").pt()*sin(j.correctedJet("Uncorrected").phi()) - j.pt()*sin(j.phi())
            cor_pt        = mu_ef*pt_raw + (jec)*(1-mu_ef)*pt_raw
            cor_pt_L1     = mu_ef*pt_raw + (jec_L1)*(1-mu_ef)*pt_raw
            if jec*pt_raw*(1-mu_ef)>15:
                shift_x  += (pt_raw - cor_pt)*cos(phi)
                shift_y  += (pt_raw - cor_pt)*sin(phi)
                shift_x_noL1  += (cor_pt_L1 - cor_pt)*cos(phi)
                shift_y_noL1  += (cor_pt_L1 - cor_pt)*sin(phi)
                print "CORRECT WITH: jet pt (raw) % 7.5f muEF % 7.5f jec % 7.5f cor_pt % 7.5f cor_pt_L1 % 7.5f jec % 7.5f jec_L1 % 7.5f" %( pt_raw, mu_ef, jec, cor_pt, cor_pt_L1, jec, jec_L1)
            else:
                print "DON'T CORRECT WITH: jet pt (raw) % 7.5f muEF % 7.5f jec % 7.5f cor_pt % 7.5f cor_pt_L1 % 7.5f jec % 7.5f jec_L1 % 7.5f" %( pt_raw, mu_ef, jec, cor_pt, cor_pt_L1, jec, jec_L1)
    m = r.event.mets[0]
   
    met_type1_px = m.uncorPt()*cos(m.uncorPhi()) + shift_x
    met_type1_py = m.uncorPt()*sin(m.uncorPhi()) + shift_y
    met_type1 = sqrt( met_type1_px**2 + met_type1_py**2)

    met_type1_noL1_px = m.uncorPt()*cos(m.uncorPhi()) + shift_x_noL1
    met_type1_noL1_py = m.uncorPt()*sin(m.uncorPhi()) + shift_y_noL1
    met_type1_noL1 = sqrt( met_type1_px**2 + met_type1_py**2)

    MEx = sum( [ -p.px() for p in r.event.pf ] )   
    MEy = sum( [ -p.py() for p in r.event.pf ] )   
    rawMET_recalc = sqrt( MEx**2 + MEy**2 ) 
    print "%i:%i:%i MET(raw) % 7.5f MET(raw, recalc) % 7.5f MET (mAOD) % 7.5f MET (uncorr) % 7.5f myMET % 7.5f myMET (t1,noL1) % 7.5f" % (r.evt[0], r.evt[1], r.evt[2], m.uncorPt(), rawMET_recalc, m.pt(), m.uncorPt(), met_type1, met_type1_noL1)
