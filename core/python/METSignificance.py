'''
Calculate the MET significance based on the standard approach, using nanoAOD input.
'''

import ROOT
import math
import os

import FWCore.ParameterSet.Config as cms
from RecoMET.METProducers.METSignificanceParams_cfi import *

class METSignificance:

    def __init__(self, JERera = 'Spring16_25nsV10_MC'):
        jetCorrParam = ROOT.JetCorrectorParameters()

        self.JERera             = JERera
        self.JetResolutionFile  = "$CMSSW_BASE/src/JetMETCorrections/Modules/src/JetResolution.cc+"
        self.JERdirectory       = "$CMSSW_BASE/src/nanoMET/core/data/JER/JRDatabase/textFiles/%s/"%JERera
        self.JetResolutionFile  = os.path.expandvars(self.JetResolutionFile)
        ROOT.gROOT.ProcessLine('.L '+self.JetResolutionFile)
        
        # Load the resolutions for data and MC, pt and phi
        self.JERdirectory   = os.path.expandvars(self.JERdirectory)
        self.res_pt         = ROOT.JME.JetResolution("%s/%s_PtResolution_AK4PFchs.txt"%(self.JERdirectory, self.JERera))
        self.res_phi        = ROOT.JME.JetResolution("%s/%s_PhiResolution_AK4PFchs.txt"%(self.JERdirectory, self.JERera))
        self.jer_SF         = ROOT.JME.JetResolutionScaleFactor("%s/%s_SF_AK4PFchs.txt"%(self.JERdirectory, self.JERera))

    def getJER(self, event):
        # calculate JER for all jets in the events
        event.Jet_dpt  = []
        event.Jet_dphi = []
        for i in range(event.nJet):
            jet = ROOT.JME.JetParameters()
            jet.setJetEta(event.Jet_eta[i]).setJetPt(event.Jet_pt[i]).setRho(event.fixedGridRhoFastjetAll)
            event.Jet_dpt    += [ self.res_pt.getResolution(jet) ]
            event.Jet_dphi   += [ self.res_phi.getResolution(jet) ]
            
    
    def calculate(self, event, args):
        
        cov_xx  = 0
        cov_xy  = 0
        cov_yy  = 0
        jet_pt  = event.Jet_pt
        i = 0
        for j in jet_pt:
            j_pt = j
            j_phi = event.Jet_phi[i]
            j_sigmapt = event.Jet_dpt[i]
            j_sigmaphi = event.Jet_dphi[i]
            index = event.Jet_etabin[i]

            cj = math.cos(j_phi)
            sj = math.sin(j_phi)
            dpt = args[index] * j_pt * j_sigmapt
            dph =               j_pt * j_sigmaphi

            dpt *= dpt
            dph *= dph

            cov_xx += dpt*cj*cj + dph*sj*sj
            cov_xy += (dpt-dph)*cj*sj
            cov_yy += dph*cj*cj + dpt*sj*sj

            i += 1

        # unclustered energy
        cov_tt = args[5]*args[5] + args[6]*args[6]*event.MET_sumEt
        cov_xx += cov_tt
        cov_yy += cov_tt

        det = cov_xx*cov_yy - cov_xy*cov_xy

        if det>0:
            ncov_xx =  cov_yy / det
            ncov_xy = -cov_xy / det
            ncov_yy =  cov_xx / det
        else:
            #print cov_xx, cov_yy, cov_xy
            ncov_xx = cov_xx if cov_xx > 0 else 1
            ncov_yy = cov_yy if cov_yy > 0 else 1
            ncov_xy = cov_xy if cov_xy > 0 else 1

        met_x = event.MET_pt * math.cos(event.MET_phi)
        met_y = event.MET_pt * math.sin(event.MET_phi)

        event.det = det
        event.MET_sig = met_x*met_x*ncov_xx + 2*met_x*met_y*ncov_xy + met_y*met_y*ncov_yy
    
