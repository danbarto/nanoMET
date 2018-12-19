'''
Calculate the MET significance based on the standard approach, using nanoAOD input.
'''

import ROOT
import math
import os

import FWCore.ParameterSet.Config as cms
from RecoMET.METProducers.METSignificanceParams_cfi import *

class JetResolution:

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
            
    
