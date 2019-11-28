'''

'''

import math

etabins = [0.8,1.3,1.9,2.5,100] # 3.1 added once

def getBin(abseta):
    for i, a in enumerate(etabins):
        if abseta < a:
            return int(i)
            break

def cartesian(pt, phi):
    return (pt*math.cos(phi), pt*math.sin(phi))

class Event:
    def __init__(self, event, jetResolution, weightModifier=1, METPtVar="MET_pt", METPhiVar="MET_phi", JetCollection="Jet_pt", isData=False, vetoEtaRegion=(10,10), jetThreshold=15., puWeight="puWeight", pTdepMetSig=True):
        jetResolution.getJER(event, JetCollection=JetCollection)

        self.nJet        = event.nJet
        self.pTdepMetSig = pTdepMetSig
        # The preliminary conclusion on jet/lepton cleaning is: use simple deltaR cleaning of jets against electrons/muons/photons
        cleanJetIndices  = [ i for i,x in enumerate(event.Jet_cleanmaskMETSig) if x>0 ]
        
        # only use jets above threshold and outside a veto region
        acceptedJetIndices  = [ x for x in cleanJetIndices if ( getattr(event, JetCollection)[x]>float(jetThreshold) and not (float(vetoEtaRegion[0]) <= abs(event.Jet_eta[x]) < float(vetoEtaRegion[1])) ) ]
        lowPtJetIndices     = [ x for x in cleanJetIndices if ( getattr(event, JetCollection)[x]<float(jetThreshold) and not (float(vetoEtaRegion[0]) <= abs(event.Jet_eta[x]) < float(vetoEtaRegion[1])) ) ] # EE stuff (if rejected) shouldn't be used for anything

        self.Jet_pt     = [ getattr(event, JetCollection)[i] for i in acceptedJetIndices ]
        self.Jet_eta    = [ event.Jet_eta[i]                 for i in acceptedJetIndices ]
        self.Jet_etabin = [ getBin(abs(x))                   for x in self.Jet_eta ]
        self.Jet_phi    = [ event.Jet_phi[i]                 for i in acceptedJetIndices ]
        self.Jet_dpt    = [ event.Jet_dpt[i]                 for i in acceptedJetIndices ]
        self.Jet_dphi   = [ event.Jet_dphi[i]                for i in acceptedJetIndices ]

        self.MET_pt     = getattr(event, METPtVar)
        self.MET_phi    = getattr(event, METPhiVar)
        self.MET_sumPt  = getattr(event, "MET_sumPt")
        
        sumPt_lowPtJet  = sum( [ getattr(event, JetCollection)[x] for x in lowPtJetIndices ] )
        self.MET_sumPt += sumPt_lowPtJet

        self.fixedGridRhoFastjetAll = event.fixedGridRhoFastjetAll
        if puWeight:
            self.weight = event.weight * getattr(event, puWeight) * weightModifier if not isData else weightModifier
        else:
            self.weight = event.weight * weightModifier if not isData else weightModifier

    def calcLL(self, args):
        # calculate the log likelihood
        self.calcMETSig(args)

        try:
            LL = self.weight * ( self.MET_sig + math.log(self.det) )
        except:
            LL = self.weight * self.MET_sig

        return LL

    def calcMETSig(self, args):
        
        cov_xx  = 0
        cov_xy  = 0
        cov_yy  = 0
        jet_pt  = self.Jet_pt

        for i,j in enumerate(jet_pt):
            j_pt        = j
            j_phi       = self.Jet_phi[i]
            j_sigmapt   = self.Jet_dpt[i]
            j_sigmaphi  = self.Jet_dphi[i]
            index       = self.Jet_etabin[i]
            jet_index   = 0 if j_pt < 40 else 1

            cj      = math.cos(j_phi)
            sj      = math.sin(j_phi)
            dpt     = args[2*index + jet_index if self.pTdepMetSig else index] * j_pt * j_sigmapt
            dph     = j_pt * j_sigmaphi

            dpt    *= dpt
            dph    *= dph

            cov_xx += dpt*cj*cj + dph*sj*sj
            cov_xy += (dpt-dph)*cj*sj
            cov_yy += dph*cj*cj + dpt*sj*sj
            

        # unclustered energy
        ind    = 10 if self.pTdepMetSig else 5
        cov_tt = args[ind]*args[ind] + args[ind+1]*args[ind+1]*self.MET_sumPt
        cov_xx += cov_tt
        cov_yy += cov_tt

        det = cov_xx*cov_yy - cov_xy*cov_xy

        if det>0:
            ncov_xx =  cov_yy / det
            ncov_xy = -cov_xy / det
            ncov_yy =  cov_xx / det
        else:
            ncov_xx = cov_xx if cov_xx > 0 else 1
            ncov_yy = cov_yy if cov_yy > 0 else 1
            ncov_xy = cov_xy if cov_xy > 0 else 1

        met_x = self.MET_pt * math.cos(self.MET_phi)
        met_y = self.MET_pt * math.sin(self.MET_phi)


        self.det = det
        self.MET_sig = met_x*met_x*ncov_xx + 2*met_x*met_y*ncov_xy + met_y*met_y*ncov_yy

