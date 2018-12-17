'''

'''

import math

from nanoMET.core.METSignificance import METSignificance

MS = METSignificance()

etabins = [0.8,1.3,1.9,2.5,100]

def getBin(abseta):
    for i, a in enumerate(etabins):
        if abseta < a:
            return int(i)
            break

def cartesian(pt, phi):
    return (pt*math.cos(phi), pt*math.sin(phi))

class Event:
    def __init__(self, event, weightModifier=1):
        #self.MS = METSignificance()
        MS.getJER(event)

        self.nJet       = event.nJet
        # use the cleanmask from nanoAOD to clean against leptons, but don't use jetId as of now
        #cleanJetIndices = [ i for i,x in enumerate(event.Jet_cleanmask) if (x>0 and event.Jet_jetId[i]>=3) ]
        cleanJetIndices = [ i for i,x in enumerate(event.Jet_cleanmask) if x>0 ]
        self.Jet_pt     = [ event.Jet_pt[i]     for i in cleanJetIndices ]
        self.Jet_eta    = [ event.Jet_eta[i]    for i in cleanJetIndices ]
        self.Jet_etabin = [ getBin(abs(x))      for x in self.Jet_eta ]
        self.Jet_phi    = [ event.Jet_phi[i]    for i in cleanJetIndices ]
        self.Jet_dpt    = [ event.Jet_dpt[i]    for i in cleanJetIndices ]
        self.Jet_dphi   = [ event.Jet_dphi[i]   for i in cleanJetIndices ]

        self.MET_pt             = event.MET_pt
        self.MET_phi            = event.MET_phi
        self.MET_sumEt          = event.MET_sumEt
        self.MET_significance   = event.MET_significance
        


        self.fixedGridRhoFastjetAll = event.fixedGridRhoFastjetAll
        self.weight = event.weight * weightModifier

    def METSignificance(self, args):

        MS.calculate(self, args)

        try:
            LL = self.weight * ( self.MET_sig + math.log(self.det) )
        except:
            LL = self.weight * self.MET_sig

        return LL


