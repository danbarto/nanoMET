'''

'''

# Standard imports
import ROOT
import os
import sys
import math
import copy
import itertools
import json

from RootTools.core.standard import *

from nanoMET.core.Event import Event
from nanoMET.core.Likelihood import minLL

class run:

    def __init__(self, samples, selection, outfile="results/tune"):
        # Need fill a list in order to do the minimization, reading from the tree is too slow

        self.eventlist = []
        self.variables = map( TreeVariable.fromString,  ['weight/F', 'nJet/I', 'fixedGridRhoFastjetAll/F', 'MET_pt/F', 'MET_phi/F', 'MET_sumEt/F', 'MET_significance/F'] )
        self.variables += [VectorTreeVariable.fromString('Jet[pt/F,eta/F,phi/F]' ) ]
        self.outfile = outfile

        for s in samples:
            s.setSelectionString(selection)
            
            nEvents = s.getYieldFromDraw()
            print "Running over %s events for sample %s."%(nEvents['val'], s.name)

            reader = s.treeReader(variables=self.variables)
            reader.start()
            i = 0

            while reader.run():
                self.eventlist += [Event(reader.event)]
                if i%10000==0:
                    print "Done with %s"%i
                i += 1


    def getLL(self, args):
        # recalculate MET Significance, Determinant and LL for given parameters
        # although it uses C methods it's still rather slow
        LLs = map(lambda x: x.METSignificance(args), self.eventlist)
        LL = sum(LLs)
        
        return LL


    def minimize(self, start=[1.0, 1.0, 1.0, 1.0, 1.0, 0., .5], step=[0.05]*7, maxSig = 9):
        gmin = ROOT.Math.Factory.CreateMinimizer("Minuit2")
        gmin.SetTolerance(10.0)
        gmin.SetStrategy(0)
        gmin.SetPrintLevel(3)
        
        variable  = ['a1','a2','a3','a4','a5','u1','u2']
        step  = step
        start = start

        LL = minLL(self)

        gmin.SetFunction(LL)

        print 'Minimizing parameters',variable
        print 'With stepsize of', step
        print 'And starting values', start
    
        for i in range(0,7):
            gmin.SetVariable(i,variable[i],start[i], step[i])
    
        gmin.Minimize()

        #filter events with high significance
        self.eventlist = [x  for x in self.eventlist if x.MET_sig < maxSig ]
        
        print
        print 'Now fitting after applying significance cut'
        print 'Total events:',len(self.eventlist)
        
        gmin.SetStrategy(1)
        gmin.Minimize()
        gmin.Hesse()
        
        pars = [gmin.X()[i] for i in range(0,7)]
        uncs = [gmin.Errors()[i] for i in range(0,7)]
        
        with open(self.outfile+'.txt', 'w') as of:
            json.dump(pars, of)
        with open(self.outfile+'_unc.txt', 'w') as of:
            json.dump(uncs, of)
    



if __name__ == '__main__':

    from StopsDilepton.samples.nanoTuples_Summer16_postProcessed import *

    preselection = "nJetGood > 0 && nGoodMuons==2 && nGoodElectrons==0 && l1_pt > 25 && l2_pt > 20 && abs(dl_mass-91.2)<10"

    DY_LO.reduceFiles( to = 1 )

    r = run([DY_LO], preselection)

    r.getLL( [1.0, 1.0, 1.0, 1.0, 1.0, 0., .5] )
    
    minimize = True
    if minimize:

        r.minimize()
