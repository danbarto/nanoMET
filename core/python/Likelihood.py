import ROOT

class minLL( ROOT.TPyMultiGenFunction ):
    def __init__( self, evlist, pTdepMetSig=True ):
        ROOT.TPyMultiGenFunction.__init__( self, self )
        self.evlist      = evlist
        self.pTdepMetSig = pTdepMetSig
    def NDim( self ):
        return 12 if self.pTdepMetSig else 7

    def DoEval( self, args ):
        LL = self.evlist.getLL(args)
        return LL
