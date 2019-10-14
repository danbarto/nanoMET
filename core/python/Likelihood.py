import ROOT

class minLL( ROOT.TPyMultiGenFunction ):
    def __init__( self, evlist ):
        ROOT.TPyMultiGenFunction.__init__( self, self )
        self.evlist = evlist
    def NDim( self ):
        return 12 # changed from 7 --> 12

    def DoEval( self, args ):
        LL = self.evlist.getLL(args)
        return LL
