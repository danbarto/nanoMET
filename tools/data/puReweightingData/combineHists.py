import ROOT

def getObjFromFile(fname, hname):
    gDir = ROOT.gDirectory.GetName()
    f = ROOT.TFile(fname)
    assert not f.IsZombie()
    f.cd()
    htmp = f.Get(hname)
    if not htmp:  return htmp
    ROOT.gDirectory.cd('PyROOT:/')
    res = htmp.Clone()
    f.Close()
    ROOT.gDirectory.cd(gDir+':/')
    return res

def writeObjsToFile(fname, objs):
    gDir = ROOT.gDirectory.GetName()
    f = ROOT.TFile(fname, 'recreate')
    for obj in objs:
        objw = obj.Clone()
        objw.Write()
    f.Close()
    ROOT.gDirectory.cd(gDir+':/')
    return

h_central   = getObjFromFile("PU_2018_58830_XSecCentral.root", "pileup")
h_plus      = getObjFromFile("PU_2018_58830_XSecUp.root", "pileup")
h_minus     = getObjFromFile("PU_2018_58830_XSecDown.root", "pileup")

h_plus.SetName("pileup_plus")
h_minus.SetName("pileup_minus")

writeObjsToFile("pileup_Cert_314472-325175_13TeV_PromptReco_Collisions18_withVar.root", [h_central, h_plus, h_minus])

