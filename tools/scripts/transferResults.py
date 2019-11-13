#!/usr/bin/python

import os

dummyFormat = "\t\t#{0}\n\t\t#{1}\n\t\t{2}: xxx\n\t\t\t'data': {3},\n\t\t\t'mc':   {4}\n\t\tyyy,\n"

filesData  = {}
filesMC    = {}
parSetData = {}
parSetMC   = {}
testKey    = []    

respath  = os.path.join(os.environ["CMSSW_BASE"], "src", "nanoMET", "run", "results" )
tunepath = os.path.join(os.environ["CMSSW_BASE"], "src", "nanoMET", "plots" )

for file in os.listdir( respath ):
        if not file.endswith(".txt"): continue
        if "_unc" in file: continue

        with open(os.path.join(respath,file), "r") as infile:
            params = infile.readlines()[0].split("\n")[0]

        fileName = file.split(".txt")[0]
        allSplit = fileName.split("_")

        isData   = "data" in fileName.lower()
        sel      = allSplit[6]
        pTdep    = "pTdep" in allSplit[-1]
        version  = allSplit[-1-int(pTdep)]
        year     = allSplit[-2-int(pTdep)]
        settings = "_".join(allSplit[7:-2-int(pTdep)])

        key = {"year":year, "selection":sel, "pTdependent":pTdep, "version":version, "settings":settings}
        if isData:
            parSetData.update( {frozenset(key.items()):params} )
            filesData.update( {frozenset(key.items()):fileName} )
        else:
            parSetMC.update( {frozenset(key.items()):params} )
            filesMC.update( {frozenset(key.items()):fileName} )

if os.path.exists( os.path.join( tunepath, "tunes.py" ) ):
    with open( os.path.join( tunepath, "tunes.py"), "r") as f:
        prevRes = f.readlines()

    # remove } from previous results file
    for i, line in enumerate(prevRes[::-1]):
        if "}" in line:
            prevRes = prevRes[:len(prevRes)-i-1]
            break
else:
    prevRes = ["\ntuneParams = {\n\n"]

# get new entries of results parameters
for key, fileName in filesData.iteritems():
    dataFile = fileName
    mcFile   = filesMC[key]
    dataList = parSetData[key]
    mcList   = parSetMC[key]
    key      = dict(key)
    sett     = "".join([ str(s) for s in key["settings"] if s.isdigit() ])
    try:
        fileKey = int("".join( [key["year"],"8" if key["pTdependent"] else "9", version.strip("v"), sett ] ))
    except:
        fileKey = key["year"]
    testKey.append(fileKey)
    prevRes.append( dummyFormat.format( dataFile,mcFile,fileKey,dataList,mcList ).replace("xxx","{").replace("yyy","}")+"\n")

prevRes.append("\t}\n")
with open( os.path.join( tunepath, "tunes.py"), "w") as f:
    for line in prevRes:
        f.write(line)

if len(testKey) == len(set(testKey)):
    print "Success: No double key entries produced!"
else:
    print "Attention: Double key entries produced!"
