#!/usr/bin/python

import os

dummyFormat = "\t\t#{0}\n\t\t#{1}\n\t\t{2}: xxx\n\t\t\t'data':         {3},\n\t\t\t'mc':           {4},\n\t\t\t'jerData':      '{5}',\n\t\t\t'jerMC':        '{6}',\n\t\t\t'jetThreshold': {7},\n\t\t\t'year':         {8},\n\t\tyyy,\n"

year       = {}
jetThresh  = {}
jerData    = {}
jerMC      = {}
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
        pTdep    = "pTdep" in file
        version  = allSplit[-1]
        y        = allSplit[-2]
        settings = "_".join(allSplit[7:-2])

        jetThreshold = int(allSplit[7].split("sumPt")[1])
        jer          = "_".join(allSplit[3:6])

        key = {"year":y, "selection":sel, "pTdependent":pTdep, "version":version, "settings":settings, "jer":jer.replace("DATA","").replace("MC","")}

        year.update( {frozenset(key.items()):y} )
        jetThresh.update( {frozenset(key.items()):jetThreshold} )
        if isData:
            parSetData.update( {frozenset(key.items()):params} )
            filesData.update( {frozenset(key.items()):fileName} )
            jerData.update( {frozenset(key.items()):jer} )
        else:
            parSetMC.update( {frozenset(key.items()):params} )
            filesMC.update( {frozenset(key.items()):fileName} )
            jerMC.update( {frozenset(key.items()):jer} )

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

#    try:
    if True:
        jThresh  = jetThresh[key]
        dataFile = fileName
        mcFile   = filesMC[key]
        dataList = parSetData[key]
        mcList   = parSetMC[key]
        jData    = jerData[key]
        jMC      = jerMC[key]
        key      = dict(key)
        y        = int(key["year"])
        sett     = "".join([ str(s) for s in key["settings"] if s.isdigit() ])

#    except:
#        print key
#        continue

    index    = len( filter( lambda k: str(y) in str(k), testKey ) ) + 1
    fileKey  = "%i%i"%(y,index)

#    try:
#        fileKey = int("".join( [y,"8" if key["pTdependent"] else "9", version.strip("v"), sett ] ))
#    except:
#        fileKey = y

#    newFileKey = fileKey
#    i = 1
#    while newFileKey in testKey:
#        newFileKey = int(str(fileKey)+str(i))
#        i += 1
#    print newFileKey
    testKey.append(fileKey)
    prevRes.append( dummyFormat.format( dataFile,mcFile,fileKey,dataList,mcList,jData,jMC,jThresh,y ).replace("xxx","{").replace("yyy","}")+"\n")

prevRes.append("\t}\n")
with open( os.path.join( tunepath, "tunes.py"), "w") as f:
    for line in prevRes:
        f.write(line)

if len(testKey) == len(set(testKey)):
    print "Success: No double key entries produced!"
else:
    print "Attention: Double key entries produced!"
