""" Copy data tune parameters if you tuned MC with multiple ttbar modifiers
"""

import os
from shutil import copyfile

dataTag = "v1"

allTunes = os.listdir("results")
dataTunes = filter(lambda x: "data" in x.lower(), allTunes)
mcTunes   = filter(lambda x: "mc"   in x.lower(), allTunes)

for tune in mcTunes:
    if "v1mod" in tune:
        tag = [ t for t in tune.split("_") if "v1mod" in t ][0]
        oldDataFile = tune.replace(tag,dataTag).replace("MC", "Data", 1).replace("MC","DATA")
        newDataFile = tune.replace("MC", "Data", 1).replace("MC","DATA")
        if oldDataFile in dataTunes and not os.path.exists(os.path.join(os.path.abspath("results"),newDataFile)):
            fromPath = os.path.join(os.path.abspath("results"),oldDataFile)
            toPath   = os.path.join(os.path.abspath("results"),newDataFile)
            copyfile( fromPath, toPath )

