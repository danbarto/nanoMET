# pT dependent 15 GeV jet thresh
python tune.py --year 2016 --ttbarModifier 1 --maxSig 25 --jetThreshold 15 --pTdependent >> /dev/null 2>>tune.err
python tune.py --year 2017 --ttbarModifier 1 --maxSig 25 --jetThreshold 15 --pTdependent >> /dev/null 2>>tune.err
python tune.py --year 2018 --ttbarModifier 1 --maxSig 25 --jetThreshold 15 --pTdependent >> /dev/null 2>>tune.err
python tune.py --year 2016 --runData         --maxSig 25 --jetThreshold 15 --pTdependent >> /dev/null 2>>tune.err
python tune.py --year 2017 --runData         --maxSig 25 --jetThreshold 15 --pTdependent >> /dev/null 2>>tune.err
python tune.py --year 2018 --runData         --maxSig 25 --jetThreshold 15 --pTdependent >> /dev/null 2>>tune.err

# pT dependent 25 GeV jet thresh
python tune.py --year 2016 --ttbarModifier 1 --maxSig 25 --jetThreshold 25 --pTdependent >> /dev/null 2>>tune.err
python tune.py --year 2017 --ttbarModifier 1 --maxSig 25 --jetThreshold 25 --pTdependent >> /dev/null 2>>tune.err
python tune.py --year 2018 --ttbarModifier 1 --maxSig 25 --jetThreshold 25 --pTdependent >> /dev/null 2>>tune.err
python tune.py --year 2016 --runData         --maxSig 25 --jetThreshold 25 --pTdependent >> /dev/null 2>>tune.err
python tune.py --year 2017 --runData         --maxSig 25 --jetThreshold 25 --pTdependent >> /dev/null 2>>tune.err
python tune.py --year 2018 --runData         --maxSig 25 --jetThreshold 25 --pTdependent >> /dev/null 2>>tune.err

# 15 GeV jet thresh
python tune.py --year 2016 --ttbarModifier 1 --maxSig 25 --jetThreshold 15 >> /dev/null 2>>tune.err
python tune.py --year 2017 --ttbarModifier 1 --maxSig 25 --jetThreshold 15 >> /dev/null 2>>tune.err
python tune.py --year 2018 --ttbarModifier 1 --maxSig 25 --jetThreshold 15 >> /dev/null 2>>tune.err
python tune.py --year 2016 --runData         --maxSig 25 --jetThreshold 15 >> /dev/null 2>>tune.err
python tune.py --year 2017 --runData         --maxSig 25 --jetThreshold 15 >> /dev/null 2>>tune.err
python tune.py --year 2018 --runData         --maxSig 25 --jetThreshold 15 >> /dev/null 2>>tune.err

# 25 GeV jet thresh
python tune.py --year 2016 --ttbarModifier 1 --maxSig 25 --jetThreshold 25 >> /dev/null 2>>tune.err
python tune.py --year 2017 --ttbarModifier 1 --maxSig 25 --jetThreshold 25 >> /dev/null 2>>tune.err
python tune.py --year 2018 --ttbarModifier 1 --maxSig 25 --jetThreshold 25 >> /dev/null 2>>tune.err
python tune.py --year 2016 --runData         --maxSig 25 --jetThreshold 25 >> /dev/null 2>>tune.err
python tune.py --year 2017 --runData         --maxSig 25 --jetThreshold 25 >> /dev/null 2>>tune.err
python tune.py --year 2018 --runData         --maxSig 25 --jetThreshold 25 >> /dev/null 2>>tune.err
