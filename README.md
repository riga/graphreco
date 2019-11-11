# graphreco
working title: Castrop-Rauxel 


### Setup

```shell
# variables
export SCRAM_ARCH="slc7_amd64_gcc820"
export CMSSW_VERSION="CMSSW_11_0_0_pre9"

# setup CMSSW
source "/cvmfs/cms.cern.ch/cmsset_default.sh" ""
cmsrel $CMSSW_VERSION
cd "$CMSSW_VERSION/src"
cmsenv

# setup the graphreco repo
git clone https://github.com/jkiesele/graphreco.git RecoHGCal/GraphReco

# compile
scram b -j
```


### Test

```shell
cmsRun RecoHGCal/GraphReco/test/windowInference_cfg.py
```
