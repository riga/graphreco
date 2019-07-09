# graphreco


### Setup

```shell
# variables
export SCRAM_ARCH="slc7_amd64_gcc700"
export CMSSW_VERSION="CMSSW_11_0_0_pre3"

# setup CMSSW
source "/cvmfs/cms.cern.ch/cmsset_default.sh" ""
scramv1 project CMSSW "$CMSSW_VERSION"
cd "$CMSSW_VERSION/src"
eval `scramv1 runtime -sh`

# setup the graphreco repo
git clone git@github.com:riga/graphreco.git RecoHGCal/GraphReco

# compile
scram b -j
```


### Test

```shell
cmsRun RecoHGCal/GraphReco/test/windowInference_cfg.py
```
