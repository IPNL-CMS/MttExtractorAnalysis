# Presentation

This is the mtt Extractor analysis.

## Get the code

First, you need to retrieve PatExtractor. For that, refer to https://github.com/IPNL-CMS/PatExtractor

```bash
src> mkdir Extractors
src> cd Extractors
Extractors> git clone https://github.com/IPNL-CMS/MttExtractorAnalysis.git
Extractors> cd ..
src> export CVSROOT=":ext:<cern-user-account>@lxplus5.cern.ch:/afs/cern.ch/user/c/cvscmssw/public/CMSSW"
src> addpkg RecoBTag/PerformanceDB V01-03-21
src> scram b -j8
```

