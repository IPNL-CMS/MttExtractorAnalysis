# Presentation

This is the mtt Extractor analysis.

## Get PatExtractor

First, you need to install PatExtractor. For that, follow this install guide: https://github.com/IPNL-CMS/PatExtractor/blob/cmssw_5_3_18/INSTALL.md

## Installation

Execute the following commands to install MttExtractorAnalysis

### CMSSW dependencies

```bash
src> git cms-addpkg TopQuarkAnalysis/TopObjectResolutions
```

### MttExtractorAnalysis
```bash
src> cd Extractors
Extractors> git clone https://github.com/IPNL-CMS/MttExtractorAnalysis.git
Extractors> cd MttExtractorAnalysis
MttExtractorAnalysis> git checkout cmssw_5_3_18
MttExtractorAnalysis> cd ..
Extractors> cd ..
```

### Apply patches
```bash
src> git am Extractors/MttExtractorAnalysis/patches/0001-Fix-resolution-creator.patch
src> git am Extractors/MttExtractorAnalysis/patches/0001-Update-top-objects-resolution-using-8-TeV-POWHEG-TT-.patch
```

### Build
```bash
src> scram b -j8
```
