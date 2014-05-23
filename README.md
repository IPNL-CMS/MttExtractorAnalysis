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
Extractors> git checkout cmssw_5_3_18
Extractors> cd ..
src> scram b -j8
```

