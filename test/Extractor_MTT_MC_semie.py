import sys, os
sys.path.insert(1, os.getcwd())

from Extractor_MTT_common import *

process = createExtractorProcess(True, False, useShiftCorrectedMET = True, globalTag = "PHYS14_25_V2")

process.source.fileNames = cms.untracked.vstring(
        'file:patTuple_1.root'
    )
