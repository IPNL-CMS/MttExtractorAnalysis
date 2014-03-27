#! /usr/bin/env python

import os, datetime, pwd

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-j", "--process", action="store", dest="cores", type="int", default=1, help="Number of core to use for launching")
parser.add_option("", "--run", action="store_true", dest="run", default=False, help="run crab")
(options, args) = parser.parse_args()

datasets = [

    ["/ElectronHad/chassera-ElectronHad_Run2012A-22Jan2013_03May13-v1-a1aaf6ddd9733b6d1a75fb37a68d50db/USER	", "ElectronHad_Run2012A-22Jan2013", "FT53_V21A_AN6"],
    ["/SingleElectron/chassera-SingleElectron_Run2012B-TOPElePlusJets-22Jan2013_03May13-v1-a1aaf6ddd9733b6d1a75fb37a68d50db/USER", "SingleElectron_Run2012B-TOPElePlusJets-22Jan2013", "FT53_V21A_AN6"],
    ["/SingleElectron/chassera-SingleElectron_Run2012C-TOPElePlusJets-22Jan2013_03May13-v1-a1aaf6ddd9733b6d1a75fb37a68d50db/USER", "SingleElectron_Run2012C-TOPElePlusJets-22Jan2013", "FT53_V21A_AN6"],

    ["/MuHad/chassera-MuHad_Run2012A-22Jan2013_03May13-v1-a1aaf6ddd9733b6d1a75fb37a68d50db/USER", "MuHad_Run2012A-22Jan2013", "FT53_V21A_AN6"],
    ["/SingleMu/chassera-SingleMu_Run2012B-TOPMuPlusJets-22Jan2013_03May13-v1-a1aaf6ddd9733b6d1a75fb37a68d50db/USER", "SingleMu_Run2012B-TOPMuPlusJets-22Jan2013", "FT53_V21A_AN6"],
    ["/SingleMu/chassera-SingleMu_Run2012C-TOPMuPlusJets-22Jan2013_03May13-v1-a1aaf6ddd9733b6d1a75fb37a68d50db/USER", "SingleMu_Run2012C-TOPMuPlusJets-22Jan2013", "FT53_V21A_AN6"],
    ["/SingleMu/chassera-SingleMu_Run2012D-TOPMuPlusJets-22Jan2013_13May13-v1-a1aaf6ddd9733b6d1a75fb37a68d50db/USER", "SingleMu_Run2012D-TOPMuPlusJets-22Jan2013", "FT53_V21A_AN6"],

    ]

# Get email address
email = "%s@ipnl.in2p3.fr" % (pwd.getpwuid(os.getuid()).pw_name)

d = datetime.datetime.now().strftime("%d%b")

version = 1

print("Creating configs for crab. Today is %s, you are %s and it's version %d" % (d, email, version))
print("")

def processDataset(dataset):
    dataset_path = dataset[0]
    dataset_name = dataset[1]
    dataset_globaltag = dataset[2]

    ui_working_dir = ("crab_data_%s_%s") % (dataset_name, d)
    output_file = "crab_data_%s_%s.cfg" % (dataset_name, d)
    output_dir = ("HTT/Extracted/data/%s/%s" % (d, dataset_name))

    python_config = "";
    if "Electron" in dataset_path:
        python_config = "Extractor_MTT_semie.py"
    else:
        python_config = "Extractor_MTT_semimu.py"

    print("Creating config file for '%s'" % (dataset_path))
    print("\tName: %s" % dataset_name)
    print("\tGlobal tag: %s" % dataset_globaltag)
    print("\tOutput directory: %s" % output_dir)
    print("")

    os.system("sed -e \"s#@datasetname@#%s#\" -e \"s#@uiworkingdir@#%s#g\" -e \"s#@email@#%s#g\" -e \"s#@globaltag@#%s#g\" -e \"s#@pset@#%s#g\" -e \"s#@outputdir@#%s#\" crab_data.cfg.template.ipnl > %s" % (dataset_path, ui_working_dir, email, dataset_globaltag, python_config, output_dir, output_file))

    if options.run:
        cmd = "crab -create -submit -cfg %s" % (output_file)
        os.system(cmd)

import multiprocessing
pool = multiprocessing.Pool(options.cores)
pool.map(processDataset, datasets)
