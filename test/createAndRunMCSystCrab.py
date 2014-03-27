#! /usr/bin/env python

import os, copy, datetime, pwd, re

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-j", "--process", action="store", dest="cores", type="int", default=1, help="Number of core to use for launching")
parser.add_option("", "--create-cfg", action="store_true", dest="create_cfg", default=False, help="create config files for crab")
parser.add_option("", "--run", action="store_true", dest="run", default=False, help="run crab")
parser.add_option("", "--status", action="store_true", dest="status", default=False, help="run crab -status")
parser.add_option("", "--get", action="store_true", dest="get", default=False, help="run crab -get")
parser.add_option("", "--resubmit", action="store_true", dest="resubmit", default=False, help="run crab -resubmit bad")
parser.add_option("", "--submit", action="store_true", dest="submit", default=False, help="run crab -submit all")
(options, args) = parser.parse_args()

if options.run:
  options.create_cfg = True

datasets = [

    # TT + jets
    ["/TTJets_matchingup_TuneZ2star_8TeV-madgraph-tauola/chassera-TTJets_matchingup_START53_V7A_02Oct13-v1-d7aa126bba28d24db7dfa79d30ae2a7a/USER", "TTJets_matchingup"],
    ["/TTJets_matchingdown_TuneZ2star_8TeV-madgraph-tauola/chassera-TTJets_matchingdown_START53_V7A_02Oct13-v1-d7aa126bba28d24db7dfa79d30ae2a7a/USER", "TTJets_matchingdown"],
    ["/TTJets_scaleup_TuneZ2star_8TeV-madgraph-tauola/chassera-TTJets_scaleup_START53_V7A_02Oct13-v1-d7aa126bba28d24db7dfa79d30ae2a7a/USER", "TTJets_scaleup"],
    ["/TTJets_scaledown_TuneZ2star_8TeV-madgraph-tauola/chassera-TTJets_scaledown_START53_V7A_02Oct13-v1-d7aa126bba28d24db7dfa79d30ae2a7a/USER", "TTJets_scaledown"],
    
    # Z + jets
    ["/DYJetsToLL_M-50_scaleup_8TeV-madgraph-tauola/sbrochet-DYJetsToLL_M-50_scaleup_START53_V7A_16Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "DYJetsToLL_M-50_scaleup"],
    ["/DYJetsToLL_M-50_scaledown_8TeV-madgraph-tauola/sbrochet-DYJetsToLL_M-50_scaledown_START53_V7A_16Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "DYJetsToLL_M-50_scaledown"],
    ["/DYJetsToLL_M-50_matchingup_8TeV-madgraph-tauola/sbrochet-DYJetsToLL_M-50_matchingup_START53_V7A_16Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "DYJetsToLL_M-50_matchingup"],
    ["/DYJetsToLL_M-50_matchingdown_8TeV-madgraph/sbrochet-DYJetsToLL_M-50_matchingdown_START53_V7A_16Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "DYJetsToLL_M-50_matchingdown"],

    # W + jets
    ["/WJetsToLNu_scaleup_8TeV-madgraph-tauola/sbrochet-WJetsToLNu_scaleup_START53_V7A_16Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "WJetsToLNu_scaleup"],
    ["/WJetsToLNu_scaledown_8TeV-madgraph-tauola/sbrochet-WJetsToLNu_scaledown_START53_V7A_16Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "WJetsToLNu_scaledown"],
    ["/WJetsToLNu_matchingup_8TeV-madgraph-tauola/sbrochet-WJetsToLNu_matchingup_START53_V7A_16Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "WJetsToLNu_matchingup"],
    ["/WJetsToLNu_matchingdown_8TeV-madgraph-tauola/sbrochet-WJetsToLNu_matchingdown_START53_V7A_16Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "WJetsToLNu_matchingdown"]

    ]



# Get email address
email = "%s@ipnl.in2p3.fr" % (pwd.getpwuid(os.getuid()).pw_name)

d = datetime.datetime.now().strftime("%d%b%y")

version = 1

print("Creating configs for crab. Today is %s, you are %s and it's version %d" % (d, email, version))
print("")

from string import Template
multicrab = Template(r"""[MULTICRAB]
cfg=crab_MC.cfg.template.ipnl

[COMMON]
USER.ui_working_dir = ${ui_working_dir}
USER.eMail = ${email}
CMSSW.datasetpath = ${dataset}
CMSSW.total_number_of_events = ${events}
"""
)

multicrab_semie = r"""
[${name}_semie]
CMSSW.pset = Extractor_MTT_MC_semie.py
USER.user_remote_dir = ${remote_dir_semie}
"""

multicrab_semimu = r"""
[${name}_semimu]
CMSSW.pset = Extractor_MTT_MC_semimu.py
USER.user_remote_dir = ${remote_dir_semimu}
"""

def processDataset(dataset):
    dataset_name = dataset[1]
    dataset_path = dataset[0]
    dataset_size = -1
    if len(dataset) > 2:
      dataset_size = dataset[2]
    #dataset_globaltag = re.search('START\d{0,2}_V\d[A-Z]?', dataset_path).group(0)

    #publish_name = "%s_%s_%s-v%d" % (dataset_name, dataset_globaltag, d, version)
    output_file = "multicrab_MC_%s_%s.cfg" % (dataset_name, d)
    ui_working_dir = ("multicrab_MC_%s") % (dataset_name)

    if options.create_cfg:
        output_dir_semie = ("HTT/Extracted/Systematics/%s/semie/%s" % (d, dataset_name))
        output_dir_semimu = ("HTT/Extracted/Systematics/%s/semimu/%s" % (d, dataset_name))

        full_template = copy.copy(multicrab)
        if "EMEnriched" in dataset_path or "BCtoE" in dataset_path:
            full_template.template += multicrab_semie
        elif "MuEnriched" in dataset_path:
            full_template.template += multicrab_semimu
        else:
            full_template.template += multicrab_semie
            full_template.template += multicrab_semimu

        print("Creating config file for '%s'" % (dataset_path))
        print("\tName: %s" % dataset_name)
        print("\tOutput directory (semi-mu): %s" % output_dir_semimu)
        print("\tOutput directory (semi-e): %s" % output_dir_semie)
        print("")

        f = open(output_file, "w")
        f.write(full_template.substitute(ui_working_dir=ui_working_dir, dataset=dataset_path, remote_dir_semie=output_dir_semie, remote_dir_semimu=output_dir_semimu, name=dataset_name, email=email, events=dataset_size))
        f.close()

    if options.run:
      cmd = "./multicrab.sh -create -submit -cfg %s" % (output_file)
      os.system(cmd)

    if options.status:
      cmd = "./multicrab.sh -status -c %s" % (ui_working_dir)
      os.system(cmd)

    if options.get:
      cmd = "./multicrab.sh -get -c %s" % (ui_working_dir)
      os.system(cmd)

    if options.resubmit:
      cmd = "./multicrab.sh -resubmit bad -c %s" % (ui_working_dir)
      os.system(cmd)

    if options.submit:
      cmd = "./multicrab.sh -submit all -c %s" % (ui_working_dir)
      os.system(cmd)

import multiprocessing
pool = multiprocessing.Pool(options.cores)
pool.map(processDataset, datasets)
