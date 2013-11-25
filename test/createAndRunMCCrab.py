#! /usr/bin/env python

import os, copy, datetime, pwd, re

from optparse import OptionParser
parser = OptionParser()
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

    # Single anti-top
    #["/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola/chassera-Tbar_t-channel_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "Tbar_t-channel"],
    #["/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/chassera-Tbar_tW-channel_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "Tbar_tW-channel"],
    #["/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola/chassera-Tbar_s-channel_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "Tbar_s-channel"],
    
    ## Single top
    #["/T_t-channel_TuneZ2star_8TeV-powheg-tauola/chassera-T_t-channel_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "T_t-channel"],
    #["/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/chassera-T_tW-channel_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "T_tW-channel"],
    #["/T_s-channel_TuneZ2star_8TeV-powheg-tauola/chassera-T_s-channel_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "T_s-channel"],

    ## TT + jets
    #["/TT_CT10_TuneZ2star_8TeV-powheg-tauola/chassera-TT_CT10_powheg_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "TT_powheg"],
    
    ## Z + jets
    #["/DY1JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/chassera-DY1JetsToLL_M-50_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "DY1JetsToLL_M-50"],
    #["/DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/chassera-DY2JetsToLL_M-50_START53_V7C_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "DY2JetsToLL_M-50"],
    #["/DY3JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/chassera-DY3JetsToLL_M-50_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "DY3JetsToLL_M-50"],
    #["/DY4JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/chassera-DY4JetsToLL_M-50_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "DY4JetsToLL_M-50"],

    ## W + jets
    #["/W1JetsToLNu_TuneZ2Star_8TeV-madgraph/chassera-W1JetsToLNu_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "W1JetsToLNu"],
    #["/W2JetsToLNu_TuneZ2Star_8TeV-madgraph/chassera-W2JetsToLNu_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "W2JetsToLNu"],
    ["/W3JetsToLNu_TuneZ2Star_8TeV-madgraph/chassera-W3JetsToLNu_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "W3JetsToLNu"],
    #["/W4JetsToLNu_TuneZ2Star_8TeV-madgraph/chassera-W4JetsToLNu_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "W4JetsToLNu"],
    
    # Signal
    #["/S0_S_i_M500_cpl1_scalar_hadronized_fastsim_16Nov13-v1/sbrochet-S0_S_i_M500_cpl1_scalar18Nov13-v1-000dcc01450dd869c68c0ab31d140828/USER", "S0_S_i_M500_cpl1_scalar"],
    #["/S0_S_i_M700_cpl1_scalar_hadronized_fastsim_16Nov13-v1/sbrochet-S0_S_i_M700_cpl1_scalar18Nov13-v1-000dcc01450dd869c68c0ab31d140828/USER	", "S0_S_i_M700_cpl1_scalar"],

    #["/S0_S_i_M500_cpl1_pseudoscalar_hadronized_fastsim_16Nov13-v1/sbrochet-S0_S_i_M500_cpl1_pseudoscalar18Nov13-v1-000dcc01450dd869c68c0ab31d140828/USER", "S0_S_i_M500_cpl1_pseudoscalar"],
    #["/S0_S_i_M700_cpl1_pseudoscalar_hadronized_fastsim_16Nov13-v1/sbrochet-S0_S_i_M700_cpl1_pseudoscalar18Nov13-v1-000dcc01450dd869c68c0ab31d140828/USER", "S0_S_i_M700_cpl1_pseudoscalar"],

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

for dataset in datasets:

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
    output_dir_semie = ("HTT/Extracted/MC/Summer12/%s/semie/%s" % (d, dataset_name))
    output_dir_semimu = ("HTT/Extracted/MC/Summer12/%s/semimu/%s" % (d, dataset_name))

    full_template = copy.copy(multicrab)
    if "EMEnriched" in dataset_path:
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
    cmd = "multicrab -create -submit -cfg %s" % (output_file)
    os.system(cmd)

  if options.status:
    cmd = "multicrab -status -c %s" % (ui_working_dir)
    os.system(cmd)

  if options.get:
    cmd = "multicrab -get -c %s" % (ui_working_dir)
    os.system(cmd)

  if options.resubmit:
    cmd = "multicrab -resubmit bad -c %s" % (ui_working_dir)
    os.system(cmd)

  if options.submit:
    cmd = "multicrab -submit all -c %s" % (ui_working_dir)
    os.system(cmd) 
