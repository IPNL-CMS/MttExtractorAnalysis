#! /usr/bin/env python

import os, copy, datetime, pwd, re

def check_output(*popenargs, **kwargs):
    import subprocess
    r"""Run command with arguments and return its output as a byte string.
 
    Backported from Python 2.7 as it's implemented as pure python on stdlib.
 
    >>> check_output(['/usr/bin/python', '--version'])
    Python 2.6.2
    """
    process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
    output, unused_err = process.communicate()
    retcode = process.poll()
    if retcode:
        cmd = kwargs.get("args")
        if cmd is None:
            cmd = popenargs[0]
        error = subprocess.CalledProcessError(retcode, cmd)
        error.output = output
        raise error
    return output

def getGitTag():
    if "tag" not in getGitTag.__dict__:
        getGitTag.tag = check_output(["git", "describe", "--tags"]).rstrip('\n')

    return getGitTag.tag

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-j", "--process", action="store", dest="cores", type="int", default=1, help="Number of core to use for launching")
parser.add_option("", "--create-cfg", action="store_true", dest="create_cfg", default=False, help="create config files for crab")
parser.add_option("", "--run", action="store_true", dest="run", default=False, help="run crab")
parser.add_option("", "--status", action="store_true", dest="status", default=False, help="run crab -status")
parser.add_option("", "--get", action="store_true", dest="get", default=False, help="run crab -get")
parser.add_option("", "--resubmit", action="store_true", dest="resubmit", default=False, help="run crab -resubmit bad")
parser.add_option("", "--submit", action="store_true", dest="submit", default=False, help="run crab -submit all")
parser.add_option("", "--kill", action="store_true", dest="kill", default=False, help="run crab -kill all")
(options, args) = parser.parse_args()

if options.run:
  options.create_cfg = True

datasets = [

    # Single anti-top
    ["/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola/chassera-Tbar_t-channel_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "Tbar_t-channel"],
    ["/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/chassera-Tbar_tW-channel_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "Tbar_tW-channel"],
    ["/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola/chassera-Tbar_s-channel_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "Tbar_s-channel"],
    
    # Single top
    ["/T_t-channel_TuneZ2star_8TeV-powheg-tauola/chassera-T_t-channel_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "T_t-channel"],
    ["/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/chassera-T_tW-channel_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "T_tW-channel"],
    ["/T_s-channel_TuneZ2star_8TeV-powheg-tauola/chassera-T_s-channel_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "T_s-channel"],

    # TT + jets
    ["/TT_CT10_TuneZ2star_8TeV-powheg-tauola/chassera-TT_CT10_powheg_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "TT_powheg"],
    ["/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/chassera-TTJets_MassiveBinDECAY_START53_V7C_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "TT_madgraph"],
    #["/TT_8TeV-mcatnlo/chassera-TT_8TeV-mcatnlo_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "TT_mcatnlo"],
    
    # Z + jets
    ["/DY1JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/chassera-DY1JetsToLL_M-50_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "DY1JetsToLL_M-50"],
    ["/DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/chassera-DY2JetsToLL_M-50_START53_V7C_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "DY2JetsToLL_M-50"],
    ["/DY3JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/chassera-DY3JetsToLL_M-50_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "DY3JetsToLL_M-50"],
    ["/DY4JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/chassera-DY4JetsToLL_M-50_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "DY4JetsToLL_M-50"],

    # W + jets
    ["/W1JetsToLNu_TuneZ2Star_8TeV-madgraph/chassera-W1JetsToLNu_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "W1JetsToLNu"],
    ["/W2JetsToLNu_TuneZ2Star_8TeV-madgraph/chassera-W2JetsToLNu_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "W2JetsToLNu"],
    ["/W3JetsToLNu_TuneZ2Star_8TeV-madgraph/chassera-W3JetsToLNu_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "W3JetsToLNu"],
    ["/W4JetsToLNu_TuneZ2Star_8TeV-madgraph/chassera-W4JetsToLNu_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "W4JetsToLNu"],

    # QCD
    #["/QCD_Pt_20_30_EMEnriched_TuneZ2star_8TeV_pythia6/sbrochet-QCD_Pt_20_30_EMEnriched_START53_V7A_08Dec13-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "QCD_Pt_20_30_EMEnriched"],
    #["/QCD_Pt_30_80_EMEnriched_TuneZ2star_8TeV_pythia6/sbrochet-QCD_Pt_30_80_EMEnriched_START53_V7A_08Dec13-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "QCD_Pt_30_80_EMEnriched"],
    #["/QCD_Pt_80_170_EMEnriched_TuneZ2star_8TeV_pythia6/sbrochet-QCD_Pt_80_170_EMEnriched_START53_V7A_08Dec13-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "QCD_Pt_80_170_EMEnriched"],
    #["/QCD_Pt_170_250_EMEnriched_TuneZ2star_8TeV_pythia6/sbrochet-QCD_Pt_170_250_EMEnriched_START53_V7A_08Dec13-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "QCD_Pt_170_250_EMEnriched"],
    #["/QCD_Pt_250_350_EMEnriched_TuneZ2star_8TeV_pythia6/sbrochet-QCD_Pt_250_350_EMEnriched_START53_V7A_08Dec13-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "QCD_Pt_250_350_EMEnriched"],
    #["/QCD_Pt_350_EMEnriched_TuneZ2star_8TeV_pythia6/sbrochet-QCD_Pt_350_EMEnriched_START53_V7A_08Dec13-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "QCD_Pt_350_EMEnriched"],

    #["/QCD_Pt_30_80_BCtoE_TuneZ2star_8TeV_pythia6/sbrochet-QCD_Pt_30_80_BCtoE_16Dec13-v1-000dcc01450dd869c68c0ab31d140828/USER", "QCD_Pt_30_80_BCtoE"],
    #["/QCD_Pt_80_170_BCtoE_TuneZ2star_8TeV_pythia6/sbrochet-QCD_Pt_80_170_BCtoE_16Dec13-v1-000dcc01450dd869c68c0ab31d140828/USER", "QCD_Pt_80_170_BCtoE"],
    #["/QCD_Pt_170_250_BCtoE_TuneZ2star_8TeV_pythia6/sbrochet-QCD_Pt_170_250_BCtoE_16Dec13-v1-000dcc01450dd869c68c0ab31d140828/USER", "QCD_Pt_170_250_BCtoE"],
    #["/QCD_Pt_250_350_BCtoE_TuneZ2star_8TeV_pythia6/sbrochet-QCD_Pt_250_350_BCtoE_16Dec13-v1-000dcc01450dd869c68c0ab31d140828/USER", "QCD_Pt_250_350_BCtoE"],
    #["/QCD_Pt_350_BCtoE_TuneZ2star_8TeV_pythia6/sbrochet-QCD_Pt_350_BCtoE_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "QCD_Pt_350_BCtoE"],

    #["/QCD_pt15to30_bEnriched_MuEnrichedPt14_TuneZ2star_8TeV-pythia6/sbrochet-QCD_pt15to30_bEnriched_MuEnrichedPt14_START53_V7A_08Dec13-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "QCD_pt15to30_bEnriched_MuEnrichedPt14"],
    #["/QCD_pt30to50_bEnriched_MuEnrichedPt14_TuneZ2star_8TeV-pythia6/sbrochet-QCD_pt30to50_bEnriched_MuEnrichedPt14_START53_V7A_08Dec13-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "QCD_pt30to50_bEnriched_MuEnrichedPt14"],
    #["/QCD_pt50to150_bEnriched_MuEnrichedPt14_TuneZ2star_8TeV-pythia6/sbrochet-QCD_pt50to150_bEnriched_MuEnrichedPt14_START53_V7A_08Dec13-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "QCD_pt50to150_bEnriched_MuEnrichedPt14"],
    #["/QCD_pt150_bEnriched_MuEnrichedPt14_TuneZ2star_8TeV-pythia6/sbrochet-QCD_pt150_bEnriched_MuEnrichedPt14_START53_V7A_08Dec13-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "QCD_pt150_bEnriched_MuEnrichedPt14"],
    
    # Signal
    ["/S0_S_i_M400_cpl1_pseudoscalar_15Dec13_START53_V7C-GEN/sbrochet-S0_S_i_M400_cpl1_pseudoscalar_START53_V7C_09Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "S0_S_i_M400_cpl1_pseudoscalar"],
    ["/S0_S_i_M500_cpl1_pseudoscalar_15Dec13_START53_V7C-GEN/sbrochet-S0_S_i_M500_cpl1_pseudoscalar_START53_V7C_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "S0_S_i_M500_cpl1_pseudoscalar"],
    ["/S0_S_i_M600_cpl1_pseudoscalar_15Dec13_START53_V7C-GEN/sbrochet-S0_S_i_M600_cpl1_pseudoscalar_START53_V7C_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "S0_S_i_M600_cpl1_pseudoscalar"],
    ["/S0_S_i_M700_cpl1_pseudoscalar_15Dec13_START53_V7C-GEN/sbrochet-S0_S_i_M700_cpl1_pseudoscalar_START53_V7C_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "S0_S_i_M700_cpl1_pseudoscalar"],
    ["/S0_S_i_M800_cpl1_pseudoscalar_15Dec13_START53_V7C-GEN/sbrochet-S0_S_i_M800_cpl1_pseudoscalar_START53_V7C_12Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "S0_S_i_M800_cpl1_pseudoscalar"],

    ["/S0_S_i_M400_cpl1_scalar_15Dec13_START53_V7C-GEN/sbrochet-S0_S_i_M400_cpl1_scalar_START53_V7C_18Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "S0_S_i_M400_cpl1_scalar"],
    ["/S0_S_i_M500_cpl1_scalar_15Dec13_START53_V7C-GEN/sbrochet-S0_S_i_M500_cpl1_scalar_START53_V7C_12Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "S0_S_i_M500_cpl1_scalar"],
    ["/S0_S_i_M600_cpl1_scalar_15Dec13_START53_V7C-GEN/sbrochet-S0_S_i_M600_cpl1_scalar_START53_V7C_18Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "S0_S_i_M600_cpl1_scalar"],
    ["/S0_S_i_M700_cpl1_scalar_15Dec13_START53_V7C-GEN/sbrochet-S0_S_i_M700_cpl1_scalar_START53_V7C_18Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "S0_S_i_M700_cpl1_scalar"],
    ["/S0_S_i_M800_cpl1_scalar_15Dec13_START53_V7C-GEN/sbrochet-S0_S_i_M800_cpl1_scalar_START53_V7C_18Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "S0_S_i_M800_cpl1_scalar"],
    
    # Z' narrow
    ["/ZPrimeToTTJets_M500GeV_W5GeV_TuneZ2star_8TeV-madgraph-tauola/apequegn-Zprime_500_Narrow_START53_V7A_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "ZPrimeToTTJets_M500GeV_W5GeV"],
    ["/ZPrimeToTTJets_M750GeV_W7p5GeV_TuneZ2star_8TeV-madgraph-tauola/apequegn-Zprime_750_Narrow_START53_V7A_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "ZPrimeToTTJets_M750GeV_W7p5GeV"],
    ["/ZPrimeToTTJets_M1000GeV_W10GeV_TuneZ2star_8TeV-madgraph-tauola/apequegn-Zprime_1000_Narrow_START53_V7A_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "ZPrimeToTTJets_M1000GeV_W10GeV"],
    ["/ZPrimeToTTJets_M1250GeV_W12p5GeV_TuneZ2star_8TeV-madgraph-tauola/apequegn-Zprime_1250_Narrow_START53_V7A_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "ZPrimeToTTJets_M1250GeV_W12p5GeV"],
    ["/ZPrimeToTTJets_M1500GeV_W15GeV_TuneZ2star_8TeV-madgraph-tauola/apequegn-Zprime_1500_Narrow_START53_V7A_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "ZPrimeToTTJets_M1500GeV_W15GeV"],
    ["/ZPrimeToTTJets_M2000GeV_W20GeV_TuneZ2star_8TeV-madgraph-tauola/apequegn-Zprime_2000_Narrow_START53_V7A_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "ZPrimeToTTJets_M2000GeV_W20GeV"],
    
    # Z' narrow ext
    ["/ZPrimeToTTJets_M500GeV_W5GeV_TuneZ2star_8TeV_ext-madgraph-tauola/apequegn-Zprime_500_Narrow_ext_START53_V7C_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "ZPrimeToTTJets_M500GeV_W5GeV_ext"],
    ["/ZPrimeToTTJets_M750GeV_W7p5GeV_TuneZ2star_8TeV_ext-madgraph-tauola/apequegn-Zprime_750_Narrow_ext_START53_V7C_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "ZPrimeToTTJets_M750GeV_W7p5GeV_ext"],
    ["/ZPrimeToTTJets_M1000GeV_W10GeV_TuneZ2star_8TeV_ext-madgraph-tauola/apequegn-Zprime_1000_Narrow_ext_START53_V7C_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "ZPrimeToTTJets_M1000GeV_W10GeV_ext"],
    ["/ZPrimeToTTJets_M1250GeV_W12p5GeV_TuneZ2star_8TeV_ext-madgraph-tauola/apequegn-Zprime_1250_Narrow_ext_START53_V7C_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "ZPrimeToTTJets_M1250GeV_W12p5GeV_ext"],
    ["/ZPrimeToTTJets_M1500GeV_W15GeV_TuneZ2star_8TeV_ext-madgraph-tauola/apequegn-Zprime_1500_Narrow_ext_START53_V7C_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "ZPrimeToTTJets_M1500GeV_W15GeV_ext"],
    ["/ZPrimeToTTJets_M2000GeV_W20GeV_TuneZ2star_8TeV_ext-madgraph-tauola/apequegn-Zprime_2000_Narrow_ext_START53_V7C_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "ZPrimeToTTJets_M2000GeV_W20GeV_ext"],
    
    # Z' large
    ["/ZPrimeToTTJets_M500GeV_W50GeV_TuneZ2star_8TeV-madgraph-tauola/apequegn-Zprime_500_Large_START53_V7A_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "ZPrimeToTTJets_M500GeV_W50GeV"],
    ["/ZPrimeToTTJets_M750GeV_W75GeV_TuneZ2star_8TeV-madgraph-tauola/apequegn-Zprime_750_Large_START53_V7A_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "ZPrimeToTTJets_M750GeV_W75GeV"],
    ["/ZPrimeToTTJets_M1000GeV_W100GeV_TuneZ2star_8TeV-madgraph-tauola/apequegn-Zprime_1000_Large_START53_V7A_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "ZPrimeToTTJets_M1000GeV_W100GeV"],
    ["/ZPrimeToTTJets_M1250GeV_W125GeV_TuneZ2star_8TeV-madgraph-tauola/apequegn-Zprime_1250_Large_START53_V7A_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "ZPrimeToTTJets_M1250GeV_W125GeV"],
    ["/ZPrimeToTTJets_M1500GeV_W150GeV_TuneZ2star_8TeV-madgraph-tauola/apequegn-Zprime_1500_Large_START53_V7A_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "ZPrimeToTTJets_M1500GeV_W150GeV"], 
    ["/ZPrimeToTTJets_M2000GeV_W200GeV_TuneZ2star_8TeV-madgraph-tauola/apequegn-Zprime_2000_Large_START53_V7A_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "ZPrimeToTTJets_M2000GeV_W200GeV"], 
    
    # Z' large ext
    ["/ZPrimeToTTJets_M500GeV_W50GeV_TuneZ2star_8TeV_ext-madgraph-tauola/apequegn-Zprime_500_Large_ext_START53_V7C_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "ZPrimeToTTJets_M500GeV_W50GeV_ext"],
    ["/ZPrimeToTTJets_M750GeV_W75GeV_TuneZ2star_8TeV_ext-madgraph-tauola/apequegn-Zprime_750_Large_ext_START53_V7C_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "ZPrimeToTTJets_M750GeV_W75GeV_ext"],
    ["/ZPrimeToTTJets_M1000GeV_W100GeV_TuneZ2star_8TeV_ext-madgraph-tauola/apequegn-Zprime_1000_Large_ext_START53_V7C_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "ZPrimeToTTJets_M1000GeV_W100GeV_ext"], 
    ["/ZPrimeToTTJets_M1250GeV_W125GeV_TuneZ2star_8TeV_ext-madgraph-tauola/apequegn-Zprime_1250_Large_ext_START53_V7C_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "ZPrimeToTTJets_M1250GeV_W125GeV_ext"],
    ["/ZPrimeToTTJets_M1500GeV_W150GeV_TuneZ2star_8TeV_ext-madgraph-tauola/apequegn-Zprime_1500_Large_ext_START53_V7C_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "ZPrimeToTTJets_M1500GeV_W150GeV_ext"],
    ["/ZPrimeToTTJets_M2000GeV_W200GeV_TuneZ2star_8TeV_ext-madgraph-tauola/apequegn-Zprime_2000_Large_ext_START53_V7C_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "ZPrimeToTTJets_M2000GeV_W200GeV_ext"], 
    
    # RS gluons
    ["/RSGluonToTT_M-700_Tune4C_8TeV-pythia8/apequegn-RSGluon700_START53_V7A_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "RSGluonToTT_M-700"],
    ["/RSGluonToTT_M-1000_Tune4C_8TeV-pythia8/apequegn-RSGluon1000_START53_V7A_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "RSGluonToTT_M-1000"],
    ["/RSGluonToTT_M-1200_Tune4C_8TeV-pythia8/apequegn-RSGluon1200_START53_V7A_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "RSGluonToTT_M-1200"],
    ["/RSGluonToTT_M-1500_Tune4C_8TeV-pythia8/apequegn-RSGluon1500_START53_V7A_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "RSGluonToTT_M-1500"],
    ["/RSGluonToTT_M-2000_Tune4C_8TeV-pythia8/apequegn-RSGluon2000_START53_V7A_08Jan14-v1-2051feb9cae877bd6999b7eb77fd9c68/USER", "RSGluonToTT_M-2000"],

    # Dibosons
    ["/WW_TuneZ2star_8TeV_pythia6_tauola/chassera-WWincl_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "WW"],
    ["/WZ_TuneZ2star_8TeV_pythia6_tauola/chassera-WZincl_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER", "WZ"],
    ["/ZZ_TuneZ2star_8TeV_pythia6_tauola/chassera-ZZincl_START53_V7A_03May13-v1-01f389e36b58797d8560cb86e692fc11/USER ", "ZZ"]
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
    output_file = "multicrab_MC_%s_%s_extractor_%s.cfg" % (dataset_name, d, getGitTag())
    ui_working_dir = ("multicrab_MC_%s_extractor") % (dataset_name)

    if options.create_cfg:
        output_dir_semie = ("HTT/Extracted/MC/Summer12/extractor_%s/%s/semie/%s" % (getGitTag(), d, dataset_name))
        output_dir_semimu = ("HTT/Extracted/MC/Summer12/extractor_%s/%s/semimu/%s" % (getGitTag(), d, dataset_name))

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

    if options.kill:
        cmd = "./multicrab.sh -kill all -c %s" % (ui_working_dir)
        os.system(cmd)

import multiprocessing
pool = multiprocessing.Pool(options.cores)
pool.map(processDataset, datasets)
