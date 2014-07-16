#! /bin/bash

# JEC
mkdir JECup
cp Extractor_MTT_*.py JECup/
cp createAndRunMCCrab.py JECup/
cp crab_MC.cfg.template.ipnl JECup/
cp createOutputListForMC.py JECup/
cp multicrab.sh JECup/

sed -i 's/jesSign = 0/jesSign = 1/' JECup/Extractor_MTT_common.py
sed -i 's/dataset_name = dataset\[1\]/dataset_name = dataset\[1\] + "_JECup"/' JECup/createAndRunMCCrab.py
sed -i 's/Extracted\/MC/Extracted\/Systematics/' JECup/createAndRunMCCrab.py

#sed -i '54,69s/^/#/' JECup/createAndRunMCCrab.py

mkdir JECdown
cp Extractor_MTT_*.py JECdown/
cp createAndRunMCCrab.py JECdown/
cp crab_MC.cfg.template.ipnl JECdown/
cp createOutputListForMC.py JECdown/
cp multicrab.sh JECdown/

sed -i 's/jesSign = 0/jesSign = -1/' JECdown/Extractor_MTT_common.py
sed -i 's/dataset_name = dataset\[1\]/dataset_name = dataset\[1\] + "_JECdown"/' JECdown/createAndRunMCCrab.py
sed -i 's/Extracted\/MC/Extracted\/Systematics/' JECdown/createAndRunMCCrab.py

#sed -i '54,69s/^/#/' JECdown/createAndRunMCCrab.py

# JER
mkdir JERup
cp Extractor_MTT_*.py JERup/
cp createAndRunMCCrab.py JERup/
cp crab_MC.cfg.template.ipnl JERup/
cp createOutputListForMC.py JERup/
cp multicrab.sh JERup/

sed -i 's/jerSign = 0/jerSign = 1/' JERup/Extractor_MTT_common.py
sed -i 's/dataset_name = dataset\[1\]/dataset_name = dataset\[1\] + "_JERup"/' JERup/createAndRunMCCrab.py
sed -i 's/Extracted\/MC/Extracted\/Systematics/' JERup/createAndRunMCCrab.py

#sed -i '54,69s/^/#/' JERup/createAndRunMCCrab.py

mkdir JERdown
cp Extractor_MTT_*.py JERdown/
cp createAndRunMCCrab.py JERdown/
cp crab_MC.cfg.template.ipnl JERdown/
cp createOutputListForMC.py JERdown/
cp multicrab.sh JERdown/

sed -i 's/jerSign = 0/jerSign = -1/' JERdown/Extractor_MTT_common.py
sed -i 's/dataset_name = dataset\[1\]/dataset_name = dataset\[1\] + "_JERdown"/' JERdown/createAndRunMCCrab.py
sed -i 's/Extracted\/MC/Extracted\/Systematics/' JERdown/createAndRunMCCrab.py

#sed -i '54,69s/^/#/' JERdown/createAndRunMCCrab.py
