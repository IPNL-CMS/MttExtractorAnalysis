#! /bin/bash

# Perform the same tash as the origin multicrab, but don't create temporary files

if [ -z "$PYTHONPATH" ]; then
  export PYTHONPATH=${CRABDBSAPIPYTHON}:${CRABDLSAPIPYTHON}:${CRABPSETPYTHON}:${CRABPYTHON}
else
  export PYTHONPATH=${CRABDBSAPIPYTHON}:${CRABDLSAPIPYTHON}:${CRABPSETPYTHON}:${PYTHONPATH}:${CRABPYTHON}
fi
#echo $PYTHONPATH

export LD_LIBRARY_PATH=${GLITE_LOCATION}/lib:${LD_LIBRARY_PATH}

# to be removed asap
if [ -z "$CMSSW_VERSION" ]; then
  echo ''
  echo 'crab Error: Please run cmsenv before setting the CRAB environment'
  echo '(see also: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideCrabHowTo#Setup_local_Environment) '
  echo ''
else
# Define working directory for temporary files.
  multicrab_ouput=$(python $CRABPYTHON/multicrab.py $*)
  echo $multicrab_ouput
  MULTICRAB_WORKDIR=`echo $multicrab_ouput | grep "working directory" | sed -e 's/\(.*\)working directory//'`
  /usr/bin/env sh $MULTICRAB_WORKDIR/multicrab.exe
fi
