#!/bin/csh
set fn=$1
set run=`echo ${fn:gas/.root//} | sed 's/.*skim4_//g'`
echo $run
mkdir $run >& /dev/null
cd $run
pwd
set field=-1
if ($run == 3932) set field=1
set script="../BSA_survey.cxx("\"$fn\"","\"ntuple_data\"",$field)"
echo $script
root -l -b -q "$script"
cd ..

