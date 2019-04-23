#!/bin/csh
set dir="~/osoto_volatile_hallb/dhoutput/skim4March"
set outfile="ALU_run.txt"
rm $outfile >& /dev/null
touch $outfile
foreach fn (`ls ${dir}/*.root`)
set run=`echo ${fn:gas/.root//} | sed 's/.*skim4_//g'`
echo $run
echo $run>> $outfile
tail -3 $run/*bin*/*.log >> $outfile
end



