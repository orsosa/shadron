#!/bin/csh
set ddir="~/osoto_volatile_hallb/dhoutput/data/skim4March"
foreach fn (`ls ${ddir}/*[3-4]???.root`)
    set cdir=`pwd`
    echo "fn: "$fn
    set log=${cdir}"/`basename $fn`.log"
    echo "log: "$log
    ./process_file.sh $fn >& $log &
    sleep 1
end

