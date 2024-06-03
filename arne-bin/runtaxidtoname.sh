#!/bin/bash -x

dir=$1

#for i in $dir/[A-[A-Za-z0-9][A-Za-z0-9][A-Za-z0-9][A-Za-z0-9][A-Za-z0-9][A-Za-z0-9]/msas/A/un*a3m;  


#j=`dirname $i` ; k=`dirname $j` ; l=`dirname $k` ; m=`basename $l` ; echo |$m ; 

for l in `grep \> $dir | gawk '{print $2}' |  sed "s/Memid://g" | sort -u `  
do 
    echo -n $l ", " 
    bin/taxidtoname.py $l 
done > $dir.taxid2name

