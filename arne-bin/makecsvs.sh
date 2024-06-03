#!/bin/bash -x



for i in $1 
do 
    for j in ptm pdockq pdockq_v21 pconsdock pdock_fd PconsFoldSeek
    do 
        n=`basename $i`
        cat $i/*_*/$j.csv | grep -v sbatch  > herpes-ppi/data/$n-$j.csv 
    done 
done

# Made a mistake when running 


#for j in multimer* ; do for i in `grep   \* herpes-ppi/data/${j}-ptm.csv | gawk -F "," '{print $2}' `; do rm $j/$i/*.csv ; done ; done
# sleep 7200 ; for j in multimer* ; do bin/makecsvs.sh $j ; for i in `grep   \* herpes-ppi/data/${j}-ptm.csv | gawk -F "," '{print $2}' `; do rm $j/$i/*.csv ; done ; done 