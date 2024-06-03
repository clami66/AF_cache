#!/bin/bash -x

for i in `find multimer/ -name "ranked_0.pdb"`
do
    j=`dirname $i`
    k=`basename $j`
    cp $i $j/$k.pdb
done
