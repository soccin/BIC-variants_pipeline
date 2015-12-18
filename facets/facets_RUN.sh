#! /bin/bash

FACETS=$1
DIR=$2
TAG=$3
FILE=$4
PC=300
C=100

x=3;

while [ $x != 0 ]
do
    echo "X: $x"
    echo "C: $C"

    $FACETS/facets doFacets -D $DIR -t $TAG -f $FILE -G T -pc $PC -c $C
     
    EXITVAL=$?;

    #if [ $EXITVAL == 0 ]
    if [ -e $DIR/$TAG_hisens.CNCF.png ]
    then
        x=1;
    fi

    echo "EXITVALUE=$EXITVAL";
    C=$[$C+$C]
    x=$[$x-1]
done

