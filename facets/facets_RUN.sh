#! /bin/bash

FACETS=$1
DIR=$2
TAG=$3
FILE=$4
GENOME=$5
PC=300
C=100

x=3;

while [ $x != 0 ]
do
    echo "X: $x"
    echo "C: $C"

    rm $DIR/*

    $FACETS/facets doFacets -D $DIR -t $TAG -f $FILE -g $GENOME -G T -pc $PC -c $C
     
    if [ -e $DIR"/"$TAG"_hisens.CNCF.png" ]
    then
        x=1;
        exit 0
    fi

    echo $DIR"/"$TAG"_hisens.CNCF.png";

    C=$[$C+$C]
    x=$[$x-1]
done

if [ ! -e $DIR"/"$TAG"_hisens.CNCF.png" ]
then
    exit 1;
fi
