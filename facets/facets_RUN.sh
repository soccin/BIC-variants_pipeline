#! /bin/bash

FACETS_SUITE=$1
FACETS_LIB=$2
DIR=$3
TAG=$4
FILE=$5
GENOME=$6
PC=$7
C=$8

x=3;

while [ $x != 0 ]
do
    echo "X: $x"
    echo "C: $C"

    rm $DIR/*

    $FACETS_SUITE/facets doFacets -D $DIR -t $TAG -f $FILE -g $GENOME -G T -pc $PC -c $C -r $FACETS_LIB
     
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
