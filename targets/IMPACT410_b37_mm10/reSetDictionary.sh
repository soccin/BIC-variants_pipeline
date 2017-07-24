#!/bin/bash
NEWDICT=$1
ILIST=$2

cat $NEWDICT
egrep -v "^@" $ILIST

