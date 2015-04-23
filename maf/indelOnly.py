#!/usr/bin/env python2.7
import sys
header=sys.stdin.readline()
print header,
for line in sys.stdin:
    F=line.strip().split("\t")
    if F[9] in ["DEL","INS"]:
        print line,