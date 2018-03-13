#!/usr/bin/python

'''
Usage: python cmp-query.py BLASTX_ALIGNMENT_RESULT OTHER_TOOL_ALIGNMENT_RESULT
'''

import sys

if len(sys.argv) != 3:
    print "Usage: python cmp-query.py BLASTX_ALIGNMENT_RESULT OTHER_TOOL_ALIGNMENT_RESUL"
    sys.exit()

overlap = 0
dic = {}
blastx_file_name = sys.argv[1]
otherTool_file_name = sys.argv[2]

file1 = file(blastx_file_name, 'r')
file2 = file(otherTool_file_name, 'r')

for line in file1.readlines():
    aln = line.split()
    key = str(aln[0])
    if key not in dic:
        dic[key] = 1

query = {}
for line in file2.readlines():
    aln1 = line.split()
    key = str(aln1[0])
    if key in dic:
        if key not in query:
            query[key] = 1

print "Sensitivity: " + str(float(len(query)) / float(len(dic)) * 100) + "%"
