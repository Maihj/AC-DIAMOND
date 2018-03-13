#!/bin/bash
DB_FAST_INDEX=acd.fast.index
DB_SEN_INDEX=acd.sensitive.index
QUERY=query.fa
UNALIGNED_QUERY=unaligned_queries.fa
OUTPUT1=acd.fast.result
OUTPUT2=acd.extra.result
ACD_FLAGS="-p 24 -z 6 -e 0.001"
./ac-diamond align -d ${DB_FAST_INDEX} -q ${QUERY} -a ${OUTPUT1} ${ACD_FLAGS}
./ac-diamond view -a ${OUTPUT1}.daa -o ${OUTPUT1}.m8
python unaligned.py ${OUTPUT1}.m8 ${QUERY} > ${UNALIGNED_QUERY}
./ac-diamond align -d ${DB_SEN_INDEX} -q ${UNALIGNED_QUERY} -a ${OUTPUT2} ${ACD_FLAGS}
./ac-diamond view -a ${OUTPUT2}.daa -o ${OUTPUT2}.m8
cat ${OUTPUT1}.m8 ${OUTPUT2}.m8