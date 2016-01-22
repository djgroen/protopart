#!/bin/bash
GMYFILE=${1}
python protopart.py  -i ${GMYFILE} -o ${GMYFILE%gmy}out -c ${GMYFILE%gmy}asc --all 4 > protopart.out 
