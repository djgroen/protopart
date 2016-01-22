#!/bin/bash
GMY=$1
RANKS=$2
PROTOPART="/home/schmie/src/protopart.derek/Python/protopart.py"
PiPiSTIENATOR="/home/schmie/src/ppstee_and_deps/ppstee-prefix/src/ppstee/\$INSTALL_DIR/bin/ppstee_hemelb_graph_data_io"
HGB2GEOSTRIP5000="/home/schmie/src/protopart.derek/Python/hgb2geoStrip5000.py"
echo "python ${PROTOPART} -i ${GMY} -o ${GMY%gmy}out -c ${GMY%gmy}hga -b ${GMY%gmy}hgb --all ${RANKS}"
echo "${PiPiSTIENATOR} -i ${GMY%gmy}hgb -o ${GMY%gmy}pipi -P ${RANKS}"
python ${HGB2GEOSTRIP5000} ${GMY%gmy}pipi
