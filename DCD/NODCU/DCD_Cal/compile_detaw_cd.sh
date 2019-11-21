#!/bin/bash
source /opt/intel/parallel_studio_2019/compilers_and_libraries/linux/bin/compilervars.sh ia32
ifort -ipo -O3 -no-prec-div -o DETAW_CD DETAW_CD.f
