module load cmake
source /opt/intel/parallel_studio_2019/parallel_studio_xe_2019.5.075/bin/psxevars.sh intel64
ifort -O3 -assume protect_parens,buffered_io,byterecl,minus0 -o DCD_kernel64 DCD_kernel.f
mv DCD_kernel64 DCD_kernel

