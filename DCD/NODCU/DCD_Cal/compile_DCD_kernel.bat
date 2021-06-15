"c:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\bin\compilervars.bat" ia32 vs2015
rem Improved to use buffered IO - NS
ifort /O3 /assume:protect_parens,buffered_io,bytrerecl,minus0 -o DCD_kernel DCD_kernel.f