# D->KKpipi from B->DK decays binned fit

This repository contains ```C++```  code for doing a binned fit of D->KKpipi from B->DK decays. The main goal is to extract gamma.

Tested with the follwing setup:
* cmake version 3
* C++ compiler with C++14
* ROOT version 6.22

To install, run
 ```shell
git clone git@github.com:MartinDuyTat/KKpipi_Binned_Fit.git
cd KKpipi_Binned_Fit
mkdir build
cmake ..
make
make install