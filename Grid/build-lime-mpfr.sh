#! /bin/bash

# Builds two packages needed by Grid

##################
# prerequisites: Lime, MPFR place in $HOME/prefix
##################
export prefix=$HOME/prefix
export grid=$HOME/GridCompile
mkdir $prefix
##################
#LIME
##################
cd $prefix
wget http://usqcd-software.github.io/downloads/c-lime/lime-1.3.2.tar.gz
tar xvzf lime-1.3.2.tar.gz
cd lime-1.3.2
./configure --prefix $prefix
make all install
##################
#MPFR - summit is badly configured
##################
cd $prefix
wget https://www.mpfr.org/mpfr-current/mpfr-4.0.2.tar.gz
tar xvzf mpfr-4.0.2.tar.gz
cd mpfr-4.0.2
./configure --prefix $prefix
make all install
 
