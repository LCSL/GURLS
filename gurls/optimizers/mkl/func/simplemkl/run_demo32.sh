#!/bin/bash

export LD_LIBRARY_PATH=$PWD/ipopt/ipopt-3.3.5-mumps-i686/lib:$LD_LIBRARY_PATH
matlab -nojvm -r "simplemkl_example; exit;"
matlab -nojvm -r "ikl_example; exit;"
