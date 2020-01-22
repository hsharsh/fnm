#!/bin/bash
array=(0.15 0.2 0.25 0.3)

for (( i=0; i<4; i++ ))
do
  echo -e "tmax = 50\ndt = 1e-2\nE = 1\nnu = 0\nrho = 1\nalpha = ${array[i]} \nsy = 100\nar_tol = 1e-4\nrf = 2\ntc = 0.2\nsrate = 10\nnlyrs = 2\ninit_c = 1\nj_tol = 1e-7" > './parameters.cfg'
  g++ -fopenmp -O3 fnm.cpp
  OMP_NUM_THREADS=1 ./a.out > "log$i.log"
  echo "Completed simulation for ${array[i]} damping."
done
