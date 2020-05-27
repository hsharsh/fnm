#!/bin/bash
array=(0 0.1 0.2)

for (( i=0; i<3; i++ ))
do
  echo -e "tmax = 500\ndt = 5e-3\nE = 1\nnu = 0\nrho = 1\nalpha = 0 \nsy = 100\nar_tol = 1e-4\nrf = 2\ntc = ${array[i]}\nsrate = 10\nnlyrs = 2\ninit_c = 1\nj_tol = 1e-7" > './parameters.cfg'
  ./a.out > "log$i.log"
  echo "Completed simulation for ${array[i]} time delay."
done

# "tmax = 200\ndt = 2e-2\nE = 1\nnu = 0.3\nrho = 1\nalpha = 1\nsy = 100\nar_tol = 1e-4\nrf = 2\ntc = 0\nsrate = 10\nnlyrs = 2\ninit_c = 1\nj_tol = 1e-8\n"
#
# "tmax = 50\ndt = 1e-2\nE = 1\nnu = 0\nrho = 1\nalpha = ${array[i]} \nsy = 100\nar_tol = 1e-4\nrf = 2\ntc = 0.2\nsrate = 10\nnlyrs = 2\ninit_c = 1\nj_tol = 1e-7"
