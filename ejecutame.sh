#!/bin/bash

make
./proyecto 0
gnuplot graphs.plot
./proyecto 10
gnuplot graphs.plot
./proyecto 100
gnuplot graphs.plot
./proyecto 1000
gnuplot graphs.plot
./proyecto 10000
gnuplot graphs.plot
