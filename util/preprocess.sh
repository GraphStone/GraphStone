#! /bin/bash

fpath=$1
nv_per_par=$2
num_vertices=$3
shuffle_size=$4

indegs-preprocess $fpath $nv_per_par $num_vertices $shuffle_size
preprocess $fpath $nv_per_par $num_vertices $shuffle_size
