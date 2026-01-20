#!/bin/bash

module load gnuplot
NUM_ATOMS=${1:-1}

(echo -e "k-points N N N 0 0 0\tTotal Energy (Ry)\tTotal Energy (eV)\tPer-Atom ($NUM_ATOMS) Energy (eV)"; grep "!" kpts-*/qe.out | awk -v n="$NUM_ATOMS" -F'[-/:=]+' '{ry=-$6; ev=ry*13.6057039763; printf "%d\t%.8f\t%.8f\t%.8f\n", $2, ry, ev, ev/n}' | sort -n) > energy-vs-kpts.dat

mv *.slurm-out *.err ./slurm-out-err

gnuplot create-plot.gnuplot
