#!/bin/bash

module load gnuplot
NUM_ATOMS=${1:-1}

(echo -e "ecutrho = N * ecutwfc (Ry)\tTotal Energy (Ry)\tTotal Energy (eV)\tPer-Atom ($NUM_ATOMS) Energy (eV)"; grep "!" ecutrho-mult-*/qe.out | awk -v n="$NUM_ATOMS" -F'[-/:=]+' '{ry=-$7; ev=ry*13.6057039763; printf "%d\t%.8f\t%.8f\t%.8f\n", $3, ry, ev, ev/n}' | sort -n) > energy-vs-ecutrhomult.dat

mv *.slurm-out *.err ./slurm-out-err

gnuplot create-plot.gnuplot
