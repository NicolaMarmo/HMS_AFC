#!/bin/bash

#indir="results/HMS-temp/"
indir="dbSOL/D+ol/" 
outdir="/" 
gnuplot -e "indir=\"${indir}\"; outdir=\"${outdir}\"" "plot-files/plot_traj_RR_HMS_tex_D+.plt"

latexmk -pdf
latexmk -c
rm *.eps *.tex *-inc-eps-converted-to.pdf
mv *.pdf ${indir}