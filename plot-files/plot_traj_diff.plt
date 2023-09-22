# fname1 = 'results/HMS-temp/20-std/fullsol.dat'
# fname2 = 'results/HMS-temp/20-std-nav-5/fullsol.dat'
filename = 'results/HMS-temp/fullsol_paste.dat'

set term pngcairo enhanced color font "Arial-Bold, 20"  size 1600, 1600 #lw 3
set output "results/HMS-temp/traj_diff.png"
set size ratio -1 # axis equal
set xrange [-1.8:1.8]
set yrange [-1.8:1.8]


SCALE = 100.

rf = 1.6507
set parametric
set trange [0:2*pi]
# set style ellipse units xx

plot cos(t), sin(t) notitle lw 2, \
     rf*cos(t), rf*sin(t) notitle lw 2 , \
     filename using 2:3 with l lc "red" lw 2 title "20-std", \
     filename using ($2)+SCALE*($9-$2):($3)+SCALE*($10-$3) with l lc "blue" lw 2 title "20-std"
     # fname1 using 2:3 with l lc "red" lw 2 title "20-std", \
     # fname2 using 2:3 with l lc "blue" lw 2 title "20-std-nav"