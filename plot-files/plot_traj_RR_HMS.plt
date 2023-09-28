filename1 = "results/HMS-temp/fullsol1.dat"
filename2 = "results/HMS-temp/fullsol2.dat"
file_opt_X = "results/HMS-temp/opt_X.dat"
fileDVs1 = "results/HMS-temp/DVs1.dat"
fileDVs2 = "results/HMS-temp/DVs2.dat"
fileELLs1 = "results/HMS-temp/covEllipses1.dat"
fileELLs2 = "results/HMS-temp/covEllipses2.dat"
trajE = "results/HMS-temp/trajEfile.dat"
trajM = "results/HMS-temp/trajMfile.dat"

chi_sq = sqrt(5.991)
# semi_major_axis(x) = sqrt(chi_sq * x)

set term pngcairo enhanced color font "Arial-Bold, 20"  size 1600, 1600 #lw 3
set output "results/HMS-temp/traj.png"
set size ratio -1 # axis equal
set xrange [-1.8:1.8]
set yrange [-1.8:1.8]

rf = 1.6507

scalaP = 1*chi_sq
scalaV = 5
e_step = 1

set style line 1 lw 2 lc 1
set style line 2 lw 2 lc 2
set style line 3 lw 2 lc "black"  

     set style arrow 1 nohead lw 4
     set arrow arrowstyle 1

set parametric
set trange [0:2*pi]
set style ellipse units xx

plot trajE using 2:3 with l notitle lw 2 lc "blue", \
     filename1 using 2:3 with l lc "dark-magenta" lw 2 notitle, \
     filename2 using 2:3 with l lc "dark-green" lw 2 notitle, \
     filename2 using 2:3 every ::1::1 with points lc 'dark-green' pt 4 lw 3 ps 4 notitle, \
     fileDVs1 using 2:3:($5)*scalaV:($6)*scalaV with vectors filled head lc "dark-magenta" lw 4 notitle, \
     fileDVs2 using 2:3:($5)*scalaV:($6)*scalaV with vectors filled head lc "dark-green" lw 4 notitle, \
     fileELLs1 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "dark-magenta" lw 2 notitle, \
     fileELLs2 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "dark-green" lw 2 notitle

# ---------------------------------------------------------------#
#                       Gain - plot                              #
# ---------------------------------------------------------------# 
reset
set term pngcairo enhanced color font "Arial-Bold, 20"  size 1600, 800 #lw 3
set output "results/HMS-temp/control.png"
# set size ratio -1 # axis equal

set boxwidth 0.5
set style fill solid

set xrange [-0.1 : 40.1]
plot for [i=10:18] file_opt_X using i notitle with boxes
# plot file_opt_X using 20 