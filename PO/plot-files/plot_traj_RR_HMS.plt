indir = "results/HMS-temp/"
outdir = indir
nLeg = 8

fileTraj1 = indir."Traj/Traj1.dat"
fileTraj2 = indir."Traj/Traj2.dat"
fileTraj3 = indir."Traj/Traj3.dat"
fileTraj4 = indir."Traj/Traj4.dat"
fileTraj5 = indir."Traj/Traj5.dat"
fileTraj6 = indir."Traj/Traj6.dat"
fileTraj7 = indir."Traj/Traj7.dat"
fileTraj8 = indir."Traj/Traj8.dat"
file_opt_X = "results/HMS-temp/opt_X.dat"

fileDVs1 = indir."DVs/DVs1.dat"
fileDVs2 = indir."DVs/DVs2.dat"
fileDVs3 = indir."DVs/DVs3.dat"
fileDVs4 = indir."DVs/DVs4.dat"
fileDVs5 = indir."DVs/DVs5.dat"
fileDVs6 = indir."DVs/DVs6.dat"
fileDVs7 = indir."DVs/DVs7.dat"
fileDVs8 = indir."DVs/DVs8.dat"

filePEXY1 = indir."PE/PEXY1.dat"
filePEXY2 = indir."PE/PEXY2.dat"
filePEXY3 = indir."PE/PEXY3.dat"
filePEXY4 = indir."PE/PEXY4.dat"
filePEXY5 = indir."PE/PEXY5.dat"
filePEXY6 = indir."PE/PEXY6.dat"
filePEXY7 = indir."PE/PEXY7.dat"
filePEXY8 = indir."PE/PEXY8.dat"

filePEXZ1 = indir."PE/PEXZ1.dat"
filePEXZ2 = indir."PE/PEXZ2.dat"
filePEXZ3 = indir."PE/PEXZ3.dat"
filePEXZ4 = indir."PE/PEXZ4.dat"
filePEXZ5 = indir."PE/PEXZ5.dat"
filePEXZ6 = indir."PE/PEXZ6.dat"
filePEXZ7 = indir."PE/PEXZ7.dat"
filePEXZ8 = indir."PE/PEXZ8.dat"

filePEYZ1 = indir."PE/PEYZ1.dat"
filePEYZ2 = indir."PE/PEYZ2.dat"
filePEYZ3 = indir."PE/PEYZ3.dat"
filePEYZ4 = indir."PE/PEYZ4.dat"
filePEYZ5 = indir."PE/PEYZ5.dat"
filePEYZ6 = indir."PE/PEYZ6.dat"
filePEYZ7 = indir."PE/PEYZ7.dat"
filePEYZ8 = indir."PE/PEYZ8.dat"

chi_sq = sqrt(5.991)
# semi_major_axis(x) = sqrt(chi_sq * x)

set term pngcairo enhanced color font "Arial-Bold, 20"  size 1600, 1600 #lw 3
set output outdir."trajXY.png"
set size ratio -1 # axis equal
set xrange [0.97:1.04]
set yrange [-0.06:0.06]
set xlabel "x"
set ylabel "y"

rf = 1.6507

scalaP = 1*chi_sq
scalaV = 1e2
e_step = 1

set grid back lt 1 lc "black" lw 0.25
set style line 1 lw 2 lc 1
set style line 2 lw 2 lc 2
set style line 3 lw 2 lc "black"  

set style arrow 1 nohead lw 4
set arrow arrowstyle 1

set parametric
set trange [0:2*pi]
set style ellipse units xx

plot fileTraj1 using 2:3 with l lc "green" lw 2 notitle, \
     fileTraj2 using 2:3 with l lc "yellow" lw 2 notitle, \
     fileTraj3 using 2:3 with l lc "dark-yellow" lw 2 notitle, \
     fileTraj4 using 2:3 with l lc "red" lw 2 notitle, \
     fileTraj5 using 2:3 with l lc "dark-red" lw 2 notitle, \
     fileTraj6 using 2:3 with l lc "purple" lw 2 notitle, \
     fileTraj7 using 2:3 with l lc "blue" lw 2 notitle, \
     fileTraj8 using 2:3 with l lc "dark-green" lw 2 notitle, \
     fileDVs1 using 2:3:($5)*scalaV:($6)*scalaV with vectors filled head lc "green" lw 4 notitle, \
     fileDVs2 using 2:3:($5)*scalaV:($6)*scalaV with vectors filled head lc "yellow" lw 4 notitle, \
     fileDVs3 using 2:3:($5)*scalaV:($6)*scalaV with vectors filled head lc "dark-yellow" lw 4 notitle, \
     fileDVs4 using 2:3:($5)*scalaV:($6)*scalaV with vectors filled head lc "red" lw 4 notitle, \
     fileDVs5 using 2:3:($5)*scalaV:($6)*scalaV with vectors filled head lc "dark-red" lw 4 notitle, \
     fileDVs6 using 2:3:($5)*scalaV:($6)*scalaV with vectors filled head lc "purple" lw 4 notitle, \
     fileDVs7 using 2:3:($5)*scalaV:($6)*scalaV with vectors filled head lc "blue" lw 4 notitle, \
     fileDVs8 using 2:3:($5)*scalaV:($6)*scalaV with vectors filled head lc "dark-green" lw 4 notitle, \
     filePEXY1 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "green" lw 2 notitle,\
     filePEXY2 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "yellow" lw 2 notitle,\
     filePEXY3 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "dark-yellow" lw 2 notitle,\
     filePEXY4 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "red" lw 2 notitle,\
     filePEXY5 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "dark-red" lw 2 notitle,\
     filePEXY6 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "purple" lw 2 notitle,\
     filePEXY7 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "blue" lw 2 notitle,\
     filePEXY8 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "dark-green" lw 2 notitle

# ---------------------------------------------------------------#
#                           XZ                                   #
# ---------------------------------------------------------------# 
set term pngcairo enhanced color font "Arial-Bold, 20"  size 1600, 1600 #lw 3
set output outdir."trajXZ.png"
set size ratio -1 # axis equal
set xrange [0.97:1.04]
set yrange [-0.2:0.05]
set xlabel "x"
set ylabel "z"

plot fileTraj1 using 2:4 with l lc "green" lw 2 notitle, \
     fileTraj2 using 2:4 with l lc "yellow" lw 2 notitle, \
     fileTraj3 using 2:4 with l lc "dark-yellow" lw 2 notitle, \
     fileTraj4 using 2:4 with l lc "red" lw 2 notitle, \
     fileTraj5 using 2:4 with l lc "dark-red" lw 2 notitle, \
     fileTraj6 using 2:4 with l lc "purple" lw 2 notitle, \
     fileTraj7 using 2:4 with l lc "blue" lw 2 notitle, \
     fileTraj8 using 2:4 with l lc "dark-green" lw 2 notitle, \
     fileDVs1 using 2:4:($5)*scalaV:($7)*scalaV with vectors filled head lc "green" lw 4 notitle, \
     fileDVs2 using 2:4:($5)*scalaV:($7)*scalaV with vectors filled head lc "yellow" lw 4 notitle, \
     fileDVs3 using 2:4:($5)*scalaV:($7)*scalaV with vectors filled head lc "dark-yellow" lw 4 notitle, \
     fileDVs4 using 2:4:($5)*scalaV:($7)*scalaV with vectors filled head lc "red" lw 4 notitle, \
     fileDVs5 using 2:4:($5)*scalaV:($7)*scalaV with vectors filled head lc "dark-red" lw 4 notitle, \
     fileDVs6 using 2:4:($5)*scalaV:($7)*scalaV with vectors filled head lc "purple" lw 4 notitle, \
     fileDVs7 using 2:4:($5)*scalaV:($7)*scalaV with vectors filled head lc "blue" lw 4 notitle, \
     fileDVs8 using 2:4:($5)*scalaV:($7)*scalaV with vectors filled head lc "dark-green" lw 4 notitle,\
     filePEXZ1 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "green" lw 2 notitle,\
     filePEXZ2 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "yellow" lw 2 notitle,\
     filePEXZ3 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "dark-yellow" lw 2 notitle,\
     filePEXZ4 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "red" lw 2 notitle,\
     filePEXZ5 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "dark-red" lw 2 notitle,\
     filePEXZ6 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "purple" lw 2 notitle,\
     filePEXZ7 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "blue" lw 2 notitle,\
     filePEXZ8 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "dark-green" lw 2 notitle

# ---------------------------------------------------------------#
#                           YZ                                   #
# ---------------------------------------------------------------# 
set term pngcairo enhanced color font "Arial-Bold, 20"  size 1600, 1600 #lw 3
set output outdir."trajYZ.png"
set xrange [-0.06:0.06]
set yrange [-0.2:0.05]
set xlabel "y"
set ylabel "z"

plot fileTraj1 using 3:4 with l lc "green" lw 2 notitle, \
     fileTraj2 using 3:4 with l lc "yellow" lw 2 notitle, \
     fileTraj3 using 3:4 with l lc "dark-yellow" lw 2 notitle, \
     fileTraj4 using 3:4 with l lc "red" lw 2 notitle, \
     fileTraj5 using 3:4 with l lc "dark-red" lw 2 notitle, \
     fileTraj6 using 3:4 with l lc "purple" lw 2 notitle, \
     fileTraj7 using 3:4 with l lc "blue" lw 2 notitle, \
     fileTraj8 using 3:4 with l lc "dark-green" lw 2 notitle, \
     fileDVs1 using 3:4:($6)*scalaV:($7)*scalaV with vectors filled head lc "green" lw 4 notitle, \
     fileDVs2 using 3:4:($6)*scalaV:($7)*scalaV with vectors filled head lc "yellow" lw 4 notitle, \
     fileDVs3 using 3:4:($6)*scalaV:($7)*scalaV with vectors filled head lc "dark-yellow" lw 4 notitle, \
     fileDVs4 using 3:4:($6)*scalaV:($7)*scalaV with vectors filled head lc "red" lw 4 notitle, \
     fileDVs5 using 3:4:($6)*scalaV:($7)*scalaV with vectors filled head lc "dark-red" lw 4 notitle, \
     fileDVs6 using 3:4:($6)*scalaV:($7)*scalaV with vectors filled head lc "purple" lw 4 notitle, \
     fileDVs7 using 3:4:($6)*scalaV:($7)*scalaV with vectors filled head lc "blue" lw 4 notitle, \
     fileDVs8 using 3:4:($6)*scalaV:($7)*scalaV with vectors filled head lc "dark-green" lw 4 notitle,\
     filePEYZ1 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "green" lw 2 notitle,\
     filePEYZ2 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "yellow" lw 2 notitle,\
     filePEYZ3 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "dark-yellow" lw 2 notitle,\
     filePEYZ4 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "red" lw 2 notitle,\
     filePEYZ5 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "dark-red" lw 2 notitle,\
     filePEYZ6 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "purple" lw 2 notitle,\
     filePEYZ7 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "blue" lw 2 notitle,\
     filePEYZ8 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "dark-green" lw 2 notitle

# ---------------------------------------------------------------#
#                       Gain - plot                              #
# ---------------------------------------------------------------# 
reset
set term pngcairo enhanced color font "Arial-Bold, 20"  size 1600, 800 #lw 3
set output outdir."control.png"
# set size ratio -1 # axis equal

set boxwidth 0.5
set style fill solid

set xrange [-0.1 : 40.1]
plot for [i=10:18] file_opt_X using i notitle with boxes
# plot file_opt_X using 20 