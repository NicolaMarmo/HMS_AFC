filefullsolMC = indir."fullsolMC.dat"
filefullsol = indir."fullsol.dat"
fileDVs = indir."DVs.dat"
fileELLs = indir."covEllipses.dat"
filedPr = indir."dPr.dat"
filedPv = indir."dPv.dat"
fileMPpMC = indir."MPpMC.dat"
fileMPp = indir."MPp.dat"
fileMCMCD = indir."MCMCD.dat"
fileMCND = indir."MCND.dat"
trajE = indir."trajEfile.dat"
trajM = indir."trajMfile.dat"
fileUMC = indir."MUMC.dat"
SVS = "SVS.dat"
SVD = "SVD.dat"
SVSP0 = "SVSP0.dat"
SVDP0 = "SVDP0.dat"

# ---------------------------------------------------------------#
#                             trajMC                             #
# ---------------------------------------------------------------# 

chi_sq = sqrt(5.991)
# semi_major_axis(x) = sqrt(chi_sq * x)

set terminal epslatex standalone size 4, 4 color colortext 10 lw 2  header \
"\\usepackage{amsmath}\n\\usepackage{siunitx}"
set output "trajMC.tex"
set size ratio -1 # axis equal
set xrange [-1.8:1.8]
set yrange [-1.8:1.8]
set xlabel "$x$"
set ylabel "$y$" 

rf = 1.6507
lwd = 2

scalaP = chi_sq
scalaV = 10
e_step = 1

set grid back lt 1 lc rgb "#e6e6e6" lw 0.25
set parametric
set trange [0:2*pi]
set style ellipse units xx

plot 0, 0 with points lc rgb 'orange' pt 7 ps 1.5 notitle, \
     filefullsolMC using 2:3 with l lc "grey" lw 2 notitle, \
     trajE using 2:3 with l notitle lw 2 lc "blue", \
     trajM using 2:3 with l notitle lw 2 lc "brown", \
     filefullsol using 2:3 with l lc "black" lw lwd notitle, \
     fileDVs every :::::0 using 2:3:($5)*scalaV:($6)*scalaV with vectors filled head lc "red" lw lwd notitle, \
     fileELLs using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "black" lw lwd notitle

# ---------------------------------------------------------------#
#                             trajn                              #
# ---------------------------------------------------------------# 
set terminal epslatex standalone size 4, 4 color colortext 10 lw 2  header \
"\\usepackage{amsmath}\n\\usepackage{siunitx}"
set output "trajn.tex"
set size ratio -1 # axis equal
set xrange [-1.8:1.8]
set yrange [-1.8:1.8]
set xlabel "$x$"
set ylabel "$y$" 

rf = 1.6507

scalaP = chi_sq
scalaV = 10
e_step = 1

set parametric
set trange [0:2*pi]
set style ellipse units xx

plot 0, 0 with points lc rgb 'orange' pt 7 ps 1.5 notitle, \
     trajE using 2:3 with l notitle lw 2 lc "blue", \
     trajM using 2:3 with l notitle lw 2 lt 7 lc "brown" , \
     filefullsol using 2:3 with l lc "black" lw lwd notitle,\
     fileDVs every :::::0 using 2:3:($5)*scalaV:($6)*scalaV with vectors filled head lc "red" lw lwd notitle, \
     fileELLs using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "black" lw lwd notitle
     #trajE every ::0::0 using 2:3 with points lc rgb 'blue' pt 4 lw 1.5 ps 1.5 notitle,\
     #trajM every ::0::0 using 2:3 with points lc rgb 'brown' pt 3 lw 1.5 ps 1.5 notitle,\
     #trajM every ::0::0 using 2:3 with points lc rgb 'brown' pt 3 lw 1.5 ps 1.5 notitle,\
     ##'-' with points lc rgb 'brown' pt 3 lw 1.5 ps 1.5 notitle
     ##-1.2 0

# ---------------------------------------------------------------#
#                       dPr, dPv                                 #
# ---------------------------------------------------------------# 
reset
set terminal epslatex standalone size 4, 3 color colortext 10 lw 2 header \
"\\usepackage{amsmath}\n\\usepackage{siunitx}"
set output "dPxy.tex"
set xlabel '$k$'
set ylabel '$\sigma^2$' 

psz = 1.5
lwd1 = 1.5
lwd2 = 1

set grid back lt 1 lc rgb "#e6e6e6" lw 0.25
set key left top inside Left
set key box spacing 1.5 
set key box width -7

set format y "$%2.1t \\times 10^{%L}$"

set xtics 2
plot fileMPp using 1:2 with points lc rgb 'black' pt 4 lw lwd1 ps psz title '$\sigma^2_x$', \
     fileMPpMC using 1:2 with points lc rgb 'black' pt 2 lw lwd1 ps psz title '$\sigma^2_{x, {\mathrm MC}}$',\
     fileMPp using 1:2 with l lc 'black' lw lwd2 dt 6 notitle, \
     fileMPp using 1:3 with points lc rgb 'red' pt 4 lw lwd1 ps psz title '$\sigma^2_y$', \
     fileMPpMC using 1:3 with points lc rgb 'red' pt 2 lw lwd1 ps psz title '$\sigma^2_{y, {\mathrm MC}}$',\
     fileMPp using 1:3 with l lc 'red' lw lwd2 dt 6 notitle

set output "dPvxvy.tex"
plot fileMPp using 1:5 with points lc rgb 'black' pt 4 lw lwd1 ps psz title '$\sigma^2_{v_x}$', \
     fileMPpMC using 1:5 with points lc rgb 'black' pt 2 lw lwd1 ps psz title '$\sigma^2_{{v_x}, {\mathrm MC}}$',\
     fileMPp using 1:5 with l lc 'black' lw lwd2 dt 6 notitle, \
     fileMPp using 1:6 with points lc rgb 'red' pt 4 lw lwd1 ps psz title '$\sigma^2_{v_y}$', \
     fileMPpMC using 1:6 with points lc rgb 'red' pt 2 lw lwd1 ps psz title '$\sigma^2_{{v_y}, {\mathrm MC}}$',\
     fileMPp using 1:6 with l lc 'red' lw lwd2 dt 6 notitle

set key right top inside Left
set output "dPzvz.tex"
plot fileMPp using 1:4 with points lc rgb 'black' pt 4 lw lwd1 ps psz title '$\sigma^2_z$', \
     fileMPpMC using 1:4 with points lc rgb 'black' pt 2 lw lwd1 ps psz title '$\sigma^2_{z, {\mathrm MC}}$',\
     fileMPp using 1:4 with l lc 'black' lw lwd2 dt 6 notitle, \
     fileMPp using 1:7 with points lc rgb 'red' pt 4 lw lwd1 ps psz title '$\sigma^2_{v_z}$', \
     fileMPpMC using 1:7 with points lc rgb 'red' pt 2 lw lwd1 ps psz title '$\sigma^2_{{v_z}, {\mathrm MC}}$',\
     fileMPp using 1:7 with l lc 'red' lw lwd2 dt 6 notitle

# ---------------------------------------------------------------#
#                       SVS, SVD                                 #
# ---------------------------------------------------------------# 
reset
set terminal epslatex standalone size 4, 3 color colortext 10 lw 2 header \
"\\usepackage{amsmath}\n\\usepackage{siunitx}"
set output "SVS.tex"
#set tics font "Times-Roman, 40" 
set xlabel '$k$'
set ylabel '$\Delta v_{\mathrm s}$ (\si{\kilo\m/\s})' 

lwd = 3

set boxwidth 2.5

set grid back lt 1 lc rgb "#e6e6e6" lw 0.25

set key right bottom inside 
set key box spacing 1.8 
set key width 2

set yrange [0:2.2]
set xtics 2
set ytics 0.4
plot SVS using 1:2 with l lc "black" lw lwd title 'null', \
     SVS using 1:3 with l lc "red" lw lwd dt 6 title 'low', \
     SVS using 1:4 with l lc "blue" lw lwd dt 4 title 'medium', \
     SVS using 1:5 with l lw lwd lt 2 dt 5 title 'high'

set output "SVSP0.tex"
set key right top inside 
plot SVSP0 using 1:2 with l lc "black" lw lwd title 'null', \
     SVSP0 using 1:3 with l lc "red" lw lwd dt 6 title 'low', \
     SVSP0 using 1:4 with l lc "blue" lw lwd dt 4 title 'medium', \
     SVSP0 using 1:5 with l lw lwd lt 2 dt 5 title 'high'

set output 'SVD.tex'
set yrange [0:12]
set ytics auto
set ylabel '$\Delta v_{\mathrm d}$ (\si{\kilo\m/\s})' 
set key right bottom inside 

plot SVD using 1:2 with l lc "black" lw lwd title 'null', \
     SVD using 1:3 with l lc "red" lw lwd dt 6 title 'low', \
     SVD using 1:4 with l lc "blue" lw lwd dt 4 title 'medium', \
     SVD using 1:5 with l lw lwd lt 2 dt 5 title 'high'

set output "SVDP0.tex"
plot SVDP0 using 1:2 with l lc "black" lw lwd title 'null', \
     SVDP0 using 1:3 with l lc "red" lw lwd dt 6 title 'low', \
     SVDP0 using 1:4 with l lc "blue" lw lwd dt 4 title 'medium', \
     SVDP0 using 1:5 with l lw lwd lt 2 dt 5 title 'high'

# ---------------------------------------------------------------#
#                       UMC                                      #
# ---------------------------------------------------------------# 
reset
set terminal epslatex standalone size 4, 3 color colortext 10 lw 2 header \
"\\usepackage{amsmath}\n\\usepackage{siunitx}"
set output "DVMC.tex"
#set tics font "Times-Roman, 40" 
set xlabel '$k$'
set ylabel '$\Delta v$ (\si{\kilo\m/\s})' 

lwd = 3

set boxwidth 2.5

set grid back lt 1 lc rgb "#e6e6e6" lw 0.25

plot fileUMC every :::1 using ($1):($2)*29.7578 with l lc "gray" lw lwd notitle,\
     fileUMC every :::::0 using ($1):($2)*29.7578 with l lc "black" lw lwd notitle,\
     0.7619 with l lc "red" lw lwd dt 2 title ""

