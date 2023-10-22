
filefullsol1 = indir."fullsol1.dat"
filefullsol2 = indir."fullsol2.dat"
fileDVs1 = indir."DVs1.dat"
fileDVs2 = indir."DVs2.dat"
fileELLs1 = indir."covEllipses1.dat"
fileELLs2 = indir."covEllipses2.dat"

trajE = indir."trajEfile.dat"

# ---------------------------------------------------------------#
#                             trajn                              #
# ---------------------------------------------------------------# 
set terminal epslatex standalone size 4, 4 color colortext 10 lw 2  header \
"\\usepackage{amsmath}\n\\usepackage{siunitx}"
set output "trajn.tex"
set size ratio -1 # axis equal

stats indir."fullsol2.dat"
N = STATS_records

set xrange [-1.2:1.7]
set yrange [-1.2:1.2]
set xlabel "$x$"
set ylabel "$y$" 

rf = 1.6507
chi_sq = 1*sqrt(5.991)
lwd = 1.5
lwDV = 1

scalaP = chi_sq
scalaV = 40
e_step = 1

set grid back lt 1 lc rgb "#e6e6e6" lw 0.25
set parametric
set trange [0:2*pi]
set style ellipse units xx

plot 0, 0 with points lc rgb 'orange' pt 7 ps 1.5 notitle, \
     trajE using 2:3 with l notitle lw lwd lc "blue", \
     filefullsol1 using 2:3 with l lc "dark-magenta" lw lwd notitle,\
     fileDVs1 using 2:3:($5)*scalaV:($6)*scalaV with vectors filled head lc "dark-magenta" lw 4 notitle, \
     fileELLs1 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "dark-magenta" lw lwDV notitle, \
     filefullsol2 using 2:3 with l lc "dark-green" lw lwd notitle,\
     fileDVs2 using 2:3:($5)*scalaV:($6)*scalaV with vectors filled head lc "dark-green" lw 4 notitle, \
     fileELLs2 using ($1)+($3)*scalaP:($2) + ($4)*scalaP w l lc "dark-green" lw lwDV notitle, \
     filefullsol1 using 2:3 every ::1::1 with points lc 'dark-magenta' pt 8 lw lwd ps 2 title 'Departure',\
     filefullsol2 using 2:3 every ::1::1 with points lc 'black' pt 4 lw lwd ps 2 title 'Flyby',\
     filefullsol2 using 2:3 every ::N-2::N-2 with points lc 'dark-green' pt 6 lw lwd ps 2 title 'Arrival'
     #trajE every ::0::0 using 2:3 with points lc rgb 'blue' pt 4 lw 1.5 ps 1.5 notitle,\
     #trajM every ::0::0 using 2:3 with points lc rgb 'brown' pt 3 lw 1.5 ps 1.5 notitle,\
     #trajM every ::0::0 using 2:3 with points lc rgb 'brown' pt 3 lw 1.5 ps 1.5 notitle,\
     ##'-' with points lc rgb 'brown' pt 3 lw 1.5 ps 1.5 notitle
     ##-1.2 0



