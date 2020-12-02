    reset
    set fontpath "/usr/share/fonts/truetype/lyx"
    set terminal pngcairo enhanced color size 900,400 #transparent
    set font 'cmmi10.ttf'
    set font 'cmr10.ttf'
    set decimalsign ','

    outfile1 = sprintf('spacetime01.png')
    outfile2 = sprintf('shot01.png')
    datfile1 = sprintf('at01.dat')

    outfile3 = sprintf('spacetimelapl01.png')
    outfile4 = sprintf('laplshot01.png')
    datfile2 = sprintf('lapl01.dat')

    set border linewidth 1.5
    set lmargin at screen 0.18
    set rmargin at screen 0.90
    set bmargin at screen 0.30
    set tmargin at screen 0.90
#
#===================SPATIAL_PROFILE==================================================
#
    set xlabel '{/Cmmi10=42 i}' offset 0.0,-2.0
    set ylabel '{/Cmmi10=42 x_{i}}' offset -3.2,0.0
    set format x "{%.0f}"
    set format y "{%.0f}"
    set xrange [0.0:900]
    set yrange [-3.0:3.0]
    set xtics  0,300,900 scale 2,1 font 'Cmr10, 42' offset 0,0 out mirror
    set xtics add ("100" 300, "200" 600, "300" 900)
    set ytics -3.,3.0,3.0 scale 2,1 font 'Cmr10, 42' offset 0,0 out mirror
    set mxtics 4
    set mytics 3

    set style line 2 lc rgb 'black' pt 7 ps 0.1   # circle

    set output outfile2

# Plot  %02.0f.dat
    plot datfile1 matrix every 3:40::0::1499 with lines lt 1. lc 'black' notitle#, \
         #datfile1 matrix every 3:40::0::1499 with points lt 6. ps 0.5 pt 7. notitle
#
#===================SPATIAL_PROFILE==================================================
#
    set xlabel '{/Cmmi10=42 i}' offset 0.0,-2.0
    set ylabel '{/Cmmi10=42 u_{i}}' offset -3.2,0.0
    set format x "{%.0f}"
    set format y "{%.0f}"
    set xrange [0.0:900]
    set yrange [-0.05:5.0]
    set xtics  0,300,900 scale 2,1 font 'Cmr10, 42' offset 0,0 out mirror
    set xtics add ("100" 300, "200" 600, "300" 900)
    set ytics -4.,4.0,4.0 scale 2,1 font 'Cmr10, 42' offset 0,0 out mirror
    set mxtics 4
    set mytics 4

    set style line 2 lc rgb 'black' pt 7 ps 0.1   # circle

    set output outfile4

# Plot  %02.0f.dat
    plot datfile2 matrix every 3:40::0::1499 with lines lt 1. lc 'black' notitle#, \
         #datfile2 matrix every 3:40::0::1499 with points lt 6. ps 0.5 pt 7. notitle