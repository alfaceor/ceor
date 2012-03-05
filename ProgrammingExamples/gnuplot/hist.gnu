reset
unset key; set xtics nomirror; set ytics nomirror; set border front;
div=1.1; bw = 0.9; h=1.0; BW=0.9; wd=10; LIMIT=255-wd; white = 0
red = "#080000"; green = "#000800"; blue = "#000008"
set auto x
set yrange [0:11]
set style data histogram
set style histogram cluster gap 1
set style fill solid
set boxwidth bw
set multiplot
plot 'hist.dat' u 1 lc rgb red, '' u 2 lc rgb green, '' u 3 lc rgb blue
unset border; set xtics format " "; set ytics format " "; set ylabel " "
call 'hist_r.gnu'
unset multiplot

