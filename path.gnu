reset
set term png
set xlabel "x"
set ylabel "y"
set title "Shortest path to all cities"
Path="plot.dat"
set output "path.png"
plot Path with linespoints pointtype 7 pointsize 2
