reset
set term png
set xlabel "Number of cities"
set ylabel "Time"
set title "TSP solve time by number of cities"
Time="timing.dat"
set output "time.png"
plot Time w l 
