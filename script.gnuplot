set output "LiZ.svg"
set key off
set pointsize 0.5
set xlabel "Z"
set ylabel "occupancy of Z"
set term svg
plot [1.9:4.1][1.9:4.1] "Zocc.out" w l
