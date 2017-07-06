set terminal pdf enhanced \
    color size 12.0cm, 9.0cm
#
set output "dos.pdf"
#
set lmargin 12
set bmargin 5
#
set style line 1 lt 1 lw 4 lc rgbcolor "#990099" 
set style line 2 lt 1 lw 2 lc rgbcolor "#009900" pt 1 ps 1.5
set style line 3 lt 1 lw 2 lc rgbcolor "#0000AA" pt 2 ps 1
#
set key at -1.5, 0.03 bottom left spacing 3.0 font 'Cmr10,24'
#
set xtics scale 3.0, 1.5 offset -0.5, -0.5 font 'Cmr10,24'
set xlabel 'Energy [{/Cmmi10 t}]' offset 0.0, -1.0 font 'Cmr10,24'
#
set ytics scale 3.0, 1.5 offset 0.0, -0.5 font 'Cmr10,24'
set ylabel 'DOS [1/{/Cmmi10 t}]' offset -1.5, 0.0 font 'Cmr10,24'
#
plot "dos40.dat" w l ls 1title "Converged", \
     "dos1_8.dat" w p ls 2 title "Linear tetra.", \
     "dos2_8.dat" w p ls 3 title "Optimized tetra"
#
set terminal pdf color enhanced \
    size 12.0cm, 9.0cm
#
set output "lindhard.pdf"
#
set style line 1 lt 1 lw 4 lc rgbcolor "#990099" 
set style line 2 lt 1 lw 2 lc rgbcolor "#009900" pt 1 ps 1.5
set style line 3 lt 1 lw 2 lc rgbcolor "#0000AA" pt 2 ps 1
#
set key at 4.0, 0.9 top right spacing 3.0 font 'Cmr10,18'
#
set xtics scale 3.0, 1.5 offset -0.5, -0.5 font 'Cmr10,24'
set xlabel '{/Cmmi10 q}/{/Cmmi10 k}_F' offset 0.0, -1.0 font 'Cmr10,24'
#
set ytics scale 3.0, 1.5 0.2font 'Cmr10,24'
set ylabel 'Lindhard function' offset -0.5, 0.0 font 'Cmr10,24'
#
plot 0.5+0.5/x*(1-0.25*x**2)*log(abs((x+2)/(x-2))) w l ls 1title "Exact", \
     "lindhard1_8.dat" w p ls 2 title "Linear tetra.", \
     "lindhard2_8.dat" w p ls 3 title "Optimized tetra"
