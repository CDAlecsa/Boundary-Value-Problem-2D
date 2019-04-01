#load "F.gnu" 
set terminal postscript portrait enhanced color 'Helvetica' 14
set out "Fexp2D.eps"
 set xlabel "x";
 set ylabel "y"
 set zlabel "F(x,y)";
 set zlabel offset 4,5,4
 set size 1,0.5
 set origin 0,0
# set key out
#set key at graph 0.999,0.999,0.5 spacing  1

set dgrid3d 10,10
#set ticslevel 2.0
set xtics; set ytics; set ztics;
set xtics 2
set ytics 2
#set ztics 10
#set grid
set view 60,60
set clabel
#set contour both
set contour
set cntrparam levels 12
splot [][][] 'f' w l lt 1 

set terminal x11
#set terminal windows
 replot
 set xlabel ; set ylabel ; set nocontour ; set key ; set nologscale
set size; set origin; set title;
 unset xtics; unset ytics; unset ztics; unset contour; unset clabel; unset key;

 set title
