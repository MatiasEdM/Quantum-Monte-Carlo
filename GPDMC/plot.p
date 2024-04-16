set   autoscale                      
unset log                           
unset label                        
set xtic auto                    
set ytic auto                     
set title  "ENERGY ALONG THE PATH"
set xlabel "Time [a.u.]"
set ylabel "Energy [a.u.]"
set terminal png
set output 'Energy.png'

plot 'energy.out' using 1:2 title 'E [Hartree]' with lines lt rgb 'black', 'energy.out' using 1:3 title 'T [Hartree]' with lines lt rgb 'red', 'energy.out' using 1:4 title 'Vee [Hartree]' with lines lt rgb 'blue', 'energy.out' using 1:5 title 'VeN [Hartree]' with lines lt rgb 'green' 
