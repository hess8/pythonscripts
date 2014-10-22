# gss_plot.gp in prepared_input/problem_I_and_II
set terminal postscript landscape color "Helvetica" 24 
set output "gss.eps"
set xlabel "Concentration of hydrogen adatoms"
set ylabel "Form. En. (eV/atom) vs graphane"
set title "Ground State Search for CHW"
set nokey
plot "gss.out" u 2:8 lt 3

!ps2pdf gss.eps && rm gss.eps


