set terminal latex
set output "tsp40it.tex"
plot "tsp40it.tsp" with lines,\
     "tsp40it.tsp" with points pt 3 lc rgb 'black'
