set key top left

fold='pentagon/data/'
set key autotitle columnheader
set style data lines

plot fold.'9_1.tsv' u 1:($5/8504) t columnhead(3)
replot fold.'9_2.tsv' u 1:($5/1567) t columnhead(3)
replot fold.'9_3.tsv' u 1:($5/205) t columnhead(3)
replot fold.'9_4.tsv' u 1:($5/35) t columnhead(3)
replot fold.'9_5.tsv' u 1:($5/5) t columnhead(3)

pause -1
