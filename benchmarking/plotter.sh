#!/bin/bash


DIR=./TMP

SDATA=${DIR}/str_speed_up.dat
SDATA_PLOT=${DIR}/sdata.plt
SDATA_FIG=${DIR}/sdata.png

WDATAQ=${DIR}/wk_efficiency_mean_scalar.dat
WDATAQ_PLOT=${DIR}/wdataq.plt
WDATAQ_FIG=${DIR}/wdataq.png

WDATAQQ=${DIR}/wk_efficiency_mean_sq_scalar.dat
WDATAQQ_PLOT=${DIR}/wdataqq.plt
WDATAQQ_FIG=${DIR}/wdataqq.png

WERRQ=${DIR}/wk_errors_mean_scalar.dat
WERRQ_PLOT=${DIR}/werrq.plt
WERRQ_FIG=${DIR}/werrq.png

WERRQQ=${DIR}/wk_errors_mean_sq_scalar.dat
WERRQQ_PLOT=${DIR}/werrqq.plt
WERRQQ_FIG=${DIR}/werrqq.png


cat <<EOF > ${SDATA_PLOT}
set terminal png  size 1200, 800   
set output  "${SDATA_FIG}" 
# set multiplot layout 1, 2 title ""

set title "tree vs map"
set xlabel "Size"
set ylabel "Time (ms)"
# set logscale y
plot	 "${SDATA}" using 1:2  \
		 with lines lw 3 lc rgb "red" title \
		 "Time for Linked List"

EOF

gnuplot ${SDATA_PLOT}

cat <<EOF > ${WDATAQ_PLOT}
set terminal png  size 1200, 800   
set output  "${WDATAQ_FIG}" 
# set multiplot layout 1, 2 title ""

set title "tree vs map"
set xlabel "Size"
set ylabel "Time (ms)"
# set logscale y
plot	 "${WDATAQ}" using 1:2  \
		 with lines lw 3 lc rgb "red" title \
		 "Time for Linked List"

EOF

gnuplot ${WDATAQ_PLOT}

cat <<EOF > ${WDATAQQ_PLOT}
set terminal png  size 1200, 800   
set output  "${WDATAQQ_FIG}" 
# set multiplot layout 1, 2 title ""

set title "tree vs map"
set xlabel "Size"
set ylabel "Time (ms)"
# set logscale y
plot	 "${WDATAQQ}" using 1:2  \
		 with lines lw 3 lc rgb "red" title \
		 "Time for Linked List"

EOF

gnuplot ${WDATAQQ_PLOT}

cat <<EOF > ${WERRQQ_PLOT}
set terminal png  size 600, 400   
set output  "${WERRQQ_FIG}" 
# set multiplot layout 1, 2 title ""

set title "Error as function of number of MC Steps"
set xlabel "# Monte Carlo Steps"
set ylabel "abs (<o> - 1)"
# set logscale y
plot	 "${WERRQQ}" using (\$1)*2000:2  \
		 with lines lw 3 lc rgb "red" title \
		 "abs (<o> - 1)"

EOF

gnuplot ${WERRQQ_PLOT}


rm ${SDATA_PLOT} ${WDATAQ_PLOT} ${WDATAQQ_PLOT} ${WERRQQ_PLOT} 