#!/bin/bash

# Average calculator definition -----------------------------------
calculate_average()
{
    # X axis
    x=6
    # Colums to average
    c1=10
    c2=12
    c3=14

    # Calculate
    awk 'BEGIN {
                v1=0;
                v2=0;
                v3=0
         } 
         {
                v1+=$c1; 
                v2+=$c2;
                v3+=$c3
         }
         END {
                printf("%d   %0.5f   %0.5f   %0.5f\n", $x,
                v1/NR,
                v2/NR,
                v3/NR)
        }
        ' "x=$x" "c1=$c1" "c2=$c2" "c3=$c3" \
        $1 1>> $2
}

# Standard Deviation calculator definition ------------------------
calculate_std()
{
    # X axis
    x=6
    # Colums to std
    c1=10
    c2=12
    c3=14
    # Calculate
    awk 'BEGIN { 
                sum1=0; sumsq1=0;
                sum2=0; sumsq2=0;
                sum3=0; sumsq3=0;
         }
         {
                sum1 += $c1; sumsq1 += ($c1)^2;
                sum2 += $c2; sumsq2 += ($c2)^2;
                sum3 += $c3; sumsq3 += ($c3)^2;
         } 
    
         END {
                printf ("%0.5f   %0.5f   %0.5f\n",
                sqrt((sumsq1-sum1^2/NR)/NR),
                sqrt((sumsq2-sum2^2/NR)/NR),
                sqrt((sumsq3-sum3^2/NR)/NR))
        }
        '  "c1=$c1" "c2=$c2" "c3=$c3" \
        $1 1>> $2
}


calculate_speed_up()
{
    
    awk 'BEGIN{ 
        t1 = 0
    } 
    NR==1 {
        t1=$4
    } 
    {
        printf("%d %0.5f\n", $1, t1/$4)
    }' $1 1>> $2 
}

calculate_efficiency()
{
    awk 'BEGIN{ 
        t1 = 0
    } 
    NR==1 {
        t1=$4
    } 
    {
        printf("%d %0.5f\n", $1, t1/($4))
    }' $1 1>> $2 
}


calculate_error()
{
    exact=$3
    awk 'NR==1 {
        printf("%d %0.5f\n", $6, sqrt(($2-exact)^2))
        }' "exact=$exact" $1 1>> $2

}

write_bench_file()
{
    rows=2
    cols=2
    mass=0
    g_cc=0
    observ=$1
    MC_steps=$2
    print_matrix=0

cat << EOF > ${INPUT_FILE}
    ${rows}
    ${cols}
    ${mass}
    ${g_cc}
    ${observ}
    ${MC_steps}
    ${print_matrix} 
EOF

}


# Common definitions -------------------------------------------
EXE=../gn
INPUT_FILE=../benchmark.inp

DIR=./TMP
DAT=.dat
RAW=.raw

REPETITIONS=5
CPUS=`seq 1 1 4`

# Strong Scaling Variable definitions --------------------------
STR_DATA=${DIR}/str_time
STR_AVG_DATA=${DIR}/str_avg_time
STR_STD_DATA=${DIR}/str_std_time
STR_MERGED_DATA=${DIR}/str_avg_std_time
STR_SPEED_UP=${DIR}/str_speed_up

# Weak Scaling Variable definitions -----------------------------
WK_DATA=${DIR}/wk_time
WK_AVG_DATA=${DIR}/wk_avg_time
WK_STD_DATA=${DIR}/wk_std_time
WK_MERGED_DATA=${DIR}/wk_avg_std_time
WK_ERROR=${DIR}/wk_errors   
WK_EFFICIENCY=${DIR}/wk_efficiency

# SIZES=(500 1000 1500 2000) 
SIZES=(2000 4000 6000 8000) 
OBSRV=(mean_scalar mean_sq_scalar)
EXPECT_OBSR=(0 1)

# Writing input file --------------------------------------------
write_bench_file mean_sq_scalar 1000

# Remove previous data ------------------------------------------
rm ${DIR}/*.dat ${DIR}/*raw 2> /dev/null

# Strong Scaling Running ----------------------------------------

for cores in ${CPUS[*]}; do
    for repetition in `seq 1 ${REPETITIONS}`; do
        mpirun -np ${cores} ${EXE} ${INPUT_FILE}                \
         1>> ${STR_DATA}_${cores}${RAW}     
    done
    
    calculate_average ${STR_DATA}_${cores}${RAW}                \
                          ${STR_AVG_DATA}${RAW}
    calculate_std     ${STR_DATA}_${cores}${RAW}                \
                          ${STR_STD_DATA}${RAW}
    rm ${STR_DATA}_${cores}${RAW}

done

paste   ${STR_AVG_DATA}${RAW} ${STR_STD_DATA}${RAW}             \
          1>> ${STR_MERGED_DATA}${DAT}

calculate_speed_up ${STR_MERGED_DATA}${DAT}                     \
                   ${STR_SPEED_UP}${DAT}
                 
rm  ${STR_AVG_DATA}${RAW} ${STR_STD_DATA}${RAW}

# Weak Scaling Running -----------------------------------------

expec_indx=0
for obsrv in ${OBSRV[*]}; do

    size_indx=0
    for cores in ${CPUS[*]}; do
        
        write_bench_file ${obsrv} ${SIZES[size_indx]}
        ((size_indx += 1))

        for repetition in `seq 1 ${REPETITIONS}`; do
            mpirun -np ${cores} ${EXE} ${INPUT_FILE}            \
             1>> ${WK_DATA}_${obsrv}_${cores}${RAW}     
        done
        
        calculate_error ${WK_DATA}_${obsrv}_${cores}${RAW}      \
                            ${WK_ERROR}_${obsrv}${DAT}          \
                            ${EXPECT_OBSR[expec_indx]}

        calculate_average ${WK_DATA}_${obsrv}_${cores}${RAW}    \
                              ${WK_AVG_DATA}_${obsrv}${RAW}
        calculate_std     ${WK_DATA}_${obsrv}_${cores}${RAW}    \
                              ${WK_STD_DATA}_${obsrv}${RAW}
        rm ${WK_DATA}_${obsrv}_${cores}${RAW}

    done
    ((expec_indx += 1))

    paste   ${WK_AVG_DATA}_${obsrv}${RAW}                        \
            ${WK_STD_DATA}_${obsrv}${RAW}                        \
              1>> ${WK_MERGED_DATA}_${obsrv}${DAT}

    calculate_efficiency ${WK_MERGED_DATA}_${obsrv}${DAT}        \
                 ${WK_EFFICIENCY}_${obsrv}${DAT}

    rm  ${WK_AVG_DATA}_${obsrv}${RAW}                            \
        ${WK_STD_DATA}_${obsrv}${RAW}
done



# -----------------------------------------------------------------