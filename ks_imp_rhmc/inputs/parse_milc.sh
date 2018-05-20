#!/bin/bash
#********* MILC OUTPUT PARSING SCRIPT*********#

usage() { echo "Usage: $0 [-l <baseline|qphix|quda>] [-m <flops|time>] [-r <result-file>" 1>&2; exit 1; }

while getopts ":l:m:r:" o; do
    case "${o}" in
        l)
            l=${OPTARG}
	    if ! [[ "$l" == "qphix" || "$l" == "quda" || "$l" == "baseline" ]];then usage; exit 1; fi
            ;;
        m)
            m=${OPTARG}
	    if ! [[ "$m" == "flops" || "$m" == "time" ]];then usage; exit 1; fi
            ;;

        r)
            r=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${l}" ] || [ -z "${m}" ] || [ -z "${r}" ]; then
    usage
fi

echo "LIBRARY = ${l}"
echo "METRIC = ${m}" 
echo "RESULT FILE = ${r}" 

#Get Number of MPI Ranks to account for FLOPs correctly#
NODES=`awk '/portable/ { nodes=$6 } END {print nodes}' ${r}`

if [ "${m}" == "time" ]; then

	echo "Parsing Time for different solvers"
	
##QPHIX##
	if [ "${l}" == "qphix" ]; then

	echo "Assuming output enabled with QPhiX"
	CGTIME_PER_ITER=`awk '/^CONGRAD5:.*multicg_offset_QPHIX/{cgtime+=$4/$12;++count}END{print cgtime/count}' ${r}`
	CGTIME_PER_ITER_NR=`awk '/^CONGRAD5-QPhiX:.*multicg_offset_QPhiX/{cgtime_noremap+=$4/$12;++count}END{print cgtime_noremap/count}' ${r}`
#	CGTIME_PER_ITER_NR=`awk '/^CONGRAD5-QPhiX:.*multicg_offset_QPHIX/{cgtime_noremap=$4;}/^CONGRAD5-QPhiX:.*multicg_offset_QPHIX/{iters=$12;cgtime_nr+=cgtime_noremap/$12;++count}END{print cgtime_nr/count}' ${r}`
 
	GFTIME_PER_CALL=`awk '/^GFTIME:/{++count;if (count > 1) gftime+=$4}END{print gftime/(count-1)}' ${r}`
	GFTIME_PER_CALL_NR=`awk '/^GFTIME_QPHIX:/{++count;if (count > 1) gftime_noremap+=$4}END{print gftime_noremap/(count-1)}' ${r}`

	FFTIME_PER_CALL=`awk '/^FFTIME:/{++count;if (count > 1) fftime+=$4}END{print fftime/(count-1)}' ${r}`
	FLTIME_PER_CALL=`awk '/^FLTIME:/{++count;if (count > 1) fltime+=$4}END{print fltime/(count-1)}' ${r}`


	echo "CGTIME_PER_ITER = $CGTIME_PER_ITER"
	echo "CGTIME_PER_ITER (NO REMAP) = $CGTIME_PER_ITER_NR"

	echo "GFTIME_PER_CALL = $GFTIME_PER_CALL"
	echo "GFTIME_PER_CALL (NO REMAP)= $GFTIME_PER_CALL_NR"

	echo "FFTIME_PER_CALL = $FFTIME_PER_CALL"
	echo "FLTIME_PER_CALL = $FLTIME_PER_CALL"

##QUDA##
	elif [ "${l}" == "quda" ]; then

	echo "Assuming output enabled with QUDA on GPU"

	CGTIME_PER_ITER=`awk '/multicg_offset_QUDA/{cgtime+=$4/$12;++count}END{print cgtime/count}' ${r}`
	GFTIME_PER_CALL=`awk '/^GFTIME:/{++count;if (count > 1) gftime+=$4}END{print gftime/(count-1)}' ${r}`
	FFTIME_PER_CALL=`awk '/^FFTIME:/{++count;if (count > 1) fftime+=$4}END{print fftime/(count-1)}' ${r}`
	FLTIME_PER_CALL=`awk '/^FLTIME:/{++count;if (count > 1) fltime+=$4}END{print fltime/(count-1)}' ${r}`

	echo "CGTIME_PER_ITER = $CGTIME_PER_ITER"
	echo "GFTIME_PER_CALL = $GFTIME_PER_CALL"
	echo "FFTIME_PER_CALL = $FFTIME_PER_CALL"
	echo "FLTIME_PER_CALL = $FLTIME_PER_CALL"

##BASELINE##
	elif [ "$l" = "baseline" ]; then

	echo "Assuming output enabled with QPhiX"
	CGTIME_PER_ITER=`awk '/multicg_offset D/{cgtime+=$4/$12;++count}END{print cgtime/count}' ${r}`
	GFTIME_PER_CALL=`awk '/^GFTIME:/{++count;if (count > 1) gftime+=$4}END{print gftime/(count-1)}' ${r}`
	FFTIME_PER_CALL=`awk '/^FFTIME:/{++count;if (count > 1) fftime+=$4}END{print fftime/(count-1)}' ${r}`
	FLTIME_PER_CALL=`awk '/^FLTIME:/{++count;if (count > 1) fltime+=$4}END{print fltime/(count-1)}' ${r}`

	echo "CGTIME_PER_ITER = $CGTIME_PER_ITER"
	echo "GFTIME_PER_CALL = $GFTIME_PER_CALL"
	echo "FFTIME_PER_CALL = $FFTIME_PER_CALL"
	echo "FLTIME_PER_CALL = $FLTIME_PER_CALL"
	#elif [ "$l" = "qopqdp" ]; then

	#else
	fi
elif [ "${m}" == "flops" ]; then

	echo "Parsing FLOPs for different solvers"
##QPHIX##
	if [ "${l}" == "qphix" ]; then

		echo "Assuming output enabled with QPhiX"
		CGFLOPS=`awk '/^CONGRAD5:.*multicg_offset_QPHIX/ { cgflops+=$NF;++count } END {print cgflops/count}' ${r}`
		CGFLOPS_NR=`awk '/^CONGRAD5-QPhiX:.*multicg_offset_QPhiX/ {cgflops_noremap+=$NF;++count } END {print cgflops_noremap/count}' ${r}`

		GFFLOPS=`awk '/^GFTIME:/ { gfflops+=$NF;++count } END {print gfflops/count}' ${r}`
		GFFLOPS_NR=`awk '/^GFTIME_QPHIX:/ { gfflops_noremap+=$NF;++count } END {print gfflops_noremap/count}' ${r}`

		FFFLOPS=`awk '/^FFTIME/ { ffflops+=$NF;++count } END {print ffflops/count}' ${r}`
		FLFLOPS=`awk '/^FLTIME/ { flflops+=$NF;++count } END {print flflops/count}' ${r}`

		CGFLOPS=`echo "$CGFLOPS*$NODES*0.001" | bc -l`
		CGFLOPS_NR=`echo "$CGFLOPS_NR*$NODES*0.001" | bc -l`

		GFFLOPS=`echo "$GFFLOPS*$NODES*0.001" | bc -l`
		GFFLOPS_NR=`echo "$GFFLOPS_NR*$NODES*0.001" | bc -l`

		FFFLOPS=`echo "$FFFLOPS*$NODES*0.001" | bc -l`
		FLFLOPS=`echo "$FLFLOPS*$NODES*0.001" | bc -l`

		echo "TOTAL RANKS = $NODES"
		echo "CG GFLOPS = $CGFLOPS"
		echo "CG GFLOPS (w/o REMAP) = $CGFLOPS_NR" 
		echo "GF GFLOPS = $GFFLOPS"
		echo "GF GFLOPS (w/o REMAP) = $GFFLOPS_NR"
		echo "FF GFLOPS = $FFFLOPS"
		echo "FL GFLOPS = $FLFLOPS"
##QUDA##
	elif [ "${l}" == "quda" ]; then

		echo "Assuming output enabled with QUDA on GPU"
		CGFLOPS=`awk '/multicg_offset_QUDA/ { cgflops+=$15;++count } END {print cgflops/count}' ${r}`
		GFFLOPS=`awk '/^GFTIME:/ { ++count;if (count > 1) gfflops+=$8} END {print gfflops/(count-1)}' ${r}`
		FLFLOPS=`awk '/^FLTIME:/ { ++count;if (count > 1) flflops+=$10} END {print flflops/(count-1)}' ${r}`

		CGFLOPS=`echo "$CGFLOPS*$NODES*0.001" | bc -l`
		GFFLOPS=`echo "$GFFLOPS*$NODES*0.001" | bc -l`
		FLFLOPS=`echo "$FLFLOPS*$NODES*0.001" | bc -l`

		echo "TOTAL RANKS = $NODES"
		echo "CG GFLOPS = $CGFLOPS"
		echo "GF GFLOPS = $GFFLOPS"
		echo "FL GFLOPS = $FLFLOPS"

##BASELINE##
	elif [ "${l}" == "baseline" ]; then

		echo "Assuming baseline MILC output"
		CGFLOPS=`awk '/multicg_offset D/ { cgflops+=$NF;++count } END {print cgflops/count}' ${r}`
		GFFLOPS=`awk '/^GFTIME:/ { gfflops+=$NF;++count } END {print gfflops/count}' ${r}`
		FFFLOPS=`awk '/^FFTIME/ { ffflops+=$NF;++count } END {print ffflops/count}' ${r}`
		FLFLOPS=`awk '/^FLTIME/ { flflops+=$NF;++count } END {print flflops/count}' ${r}`

		CGFLOPS=`echo "$CGFLOPS*$NODES*0.001" | bc -l`
		GFFLOPS=`echo "$GFFLOPS*$NODES*0.001" | bc -l`
		FFFLOPS=`echo "$FFFLOPS*$NODES*0.001" | bc -l`
		FLFLOPS=`echo "$FLFLOPS*$NODES*0.001" | bc -l`

		echo "TOTAL RANKS = $NODES"
		echo "CG GFLOPS = $CGFLOPS"
		echo "GF GFLOPS = $GFFLOPS"
		echo "FF GFLOPS = $FFFLOPS"
		echo "FL GFLOPS = $FLFLOPS"

	#elif [ "$l" = "qopqdp" ]; then

	#else

	fi

fi

