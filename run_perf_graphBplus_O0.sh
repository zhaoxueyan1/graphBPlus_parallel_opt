#!/bin/bash
root_dir=`pwd`
work_dir="${root_dir}/build"
# source ${root_dir}/../env.sh

bench_lang=cpp
bench_type=perf
declare -a input_files=("graph.csv" "amazonVideo_core5_edges2.csv" "Amazon_Instruments_edges2.csv" "Amazon_Video_edges2.csv" "Amazon_Garden_edges2.csv" "Amazon_Games_edges2.csv" "Amazon_Automotive_edges2.csv" "Amazon_Android_edges2.csv" "Amazon_Outdoors_edges2.csv" "Amazon_Vinyl_edges2.csv" "Amazon_TV_edges2.csv" "Amazon_Jewelry_edges2.csv" "Amazon_Electronics_edges2.csv" "Amazon_Books_edges2.csv")
declare -a output_files=("out.csv" "amazonVideo_core5_out.csv" "Instruments_out.csv" "Video_out.csv" "Garden_out.csv" "Games_out.csv" "Automotive_out.csv" "Android_out.csv" "Outdoors_out.csv" "Vinyl_out.csv" "TV_out.csv" "Jewelry_out.csv" "Electronics_out.csv" "Amazon_Books_out.csv")

for i in "${!input_files[@]}"
do
	log_dir=$root_dir/results_perf/graphBplus_O0_${bench_lang}/${bench_type}/$((i+1))/
    if [ -d $log_dir ]; then
        echo -e "WARN: $log_dir will removed after 5 seconds *******"
        sleep 5
    fi
    rm -rf $log_dir && mkdir -p $log_dir/   
    set -x
    input_file="${input_files[i]}"
    output_file="${output_files[i]}"
    cmd_line="${work_dir}/graphBplus_O0 ${input_file} 100 ${output_file}"
    cd $work_dir
    time perf stat -r 3 \
        -e 'L1-dcache-load-misses,L1-dcache-loads,L1-icache-load-misses,L1-icache-loads' \
        -e 'LLC-load-misses,LLC-loads' \
        -e 'dTLB-load-misses,dTLB-loads,iTLB-load-misses,iTLB-loads' \
        -e 'branch-load-misses,branch-loads,br_mis_pred,br_mis_pred_retired' \
        -e 'instructions,cycles' \
        -o $log_dir/perf.out \
        ${cmd_line} \
            > $log_dir/run.perf.log 2> $log_dir/run.perf.err 
    set +x
done 
