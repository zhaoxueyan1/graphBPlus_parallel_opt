#!/bin/bash
root_dir=`pwd`
work_dir="${root_dir}/build"
# source ${root_dir}/../env.sh

bench_lang=cpp
bench_type=perf
graph_type=Amazon_Books
declare -a input_files=("graph.csv" "amazonVideo_core5_edges2.csv" "Amazon_Instruments_edges2.csv" "Amazon_Video_edges2.csv" "Amazon_Garden_edges2.csv" "Amazon_Games_edges2.csv" "Amazon_Automotive_edges2.csv" "Amazon_Android_edges2.csv" "Amazon_Outdoors_edges2.csv" "Amazon_Vinyl_edges2.csv" "Amazon_TV_edges2.csv" "Amazon_Jewelry_edges2.csv" "Amazon_Electronics_edges2.csv" "Amazon_Books_edges2.csv")
declare -a output_files=("out.csv" "amazonVideo_core5_out.csv" "Instruments_out.csv" "Video_out.csv" "Garden_out.csv" "Games_out.csv" "Automotive_out.csv" "Android_out.csv" "Outdoors_out.csv" "Vinyl_out.csv" "TV_out.csv" "Jewelry_out.csv" "Electronics_out.csv" "Amazon_Books_out.csv")

for i in "${!input_files[@]}"
do
    input_file="${input_files[i]}"
    output_file="${output_files[i]}"
    cmd_line="${work_dir}/graphBplus_serial ${root_dir}/data/${input_file} 100 ${root_dir}/data/serial_results/${output_file}"
    cd $work_dir
    time ${cmd_line} > run$((i+1))_serial.O0.out
    set +x
done 

