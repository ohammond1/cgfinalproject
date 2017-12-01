#!/usr/bin/env bash

if [[ $# < 3 ]]; then
    echo "Usage: $0 program_name genome_file sequence"
    exit 1
fi

program_name=$1
genome_file=$2
sequence=$3



for i in {1..3..1}; 
do ./$program_name $genome_file $sequence; done


