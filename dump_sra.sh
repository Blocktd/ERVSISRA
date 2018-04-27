#!/bin/bash

## Author: Spencer Richman
## Download an SR archive, record version number, time, etc.

acc_num=$1

file_dir="files_$acc_num"
[ -d $file_dir ] || mkdir $file_dir

echo "Dumping fasta files, this will take a while..."

fastq-dump -v -B -W --fasta --split-files --skip-technical --read-filter pass -O $file_dir $acc_num 

echo "...done"
