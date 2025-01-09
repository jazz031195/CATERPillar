#!/bin/bash -l

cd "/home/localadmin/Documents/CATERPillar/src/"

src_files=''
for entry in *.cpp
do
  src_files="$src_files"" ""$entry"
done

#  sbatch "$entry"

command="g++ -O3 -std=c++11 -lpthread -std=c++0x -pthread -w -I.$src_files -o Growth-Animation.exe"


eval "$command"