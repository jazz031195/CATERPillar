#!/bin/bash -l

cd "/Users/melina/Desktop/EPFL/BachelorProject/Sim_Growth/src/"

search_dir=''
src_files=''
for entry in *.cpp
do
  src_files="$src_files"" ""$entry"
done

command="g++ -I/usr/include -O3 -std=c++11 -std=c++0x -o Growth-Animation.exe $src_files -I. -L/usr/local/lib -lsfml-graphics -lsfml-window -lsfml-system -lGL -lglut -lGLU"

eval "$command"