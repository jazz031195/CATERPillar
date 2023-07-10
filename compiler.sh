#!/bin/bash -l

cd "/Users/melina/Desktop/EPFL/BachelorProject/Sim_Growth/src/"



search_dir=''
src_files=''
for entry in *.cpp
do
  src_files="$src_files"" ""$entry"
done

command="g++ -O3 -std=c++11 -o Growth-Animation $src_files -I. -I/System/Library/Frameworks/OpenGL.framework/Headers -I/usr/X11/include/ -lsfml-graphics -lsfml-window -lsfml-system -framework OpenGL -framework GLUT"
eval "$command"