#!/bin/bash -l
chmod u+x compiler.sh
./compiler.sh
cd /home/localadmin/Documents/CATERPillar/src/
chmod u+x Growth-Animation.exe
./Growth-Animation.exe "/home/localadmin/Documents/CATERPillar/args.conf"
