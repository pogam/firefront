#!/bin/bash
#
python /home/paugam/Src/ForeFire/tools/create_fire_exseg_namelist_dirs.py || exit 1
./run_mesonh_xyz
