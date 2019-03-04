#!/bin/bash
set -o xtrace
let sqr=$(($4**2))
let PENUM=$(($1*$2))
rm *.jpg
rm log/* &
rm out.npy &
make
mpirun --oversubscribe -np $PENUM bin/SPEC $1 $2 $sqr $3 $5
python build_images.py out.npz $3
ffmpeg -framerate 10 -pattern_type glob -i '*water_dummy.jpg'   -c:v libx264 -r 30 -pix_fmt yuv420p -y out.mp4 && vlc out.mp4
