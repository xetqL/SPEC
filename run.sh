#!/bin/bash
set -o xtrace
let sqr=$(($5**2))
rm *.jpg
rm log/* &
rm out.npy &
make
mpirun --oversubscribe -np $1 bin/SPEC $2 $3 $sqr $4 $6
python build_images.py out.npz $4
ffmpeg -framerate 16 -pattern_type glob -i '*water_dummy.jpg'   -c:v libx264 -r 30 -pix_fmt yuv420p -y out.mp4 && vlc out.mp4