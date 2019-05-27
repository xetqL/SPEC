#!/bin/bash
set -o xtrace
for N in "$@"
do
sed "s/?/$N/g" runall_base.sh > runall_n"$N".sh

sed "s/?/$N/g" runall_base_debug.sh > runall_n"$N"_debug.sh

chmod +x runall_n"$N".sh
chmod +x runall_n"$N"_debug.sh
done
