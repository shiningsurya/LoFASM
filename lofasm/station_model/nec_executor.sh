#!/usr/bin/env bash
set -x

if [[ $# -ne 1 ]]; then
		echo "$(basename $0) <freq>"
		exit 0
fi

FREQ=$1

./lwa1_generator.sh "alpha_x_p_$FREQ.nec" $FREQ X P
./lwa1_generator.sh "alpha_y_p_$FREQ.nec" $FREQ Y P
./lwa1_generator.sh "alpha_x_t_$FREQ.nec" $FREQ X T
./lwa1_generator.sh "alpha_y_t_$FREQ.nec" $FREQ Y T

parallel -j4  --bar nec2++ -i {} -o {.}.out ::: alpha_?_?"_$FREQ.nec"
