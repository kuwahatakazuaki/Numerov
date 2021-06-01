#!/bin/bash

if [ "$1" == "" ]; then
  echo "Specify pes file"
  exit 1
fi

input="$1"
output="pes.inp"
var="var.dat"
nx=300 # the number of plots

gnuplot << EOF
  set table "$output"
  set samples $((nx+1))

  stats "$input"
  save variables "$var"

  set format x "%12.5f"
  set format y "%16.9f"
  plot "$input" smooth csplines
EOF

sed -i '1,4d' $output

#xmin=`grep "STATS_min_x" $var | awk '{print $3}'`
#xmax=`grep "STATS_max_x" $var | awk '{print $3}'`
#eI=`grep "STATS_min_y" $var   | awk '{print $3}'`



