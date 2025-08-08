#!/bin/bash

module load nco || { echo "Failed to load NCO"; exit 1; }

infile="input/data_mlev.nc"
outdir="split_days"
mkdir -p "$outdir"

# Extract total number of timesteps from the dimension info
nt=$(ncks --trd -m -v valid_time "$infile" | grep 'valid_time dimension' | head -1 | sed -n 's/.*size = \([0-9]\+\).*/\1/p')

if [[ -z "$nt" ]]; then
    echo "Could not determine number of timesteps; aborting."
    exit 1
fi

echo "Number of timesteps: $nt"


for ((i=0; i<nt; i+=4)); do
    epoch=$(ncks -C -v valid_time -d valid_time,$i --trd "$infile" | grep -oP 'valid_time\[\d+\]=\K\d+')
    if [[ -z "$epoch" ]]; then
        echo "Could not extract epoch for timestep $i, skipping."
        continue
    fi
    date_str=$(date -u -d @"$epoch" +%Y%m%d)
    echo "Extracting $date_str ..."
    ncks -O -d valid_time,$i,$((i+3)) "$infile" "$outdir/CAMS_${date_str}.nc"

done
