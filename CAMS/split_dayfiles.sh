#!/bin/bash

module load nco || { echo "Failed to load NCO"; exit 1; }

infile_mlev="input/data_mlev.nc"
infile_sfc="input/data_sfc.nc"
outdir="split_days"
mkdir -p "$outdir"

# Extract total number of timesteps (assumed same for both files)
nt=$(ncks --trd -m -v valid_time "$infile_mlev" | grep 'valid_time dimension' | head -1 | sed -n 's/.*size = \([0-9]\+\).*/\1/p')

if [[ -z "$nt" ]]; then
    echo "Could not determine number of timesteps; aborting."
    exit 1
fi

echo "Number of timesteps: $nt"

for ((i=0; i<nt; i+=4)); do
    # Extract epoch time for naming
    epoch=$(ncks -C -v valid_time -d valid_time,$i --trd "$infile_mlev" | grep -oP 'valid_time\[\d+\]=\K\d+')
    if [[ -z "$epoch" ]]; then
        echo "Could not extract epoch for timestep $i, skipping."
        continue
    fi
    date_str=$(date -u -d @"$epoch" +%Y%m%d)

    echo "Processing $date_str ..."

    # Temporary split files
    tmp_mlev="${outdir}/tmp_CAMS_${date_str}_mlev.nc"
    tmp_sfc="${outdir}/tmp_CAMS_${date_str}_sfc.nc"
    final_out="${outdir}/CAMS_${date_str}.nc"

    # Extract 4 timesteps from mlev and sfc files
    ncks -O -d valid_time,$i,$((i+3)) "$infile_mlev" "$tmp_mlev"
    ncks -O -d valid_time,$i,$((i+3)) "$infile_sfc" "$tmp_sfc"

    # Merge them into one file
    ncks -A "$tmp_sfc" "$tmp_mlev"  # Append sfc vars into mlev file

    # Move merged file to final output
    mv "$tmp_mlev" "$final_out"

    # Remove temp sfc split file
    rm -f "$tmp_sfc"

    echo "Created $final_out"
done
