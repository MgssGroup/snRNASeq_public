#! /bin/bash

## make a folder for bsub_log if not exist
folder="sbatch_log"


if [ ! -d "$folder" ]; then
    mkdir -p "$folder"
fi

snakemake -j 99999 \
    -k \
    --rerun-incomplete \
    --latency-wait 240 \
    --cluster-config cluster.json \
    --cluster "sbatch --mem-per-cpu {cluster.mem} --gres gpu:k80:{cluster.G} --cpus-per-task {cluster.n} --time {cluster.time} --job-name {cluster.name} --output {cluster.output} --error {cluster.error} --account def-cnagy" \
    "$@"
	 
