snakemake all --jobname "{jobid}" --jobs 100 \
	--rerun-triggers mtime \
	--keep-going \
	--rerun-incomplete \
	--snakefile Snakefile \
	--use-conda \
	--use-envmodules \
	--printshellcmds \
	--latency-wait 20 \
	--cluster-config config_cluster.yaml \
	--cluster "sbatch --job-name {cluster.name} --output {cluster.output} --time {cluster.time} --mem {cluster.mem} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --parsable --partition {cluster.partition} --gres {cluster.gres}" \
	> logs/snakemake.log 2>&1 &
