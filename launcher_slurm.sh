#! /bin/bash
#SBATCH -A ghwtcmr
#SBATCH -p defq
#SBATCH --mem=8G
eval "$(micromamba shell hook --shell bash)"
micromamba activate nextflow
export NFX_OPTS="-Xms=512m -Xmx=8g -Djava.io.tmpdir=${TMPDIR}"
nextflow run kallistoViral.nf -config rocket.config -params-file params.json