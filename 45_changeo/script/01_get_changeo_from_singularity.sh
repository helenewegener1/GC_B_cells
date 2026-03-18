
# singularity is an alternative to docker that is HPC friendly

module load tools singularity/4.3.0

export SINGULARITY_CACHEDIR=/scratch/helweg/singularity/.singularity_cache
export SINGULARITY_TMPDIR=/scratch/helweg/singularity/.singularity_tmp

mkdir -p $SINGULARITY_CACHEDIR $SINGULARITY_TMPDIR


export SINGULARITY_BIND="/home/projects/dtu_00062:/home/projects/dtu_00062,/scratch/helweg:/scratch/helweg"

IMAGE="/scratch/helweg/singularity/immcantation_suite-4.7.0.sif"
singularity pull $IMAGE docker://immcantation/suite:4.7.0
