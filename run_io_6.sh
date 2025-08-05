#!/bin/bash
#SBATCH --job-name=hemelb_io_test
#SBATCH --output=hemelb-io-%j.out
#SBATCH --time=01:00:00
#SBATCH --exclusive
#SBATCH --nodes=6
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --account=m24oc-s2450341
#SBATCH --partition=standard
#SBATCH --qos=standard


# environment setup
module use /mnt/lustre/e1000/home/y07/shared/cirrus-modulefiles
module load epcc/setup-env
module load intel-compilers
# source /work/m24oc/m24oc/s2450341/hemelb/setup_hemeLB.sh > /dev/null 2>&1 && echo ">>> Environment loaded!"
echo ">>> Configuring environment..."
source /work/m24oc/m24oc/s2450341/hemelb/setup_hemeLB.sh
# source /work/m24oc/m24oc/s2450341/hemelb/setup_hemeLB.sh > /dev/null 2>&1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/compilers_and_libraries/linux/lib/intel64

echo ">>> Environment loaded (errors ignored)"
set -e
# source /work/m24oc/m24oc/s2450341/hemelb/setup_hemeLB.sh > /dev/null 2>&1 && echo ">>> Environment loaded!" || echo "Environment failed!"


# parameter settings
EXEC=/work/m24oc/m24oc/s2450341/hemelb/build/hemelb
CONFIG_DIR=/work/m24oc/m24oc/s2450341/hemelbTest-data/tests/gmy/bifurcation-40um
CONFIG=$CONFIG_DIR/config.xml

cd $CONFIG_DIR

# clean old outputs and create new output directories
# rm -rf outputs/
for N in 6; do
  for P in 6 8 12 18 16 32 36 64; do
    [[ $N -eq 8 && $P -eq 4 ]] && continue
    mkdir -p outputs/n${N}_p${P}
  done
done

# loop through different configurations
for N in 6 ; do
  for P in 6 8 12 16 18 32 64; do
    [[ $N -eq 8 && $P -eq 4 ]] && continue

    echo ">>> Running n=${N} p=${P} at $(date)"
    rm -rf results

    srun --partition=standard --nodes=${N} --ntasks=${P} $EXEC -in $CONFIG > outputs/n${N}_p${P}/log.txt 2>&1

    if [ -f results/Extracted/flow_snapshot.h5 ]; then
      cp results/Extracted/flow_snapshot.h5 outputs/n${N}_p${P}/
      cp results/Extracted/flow_snapshot.xmf outputs/n${N}_p${P}/
      cp results/report.* outputs/n${N}_p${P}/
    else
      echo "Output not found for n=${N} p=${P}" >> outputs/n${N}_p${P}/log.txt
    fi
  done
done

# extract all writing times
grep "Extraction writing" outputs/*/report.txt > outputs/summary.txt || echo " No report data found."

echo ">>> Finished all runs"
