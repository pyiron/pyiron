#!/bin/bash
#SBATCH --partition=p.cmfe
#SBATCH --ntasks=40
#SBATCH --time=5760
#SBATCH --output=time.out
#SBATCH --job-name=pi_4639938
#SBATCH --workdir=/cmmc/u/samsstud/RESEARCH/2019/1111/VACANCY/FCC/spx_b_Fe32_3c5_adddf_hdf5/spx_b_Fe32_3c5_adddf
#SBATCH --get-user-env=L

python -m pyiron.base.job.wrappercmd -p /cmmc/u/samsstud/RESEARCH/2019/1111/VACANCY/FCC/spx_b_Fe32_3c5_adddf_hdf5/spx_b_Fe32_3c5_adddf -j 4639938