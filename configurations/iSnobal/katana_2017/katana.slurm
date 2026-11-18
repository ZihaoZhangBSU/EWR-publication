#!/bin/bash
#SBATCH --job-name=katana-2017
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=12
#SBATCH --mem=512
#SBATCH --partition=short           

#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=zihaozhang@u.boisestate.edu
#SBATCH --output=slurm-%j.out-%N
#SBATCH --error=slurm-%j.err-%N
#SBATCH --chdir=/bsuhome/zihaozhang

module use /bsushare/hpmarshall-shared/spack/tcl/linux-centos7-cascadelake
module load windninja/3.11.2-gcc-12.1.0-up4e7sr

START_YEAR=2017
START_MONTH=10


KATANA_OUTPUT_DIR=/bsuscratch/zihaozhang/iSnobal/erw_2nd/katana_2017/data
KATANA_INI=/bsuscratch/zihaozhang/iSnobal/erw_2nd/katana_2017/erw_katana.ini

for MONTH in {0..2}; do  # for the next year to 9/30 should be 0..8
  START=($(date -d "${START_MONTH}/01/${START_YEAR} + ${MONTH} month" "+%m %Y"))
  END=($(date -d "${START_MONTH}/01/${START_YEAR} + $((MONTH + 1)) month" "+%m %Y"))
  MONTH=${START[0]}
  YEAR=${START[1]}
  MONTH_END=${END[0]}
  YEAR_END=${END[1]}
  echo "Processing: ${YEAR}-${MONTH}"
  # Update the config to the month processed
  sed -i -e "s/^start_date: .* UTC/start_date: ${YEAR}-${MONTH}-01 00:00:00 UTC/g" ${KATANA_INI}
  sed -i -e "s/^end_date: .* UTC/end_date: ${YEAR_END}-${MONTH_END}-01 00:00:00 UTC/g" ${KATANA_INI}
  
  
  
  cp ${KATANA_INI} ${KATANA_OUTPUT_DIR}katana_${YEAR}${MONTH}.ini

  ${HOME}/miniforge3/bin/mamba run -p "/bsushare/hpmarshall-shared/conda-envs/katana" \
    run_katana ${KATANA_INI}

  if [ $? -ne 0 ]; then
    >&2 echo "ERROR processing ${YEAR}${MONTH}"
    continue
  fi
done

find ${KATANA_OUTPUT_DIR} -name *.prj -type f -delete
find ${KATANA_OUTPUT_DIR} -name *_cld.asc -type f -delete

find ${KATANA_OUTPUT_DIR} -name *.grib2 -type f -delete
find ${KATANA_OUTPUT_DIR} -name hrrr.* -type d -empty -delete
