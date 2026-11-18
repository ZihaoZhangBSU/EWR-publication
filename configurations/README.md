# Model configurations

Configuration files for the East River Noah-MP/iSnobal comparison, retrieved
from the author's Boise State HPC account on 18 November 2026. This directory
contains text settings, an archived launch script, and provenance records.
No simulated data are included.

## Contents and scope

| Location | Purpose and dates |
| --- | --- |
| `Noah-MP/NoahMP_f06/namelist.hrldas` | The requested `NoahMP_f06/full_run` namelist: start 2016-10-01 00:00, hourly forcing/model/output steps, `KHOUR = 85439`. |
| `iSnobal/katana_2017/datakatana_201710.ini` | HRRR/Katana forcing preparation, October 2017. |
| `iSnobal/katana_2017/datakatana_201711.ini` | HRRR/Katana forcing preparation, November 2017. |
| `iSnobal/katana_2017/datakatana_201712.ini` | HRRR/Katana forcing preparation, December 2017. |
| `iSnobal/katana_2017/erw_katana.ini` | Working Katana control file, saved with December 2017 dates; identical to `datakatana_201712.ini`. |
| `iSnobal/katana_2017/katana_erw_run.ini` | WindNinja settings, including a 200 m wind mesh, 5 m output wind height, and a saved 2018-01-01 forecast path. |
| `iSnobal/katana_2017/katana.slurm` | Original Katana SLURM launch script for October-December 2017. |
| `iSnobal/katana_2017/snapshots/run20161001/awsm_config_backup.ini` | Expanded AWSM/iSnobal settings from the first daily WY2017 run, 2016-10-01 00:00-23:00 UTC. |
| `iSnobal/katana_2017/snapshots/run20161001/config.ini` | The corresponding SMRF configuration saved by that daily run. |

The top-level Katana files control forcing preparation and WindNinja; they are
not a complete AWSM/iSnobal physics configuration. The two files under
`snapshots/run20161001/` provide the saved physics and meteorological distribution
settings. Their headers identify AWSM `20260127.dev12+g6ebc94824` and SMRF
`20260310.dev14+gbb8e59064`, generated on 2026-08-08.

The AWSM snapshot uses a project described as 50 m spatial resolution. The
separate 200 m setting in the Katana/WindNinja files refers to the wind mesh.

Only the requested `katana_2017` directory and the Noah-MP f06 namelist are
represented here. The daily AWSM/SMRF files are one reference snapshot, not a
collection of every daily restart configuration or the settings for all later
water years. No model run was performed or scientifically revalidated for this
configuration release.

## Original locations

- iSnobal/Katana: `/bsuscratch/zihaozhang/iSnobal/erw_2nd/katana_2017/`
- Daily AWSM/SMRF snapshot: the same root followed by
  `model_out/erw/wy2017/erw_res_forcing/run20161001/`
- Noah-MP: `/bsuscratch/zihaozhang/Noah-MP-HRRR/NoahMP_f06/full_run/namelist.hrldas`

The source settings, filenames, comments, and original paths are preserved.
The snapshot files are relocated under `snapshots/` for clarity. The byte count,
SHA-256 of each released file, source path, and source-text integrity check are
recorded in [source_manifest.json](source_manifest.json).

## Using these settings elsewhere

Copy the settings to a working run directory and adapt the absolute HPC paths,
scheduler resources, notification address, module/conda paths, and run dates for
that environment. The archived `katana.slurm` retains its original cleanup
commands and mutable working-INI logic; inspect and adapt the script before
running a copy. Its scheduler settings are preserved as provenance, not as a
portable launch prescription.

Supply the required software and external inputs separately: model executables
and parameter tables, HRRR/LDASIN forcing, HRLDAS setup data, topography and
vegetation, WindNinja elevation inputs, and any initial/restart states needed by
the chosen run. Paths to those files do not include their contents in this
release. In particular, `topo.nc`, `topo_windninja_topo.asc`, forcing directories,
and model outputs such as `snow.nc`, `em.nc`, and radiation NetCDF files are
excluded. Output-variable names inside an INI are settings, not released data.

For checksums, run from this directory:

```sh
sha256sum -c SHA256SUMS
```

No new license is assigned by this configuration release; see the repository's
existing licensing statement.
