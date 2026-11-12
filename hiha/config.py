import os
os.environ["R_HOME"] = "R/environment/location"

from pathlib import Path

## Directories
output_dir = Path("output")
if not output_dir.exists():
    output_dir.mkdir()

image_dir = output_dir.joinpath("images")
if not image_dir.exists():
    image_dir.mkdir()

## Filenames
density_fname = output_dir.joinpath("hiha_density.RDS")
phm_fname = output_dir.joinpath("hiha_phm.pkl")
phm_time_fname = output_dir.joinpath("hiha_phm_time.pkl")
hiha_fname = output_dir.joinpath("human_immune_health_atlas_full.h5ad")