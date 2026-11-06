## Set up symlinks
# ln -s ../mnist/PHM.py PHM.py
# ln -s ../mnist/PHMCluster.py PHMCluster.py
# ln -s ../mnist/distbn_helpers.py distbn_helpers.py
# ln -s ../mnist/model_fitting.py model_fitting.py

## File management
import dill

## Time/Memory Profiling
from time import time
import tracemalloc
import gc

## PHM from symlinks
from PHM import PHM
from PHMCluster import deltaPmc
from helpers import load_params 
from config import output_dir, density_fname, phm_fname, phm_time_fname

########################################
### Retrieve R parameters to run PHM ###
########################################
if not density_fname.exists():
    raise RuntimeError("Make sure to run density_estimation.R before this file!")

clusters, _ = load_params(str(density_fname))

#########################################
### Time and memory profiling for PHM ###
#########################################
print("==== Time and Memory Profiling ====")
print("\t# Components:", len(clusters))
if phm_time_fname.exists():
    with phm_fname.open("rb") as sf:
        phm = dill.load(sf)
    with phm_time_fname.open("rb") as sf:
        stats_dict = dill.load(sf)
else:
    print("\tComputing deltaPmc...")
    gc.collect()
    tracemalloc.start()
    start = time()
    phm = PHM(clusters, "alpha")
    phm.computeDeltaPmc(mcSamples=10**6, batchSize=10**5)
    elapsed_deltapmc = time() - start
    _, peak_deltapmc = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    print("\tRunning PHM...")
    gc.collect()
    tracemalloc.start()
    start = time()
    phm.runv2()
    elapsed_phm = time() - start
    _, peak_phm = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    with phm_fname.open("wb") as sf:
        dill.dump(phm, sf)

    stats_dict = {
        "time_deltapmc": elapsed_deltapmc, 
        "time_phm": elapsed_phm,
        "mem_deltapmc": peak_deltapmc,
        "mem_phm": peak_phm
    }
    with phm_time_fname.open("wb") as sf:
        dill.dump(stats_dict, sf)

## Time and Memory Results
base_str = "{name}\n\tTime: {time:.2f} seconds | Memory: {mem:.2f} MB"
print(base_str.format(
    name="DeltaPmc Calculation",
    time=stats_dict["time_deltapmc"],
    mem=stats_dict["mem_deltapmc"] / (1024**2)
))
print(base_str.format(
    name="PHM Algorithm",
    time=stats_dict["time_phm"],
    mem=stats_dict["mem_phm"] / (1024**2)
))