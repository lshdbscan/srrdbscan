# brings the results in a somewhat more readable format
import os 
import sys
import glob
import shutil

label = "Labels_"

source="build/"
target = "results/"

files = list(glob.glob(os.path.join(source, label + "*")))

new_files = 0

for fn in files:
    original_fn = fn

    fn = fn.split("/")[-1]
    fn = fn[:-3]
    log_file = "Benchmark_" + "_".join(fn.split("_")[1:]) + ".txt"

    #print(fn)
    values = fn.split("_")[2::2]
    #print(values)
    dataset = values.pop(0)[:-5]
    eps, minPts = float(values[0]), int(values[1])
    delta, memLimit = float(values[2]), float(values[3])
    level, shrinkage = int(values[4]), float(values[5])
    result_dir = os.path.join(target, dataset, str(eps), str(minPts))
    os.makedirs(result_dir, exist_ok=True)
    new_fn = f"delta={delta}_memLimit={memLimit}_level={level}_shrinkage={shrinkage}"
    shutil.move(original_fn, os.path.join(result_dir, f"{new_fn}.hdf5"))
    new_files += 1
    if os.path.exists(os.path.join(source, log_file)):
        shutil.move(os.path.join(source, log_file), os.path.join(result_dir, f"log_{new_fn}.txt"))

print(f"Organized result files. Moved {new_files} files.")



