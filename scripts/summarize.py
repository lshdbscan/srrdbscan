import numpy as np 
from collections import Counter
import h5py
import sys

fn = sys.argv[1]
f = h5py.File(fn, 'r')
arr = np.array(f['grp']['labels'])
f.close()

c = Counter(arr)

n = len(arr)
c = dict(c)
print(c)
max_size = max(c.values())
num_clusters = len(c)
num_border = c.get(-1, 0)
frac_largest = max_size / n

print(f'Summary for {fn}')
print(f"n = {n}, number of clusters = {num_clusters}, number of border points = {num_border}")
print(f"Largest cluster has {max_size} points, accounting for {100 * frac_largest}%")
