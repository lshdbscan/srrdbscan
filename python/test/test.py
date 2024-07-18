import dbscan_srr as dbscan
import numpy as np
import sys
import h5py

f = h5py.File(sys.argv[1], "r")
eps = float(sys.argv[2])
minPts = int(sys.argv[3])
X = np.array(f['data'])

dbscan.SRR().fit_predict(X, 0.5, 3, True, "test", 56, -1, 1, eps, minPts)






