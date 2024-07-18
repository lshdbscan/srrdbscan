# SRRDBSCAN

This is the source code for SRRDBSCAN. The code is available in `src`, an example program is provided in `apps`. A python wrapper is provided in `python`. 

The experimental evaluation uses the Python wrapper. The evaluation framework can be found at <https://github.com/cikm24-1183/evaluation-framework>.

# Requirements

On a standard ubuntu 22.04, 

```
apt update -y && apt install -y libtbb-dev python3-numpy python3-scipy python3-pip build-essential git
```

covers all requirements to build the code.


How to compile:
```bash
    mkdir build; cd build; cmake ..; make
```

Building the python package works as follows

```bash 
   python setup.py build
```

Input files are supposed to be HDF5 files.
All data processing is provided with the evaluation framework <https://github.com/cikm24-1183/evaluation-framework>.


# Acknowledgements

The code is based on: 

> Amir Keramatian, Vincenzo Gulisano, Marina Papatriantafilou, Philippas Tsigas: IP.LSH.DBSCAN: Integrated Parallel Density-Based Clustering Through Locality-Sensitive Hashing. Euro-Par 2022: 268-284.

