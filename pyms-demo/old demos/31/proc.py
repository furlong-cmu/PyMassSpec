"""proc.py
"""
# This file has been replaced by jupyter/IntensityMatrix.ipynb

import pathlib

data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
# Change this if the data files are stored in a different location

from pyms.GCMS.IO.JCAMP import JCAMP_reader
from pyms.IntensityMatrix import build_intensity_matrix

# from pyms.Utils.IO import save_data

# read the raw data
jcamp_file = data_directory / "gc01_0812_066.jdx"
data = JCAMP_reader(jcamp_file)

# IntensityMatrix
# must build intensity matrix before accessing any intensity matrix methods.

# default, float masses with interval (bin interval) of one from min mass
print("default intensity matrix, bin interval = 1, boundary +/- 0.5")
im = build_intensity_matrix(data)

#
# MassSpectrum
#

ms = im.get_ms_at_index(0)
print(ms)

# attributes and properties
print(f"Length of the Mass Spectrum: {len(ms)}")
print(f"Length of the Mass List: {len(ms.mass_list)}")
print(f"Length of the Intensity List: {len(ms.intensity_list)}")

#
# IonChromatogram
#

# TIC from raw data
tic = data.tic
print(tic)

# get the first ion chromatogram of the IntensityMatrix
ic0 = im.get_ic_at_index(0)
print(ic0)

# get the ion chromatogram for m/z = 73
ic73 = im.get_ic_at_mass(73)
print(ic73)

# some tests on ion chromatogram objects
print("'tic' is a TIC:", tic.is_tic())
print("'ic' is a TIC:", ic0.is_tic())
print("'ic' is a TIC:", ic73.is_tic())

# save TIC to a file
tic.write("output/tic.dat", minutes=True)

# save IC to a file
ic0.write("output/ic_index_0.dat", minutes=True)
ic73.write("output/ic_mass_73.dat", minutes=True)
