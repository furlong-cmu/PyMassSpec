"""proc.py
"""
# This file has been replaced by jupyter/IntensityMatrix.ipynb

import pathlib
data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
# Change this if the data files are stored in a different location

from pyms.GCMS.IO.JCAMP import JCAMP_reader
from pyms.IntensityMatrix import build_intensity_matrix_i

# read the raw data
jcamp_file = data_directory / "gc01_0812_066.jdx"
data = JCAMP_reader(jcamp_file)

# IntensityMatrix
# must build intensity matrix before accessing any intensity matrix methods.

# integer intensity matrix, integer masses, in one unit steps
print("intensity matrix with integer mass and bin interval = 1, "
      "using default boundary -0.3, +0.7")
im = build_intensity_matrix_i(data)
print(im)

print(" -> size of intensity matrix (#scans, #bins):", im.size)

print(im.mass_list[:10])

print(" -> start mass:", im.min_mass)
print(" -> end mass:", im.max_mass)

index = im.get_index_of_mass(73.3)
print(" -> the index of the nearest mass to 73.3 m/z is:", index)
print(" -> the nearest mass to 73.3 m/z is:", im.get_mass_at_index(index))
