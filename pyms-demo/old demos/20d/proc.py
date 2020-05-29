"""proc.py
"""

# This file has been replaced by jupyter/comparing_datasets.ipynb

import pathlib
data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
# Change this if the data files are stored in a different location

from pyms.GCMS.IO.JCAMP import JCAMP_reader
from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.GCMS.Function import diff

# read the raw data
andi_file = data_directory / "gc01_0812_066.cdf"
jcamp_file = data_directory / "gc01_0812_066.jdx"

data1 = ANDI_reader(andi_file)
data2 = JCAMP_reader(jcamp_file)

# trim data2 between scans 1000 and 2000
# data2.trim(begin=1000,end=2000)

diff(data1, data2)

print("")

# trim data2 between scans 1000 and 2000
data2.trim(begin=1000,end=2000)

diff(data1, data2)
