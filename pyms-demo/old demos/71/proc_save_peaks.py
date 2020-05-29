#################################################################
#       Creates Peak list for use with Display module
#################################################################

# This file has been replaced by jupyter/Displaying_Detected_Peaks.ipynb

from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.IntensityMatrix import build_intensity_matrix_i
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.TopHat import tophat

from pyms.Peak.List.IO import store_peaks

from pyms.BillerBiemann import (
	BillerBiemann,
	rel_threshold, num_ions_threshold,
	)

# read in raw data
andi_file = "gc01_0812_066.cdf"
data = ANDI_reader(andi_file)

data.trim("500s", "2000s")
# Build Intensity Matrix
im = build_intensity_matrix_i(data)
n_scan, n_mz = im.size

# perform necessary pre filtering
for ii in range(n_mz):
	ic = im.get_ic_at_index(ii)
	ic_smooth = savitzky_golay(ic)
	ic_bc = tophat(ic_smooth, struct="1.5m")
	im.set_ic_at_index(ii, ic_bc)

# Detect Peaks
peak_list = BillerBiemann(im, points=9, scans=2)

print("Number of peaks found: ", len(peak_list))

# Filter peaks
# Filter the peak list,
# first by removing all intensities in a peak less than a given relative
# threshold,
# then by removing all peaks that have less than a given number of ions above
# a given value

# Parameters
# percentage ratio of ion intensity to max ion intensity
r = 2

# minimum number of ions, n
n = 3
# greater than or equal to threshold, t
t = 10000

# trim by relative intensity
pl = rel_threshold(peak_list, r)

# trim by threshold
new_peak_list = num_ions_threshold(pl, n, t)

print("Number of filtered peaks: ", len(new_peak_list))

# store peak list
store_peaks(new_peak_list, 'output/peaks.bin')
