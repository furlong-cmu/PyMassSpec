#############################################################################
#                                                                           #
#    PyMassSpec software for processing of mass-spectrometry data           #
#    Copyright (C) 2019 Dominic Davis-Foster                                #
#                                                                           #
#    This program is free software; you can redistribute it and/or modify   #
#    it under the terms of the GNU General Public License version 2 as      #
#    published by the Free Software Foundation.                             #
#                                                                           #
#    This program is distributed in the hope that it will be useful,        #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of         #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
#    GNU General Public License for more details.                           #
#                                                                           #
#    You should have received a copy of the GNU General Public License      #
#    along with this program; if not, write to the Free Software            #
#    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.              #
#                                                                           #
#############################################################################

import csv
import pickle
from copy import deepcopy

import pytest

from tests.constants import *

from pyms.GCMS.Class import GCMS_data
from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.GCMS.Function import diff, ic_window_points
from pyms.IonChromatogram import IonChromatogram
from pyms.Spectrum import Scan
from pyms.TopHat import tophat
from pyms.BillerBiemann import BillerBiemann, rel_threshold, num_ions_threshold
from pyms.IntensityMatrix import build_intensity_matrix, build_intensity_matrix_i
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.Peak.Function import peak_sum_area, peak_top_ion_areas
from pyms.Peak.Class import Peak
from pyms.Experiment import Experiment


@pytest.fixture(scope="module")
def andi(datadir):
	print("data")
	return ANDI_reader(datadir / "gc01_0812_066.cdf")


"""
@pytest.fixture(scope="module")
def im_andi(data):
	# build an intensity matrix object from the data
	return build_intensity_matrix(data)


@pytest.fixture(scope="module")
def tic_andi(data):
	# get the TIC
	return deepcopy(data.tic)


@pytest.fixture(scope="module")
def im_i_andi(data):
	# build an intensity matrix object from the data
	return build_intensity_matrix_i(data)


@pytest.fixture(scope="function")
def peak_list(im_i):
	im_i = deepcopy(im_i)
	
	# Intensity matrix size (scans, masses)
	n_scan, n_mz = im_i.size
	
	# noise filter and baseline correct
	for ii in range(n_mz):
		ic = im_i.get_ic_at_index(ii)
		ic_smooth = savitzky_golay(ic)
		ic_bc = tophat(ic_smooth, struct="1.5m")
		im_i.set_ic_at_index(ii, ic_bc)
	
	# Use Biller and Biemann technique to find apexing ions at a scan
	# default is maxima over three scans and not to combine with any neighbouring
	# scan.
	peak_list = BillerBiemann(im_i, points=9, scans=2)
	return peak_list


@pytest.fixture(scope="function")
def filtered_peak_list(im_i, peak_list):
	# peak_list = deepcopy(peak_list)
	# do peak detection on pre-trimmed data
	# trim by relative intensity
	apl = rel_threshold(peak_list, 2, copy_peaks=False)
	
	# trim by threshold
	new_peak_list = num_ions_threshold(apl, 3, 30, copy_peaks=False)
	
	# ignore TMS ions and set mass range
	for peak in new_peak_list:
		peak.crop_mass(50, 400)
		peak.null_mass(73)
		peak.null_mass(147)
		
		# find area
		area = peak_sum_area(im_i, peak)
		peak.area = area
		area_dict = peak_top_ion_areas(im_i, peak)
		peak.ion_areas = area_dict
	
	return new_peak_list


@pytest.fixture(scope="session")
def peak(im_i):
	scan_i = im_i.get_index_at_time(31.17 * 60.0)
	ms = im_i.get_ms_at_index(scan_i)
	return Peak(12.34, ms)


@pytest.fixture(scope="session")
def ms(im_i):
	return deepcopy(im_i.get_ms_at_index(0))


@pytest.fixture(scope="session")
def scan(data):
	# return deepcopy(im_i.get_scan_at_index(0))
	return deepcopy(data.scan_list[0])


@pytest.fixture(scope="function")
def expr(filtered_peak_list):
	# create an experiment
	return Experiment("ELEY_1_SUBTRACT", filtered_peak_list)"""


def test_ANDI_reader(datadir):
	# Errors
	for type in [test_float, test_int, test_list_strs, test_dict, test_list_ints, test_tuple]:
		with pytest.raises(TypeError):
			ANDI_reader(type)
	
	with pytest.raises(FileNotFoundError):
		ANDI_reader(test_string)


# def test_ANDI_OpenChrom_reader(datadir):
# todo


def test_GCMS_data(andi):
	assert isinstance(andi, GCMS_data)
	
	GCMS_data(andi.time_list, andi.scan_list)
	
	# Errors
	for type in [test_string, test_float, test_int, test_list_strs, test_dict]:
		with pytest.raises(TypeError):
			GCMS_data(type, andi.scan_list)
	
	for type in [test_string, test_float, test_int, test_list_ints, test_list_strs, test_dict, test_tuple]:
		with pytest.raises(TypeError):
			GCMS_data(andi.time_list, type)


def test_len(andi):
	assert len(andi) == 9865


def test_equality(andi):
	assert andi == GCMS_data(andi.time_list, andi.scan_list)
	assert andi != GCMS_data(list(range(len(andi.scan_list))), andi.scan_list)
	assert andi != test_string
	assert andi != test_int
	assert andi != test_float
	assert andi != test_list_ints
	assert andi != test_list_strs
	assert andi != test_tuple
	assert andi != test_dict


def test_get_scan_list(andi):
	with pytest.warns(DeprecationWarning):
		andi.get_scan_list()


def test_get_tic(andi):
	with pytest.warns(DeprecationWarning):
		andi.get_tic()


def test_info(capsys, andi):
	andi.info()
	captured = capsys.readouterr()
	assert captured.out == """ Data retention time range: 5.093 min -- 66.795 min
 Time step: 0.375 s (std=0.000 s)
 Number of scans: 9865
 Minimum m/z measured: 50.000
 Maximum m/z measured: 599.900
 Mean number of m/z values per scan: 56
 Median number of m/z values per scan: 40
"""


def test_scan_list(andi):
	# raw scans
	scans = andi.scan_list
	
	assert isinstance(scans, list)
	assert isinstance(scans[0], Scan)
	assert len(scans[0]) == 622
	assert isinstance(scans[0].mass_list, list)
	# 1st mass value for 1st scan
	assert isinstance(scans[0].mass_list[0], float)
	assert scans[0].mass_list[0] == 50.099998474121094
	
	assert isinstance(scans[0].intensity_list, list)
	# 1st intensity value for 1st scan
	assert isinstance(scans[0].intensity_list[0], float)
	assert scans[0].intensity_list[0] == 22128.0
	
	# minimum mass found in 1st scan
	assert isinstance(scans[0].min_mass, float)
	assert scans[0].min_mass == 50.099998474121094
	
	# maximum mass found in 1st scan
	assert isinstance(scans[0].max_mass, float)
	assert scans[0].max_mass == 599.4000244140625


def test_tic(andi):
	tic = andi.tic
	assert isinstance(tic, IonChromatogram)
	# number of scans in TIC
	assert len(tic) == 9865
	assert len(tic) == len(andi.time_list)
	
	# start time of TIC
	assert isinstance(tic.get_time_at_index(0), float)
	assert tic.get_time_at_index(0) == 305.582
	assert isinstance(tic.get_index_at_time(305.6), int)
	assert tic.get_index_at_time(305.6) == 0
	assert tic.get_index_at_time(306) == 1
	
	assert isinstance(tic.get_intensity_at_index(44), float)
	assert tic.get_intensity_at_index(44) == 21685482.0
	
	assert isinstance(tic.time_list, list)
	assert tic.time_list[0] == 305.582
	
	assert isinstance(tic.time_step, float)
	
	assert isinstance(tic.is_tic(), bool)
	assert tic.is_tic()


def test_trim(andi):
	# time
	trimmed = deepcopy(andi)
	trimmed.trim("6.5m", "21m")
	
	assert trimmed.min_mass == 50.099998474121094
	assert trimmed.max_mass == 542.0
	
	time = trimmed.time_list
	assert len(time) == 2319
	assert time[0] == 390.404
	assert trimmed.get_index_at_time(400.0) == 26
	
	scans = trimmed.scan_list
	assert scans[0].mass_list[0] == 50.099998474121094
	
	# Scans
	trimmed = deepcopy(andi)
	trimmed.trim(10, 2000)
	
	assert trimmed.min_mass == 50.0
	assert trimmed.max_mass == 599.9000244140625
	
	time = trimmed.time_list
	assert len(time) == 1992
	assert time[0] == 308.96000000000004
	assert trimmed.get_index_at_time(1000.0) == 1841
	
	scans = trimmed.scan_list
	assert scans[0].mass_list[0] == 50.099998474121094
	
	trimmed.trim(end=1000)
	assert len(trimmed.time_list) == 1001
	
	trimmed.trim(begin=2)
	assert len(trimmed.time_list) == 1000
	
	# Errors
	with pytest.raises(SyntaxError):
		trimmed.trim()
	
	for type in [test_list_strs, test_dict, test_list_ints, test_tuple]:
		with pytest.raises(TypeError):
			trimmed.trim(begin=type)
		with pytest.raises(TypeError):
			trimmed.trim(end=type)


def test_write(andi, outputdir):
	andi.write(outputdir / "andi_gcms_data")
	
	# Errors
	for type in [test_list_strs, test_dict, test_list_ints, test_tuple, test_int, test_float]:
		with pytest.raises(TypeError):
			andi.write(type)
	
	# Read .I.csv and check values
	assert (outputdir / "andi_gcms_data.I.csv").exists()
	i_csv = list(csv.reader((outputdir / "andi_gcms_data.I.csv").open()))
	assert [float(x) for x in i_csv[5]] == [21472.0, 10271.0, 30952.0, 26112.0, 65952.0, 53528.0, 72520.0, 108064.0,
											151104.0, 21344.0, 22112.0, 25720.0, 89352.0, 9893.0, 12390.0, 8951.0,
											14378.0, 107008.0, 277120.0, 8388096.0, 974656.0, 423872.0, 1135616.0,
											98400.0, 44240.0, 30224.0, 4480.0, 42424.0, 29952.0, 72744.0, 24256.0,
											10265.0, 8244.0, 9253.0, 12053.0, 3416.0, 2243.0, 1179.0, 2202.0, 685696.0,
											77024.0, 33320.0, 18384.0, 27776.0, 15649.0, 4037.0, 3719.0, 1789.0, 2078.0,
											2872.0, 1392.0, 9969.0, 6838.0, 20824.0, 12440.0, 20264.0, 7201.0, 7010.0,
											4439.0, 7556.0, 18160.0, 4085.0, 2164.0, 1296.0, 1668.0, 1732.0, 2246.0,
											63072.0, 16036.0, 35048.0, 11699.0, 16896.0, 3640.0, 2643.0, 1496.0, 5687.0,
											1627.0, 1495.0, 4803.0, 3933.0, 229504.0, 41920.0, 24968.0, 14807.0,
											31864.0, 5795.0, 3501.0, 1087.0, 2425.0, 2223.0, 1279.0, 4605.0, 1956.0,
											2815.0, 1720.0, 3909.0, 2960.0, 1639.0, 1893.0, 55824.0, 7833.0, 3914.0,
											18240.0, 5855.0, 4115.0, 55400.0, 11178.0, 5628.0, 1226.0, 1066.0, 1543.0,
											3053056.0, 603456.0, 300480.0, 2014208.0, 365504.0, 162560.0, 21664.0,
											4818.0, 949.0, 533.0, 560.0, 609.0, 485.0, 496.0, 473.0, 395.0, 411.0,
											542.0, 1601.0, 827.0, 965.0, 539.0, 576.0, 564.0, 427.0, 388.0, 441.0,
											965.0, 422.0, 414.0, 544.0, 795.0, 819.0, 20704.0, 5872.0, 3512.0, 1290.0,
											1094.0, 911.0, 852.0, 1260.0, 687.0, 418.0, 379.0, 587.0, 559.0, 725.0,
											1395.0, 220032.0, 49920.0, 29592.0, 381952.0, 75400.0, 32552.0, 4710.0,
											1351.0, 571.0, 320.0, 452.0, 681.0, 613.0, 740.0, 888.0, 1366.0, 192256.0,
											46840.0, 19768.0, 4863.0, 1148.0, 1260.0, 566.0, 517.0, 472.0, 411.0, 385.0,
											411.0, 515.0, 499.0, 505.0, 369.0, 359.0, 355.0, 356.0, 509.0, 316.0, 352.0,
											405.0, 315.0, 357.0, 383.0, 332.0, 458.0, 344.0, 260.0, 328.0, 400.0, 437.0,
											306.0, 429.0, 426.0, 362.0, 296.0, 439.0, 407.0, 465.0, 537.0, 451.0, 380.0,
											387.0, 376.0, 437.0, 302.0, 358.0, 336.0, 411.0, 430.0, 339.0, 395.0, 434.0,
											339.0, 535.0, 305.0, 472.0, 406.0, 383.0, 494.0, 551.0, 303.0, 486.0, 362.0,
											402.0, 356.0, 363.0, 291.0, 452.0, 361.0, 648.0, 381.0, 484.0, 504.0, 508.0,
											445.0, 405.0, 368.0, 432.0, 329.0, 287.0, 354.0, 321.0, 345.0, 381.0, 412.0,
											410.0, 440.0, 498.0, 444.0, 420.0, 627.0, 491.0, 471.0, 417.0, 583.0, 445.0,
											417.0, 402.0, 319.0, 229.0, 404.0, 394.0, 355.0, 601.0, 313.0, 356.0, 403.0,
											372.0, 304.0, 291.0, 313.0, 493.0, 443.0, 312.0, 397.0, 447.0, 390.0, 466.0,
											285.0, 391.0, 402.0, 402.0, 387.0, 363.0, 347.0, 469.0, 375.0, 297.0, 440.0,
											362.0, 278.0, 454.0, 408.0, 379.0, 311.0, 436.0, 447.0, 224.0, 388.0, 391.0,
											426.0, 446.0, 401.0, 376.0, 409.0, 319.0, 288.0, 340.0, 326.0, 371.0, 246.0,
											281.0, 326.0, 319.0, 483.0, 433.0, 291.0, 344.0, 306.0, 406.0, 395.0, 310.0,
											388.0, 436.0, 414.0, 344.0, 325.0, 397.0, 398.0, 265.0, 421.0, 494.0, 483.0,
											441.0, 374.0, 316.0, 431.0, 431.0, 356.0, 398.0, 370.0, 466.0, 502.0, 467.0,
											453.0, 349.0, 514.0, 379.0, 376.0, 563.0, 475.0, 342.0, 393.0, 370.0, 472.0,
											476.0, 532.0, 234.0, 376.0, 421.0, 525.0, 388.0, 439.0, 431.0, 464.0, 433.0,
											394.0, 545.0, 494.0, 412.0, 391.0, 387.0, 414.0, 379.0, 457.0, 519.0, 337.0,
											348.0, 317.0, 439.0, 313.0, 324.0, 434.0, 440.0, 399.0, 277.0, 293.0, 391.0,
											376.0, 585.0, 223.0, 405.0, 440.0, 449.0, 396.0, 354.0, 474.0, 423.0, 425.0,
											329.0, 359.0, 517.0, 449.0, 503.0, 329.0, 401.0, 375.0, 511.0, 462.0, 460.0,
											506.0, 336.0, 334.0, 425.0, 465.0, 370.0, 372.0, 377.0, 444.0, 315.0, 360.0,
											462.0, 514.0, 467.0, 344.0, 468.0, 433.0, 516.0, 361.0, 368.0, 277.0, 361.0,
											456.0, 390.0, 391.0, 365.0, 513.0, 437.0, 376.0, 472.0, 452.0, 492.0, 399.0,
											385.0, 428.0, 302.0, 422.0, 389.0, 369.0, 361.0, 616.0, 451.0, 439.0, 427.0,
											439.0, 546.0, 376.0, 366.0, 421.0, 561.0, 363.0, 488.0, 391.0, 429.0, 274.0,
											345.0, 367.0, 469.0, 293.0, 291.0, 395.0, 498.0, 391.0, 368.0, 333.0, 358.0,
											450.0, 316.0, 439.0, 479.0, 551.0, 445.0, 449.0, 409.0, 407.0, 434.0, 365.0,
											611.0, 446.0, 328.0, 458.0, 643.0, 405.0, 410.0, 405.0, 536.0, 290.0, 300.0,
											307.0, 286.0, 430.0, 390.0, 327.0, 446.0, 439.0, 316.0, 395.0, 422.0, 650.0,
											324.0, 451.0, 362.0, 597.0, 400.0, 618.0, 377.0, 421.0, 414.0, 407.0, 327.0,
											390.0, 311.0, 411.0, 441.0, 468.0, 388.0, 425.0, 323.0, 557.0, 517.0, 442.0,
											432.0, 431.0, 355.0, 441.0, 434.0, 348.0, 310.0, 276.0, 401.0, 430.0, 485.0,
											478.0, 418.0, 371.0, 443.0, 391.0, 343.0, 419.0, 379.0, 433.0, 450.0, 411.0,
											587.0, 475.0, 509.0, 585.0, 581.0, 469.0, 411.0, 399.0, 447.0, 438.0, 462.0,
											446.0, 360.0, 381.0, 417.0, 431.0, 501.0, 310.0, 462.0, 611.0, 353.0, 338.0,
											396.0, 503.0, 330.0, 327.0, 522.0, 313.0, 371.0, 518.0, 539.0, 509.0, 519.0,
											480.0, 423.0, 382.0, 472.0, 410.0, 546.0, 585.0, 666.0, 483.0, 503.0, 518.0,
											381.0, 355.0, 428.0, 338.0, 558.0, 452.0, 379.0, 453.0, 554.0, 415.0, 365.0,
											586.0, 524.0, 428.0, 653.0, 341.0]
	assert [float(x) for x in i_csv[50]] == [22336.0, 10586.0, 30328.0, 27512.0, 62848.0, 52512.0, 70736.0, 109832.0,
											 144896.0, 19272.0, 21608.0, 28632.0, 88808.0, 9246.0, 13705.0, 7454.0,
											 15300.0, 105096.0, 278144.0, 8388096.0, 898816.0, 388672.0, 1112576.0,
											 98208.0, 43512.0, 28296.0, 5082.0, 42040.0, 31184.0, 70240.0, 23464.0,
											 9813.0, 7986.0, 9778.0, 11460.0, 3467.0, 2691.0, 1218.0, 2170.0, 3192.0,
											 679680.0, 75504.0, 31200.0, 18848.0, 24904.0, 13627.0, 3475.0, 4193.0,
											 1919.0, 1460.0, 1911.0, 2769.0, 1394.0, 8898.0, 6042.0, 18224.0, 11870.0,
											 19120.0, 6370.0, 6647.0, 4147.0, 6576.0, 17864.0, 3568.0, 1668.0, 1354.0,
											 1513.0, 1441.0, 2042.0, 60624.0, 14000.0, 30224.0, 10249.0, 14813.0,
											 3215.0, 3265.0, 1540.0, 6064.0, 1496.0, 1246.0, 5125.0, 3310.0, 207552.0,
											 39640.0, 23168.0, 14739.0, 29552.0, 5542.0, 3638.0, 1301.0, 2074.0, 2280.0,
											 1149.0, 4321.0, 1818.0, 2188.0, 1273.0, 3351.0, 2519.0, 1432.0, 2240.0,
											 68272.0, 9274.0, 5247.0, 21008.0, 6773.0, 4355.0, 54544.0, 11492.0, 6183.0,
											 1503.0, 1543.0, 1381.0, 2587.0, 2929152.0, 579904.0, 282816.0, 1942528.0,
											 360704.0, 157952.0, 19120.0, 4547.0, 1219.0, 446.0, 522.0, 677.0, 727.0,
											 391.0, 599.0, 326.0, 486.0, 751.0, 1703.0, 893.0, 660.0, 689.0, 690.0,
											 685.0, 427.0, 449.0, 480.0, 486.0, 579.0, 341.0, 840.0, 790.0, 576.0,
											 647.0, 577.0, 531.0, 666.0, 19872.0, 5350.0, 2884.0, 1143.0, 788.0, 619.0,
											 1282.0, 1111.0, 513.0, 500.0, 663.0, 1030.0, 1046.0, 919.0, 1765.0,
											 282432.0, 64224.0, 35376.0, 404352.0, 80072.0, 33416.0, 4932.0, 1442.0,
											 523.0, 621.0, 477.0, 501.0, 680.0, 624.0, 1078.0, 1903.0, 179136.0,
											 45688.0, 16712.0, 4468.0, 1160.0, 481.0, 406.0, 499.0, 344.0, 414.0, 306.0,
											 473.0, 300.0, 478.0, 424.0, 372.0, 489.0, 304.0, 371.0, 407.0, 271.0,
											 485.0, 406.0, 489.0, 270.0, 492.0, 412.0, 498.0, 325.0, 464.0, 682.0,
											 332.0, 390.0, 340.0, 493.0, 401.0, 628.0, 434.0, 390.0, 425.0, 478.0,
											 390.0, 392.0, 465.0, 377.0, 393.0, 450.0, 437.0, 321.0, 304.0, 400.0,
											 434.0, 306.0, 387.0, 440.0, 259.0, 346.0, 437.0, 426.0, 521.0, 451.0,
											 488.0, 476.0, 365.0, 394.0, 362.0, 402.0, 343.0, 331.0, 372.0, 395.0,
											 434.0, 464.0, 340.0, 411.0, 515.0, 446.0, 529.0, 339.0, 437.0, 434.0,
											 537.0, 539.0, 497.0, 277.0, 474.0, 473.0, 433.0, 350.0, 306.0, 453.0,
											 312.0, 384.0, 457.0, 523.0, 544.0, 477.0, 427.0, 419.0, 431.0, 388.0,
											 562.0, 396.0, 539.0, 444.0, 382.0, 414.0, 386.0, 404.0, 356.0, 398.0,
											 392.0, 376.0, 560.0, 357.0, 308.0, 429.0, 415.0, 436.0, 409.0, 277.0,
											 271.0, 346.0, 380.0, 355.0, 310.0, 292.0, 379.0, 369.0, 416.0, 409.0,
											 406.0, 421.0, 413.0, 396.0, 496.0, 283.0, 501.0, 302.0, 289.0, 365.0,
											 453.0, 405.0, 504.0, 410.0, 384.0, 440.0, 424.0, 416.0, 404.0, 461.0,
											 412.0, 323.0, 356.0, 420.0, 464.0, 336.0, 472.0, 399.0, 363.0, 362.0,
											 497.0, 424.0, 450.0, 357.0, 370.0, 393.0, 251.0, 302.0, 327.0, 379.0,
											 368.0, 233.0, 365.0, 398.0, 285.0, 416.0, 339.0, 355.0, 485.0, 396.0,
											 398.0, 308.0, 467.0, 434.0, 348.0, 354.0, 342.0, 233.0, 399.0, 503.0,
											 474.0, 548.0, 318.0, 344.0, 427.0, 278.0, 317.0, 377.0, 511.0, 314.0,
											 344.0, 398.0, 391.0, 363.0, 568.0, 365.0, 264.0, 364.0, 438.0, 478.0,
											 448.0, 435.0, 415.0, 318.0, 430.0, 351.0, 498.0, 279.0, 430.0, 474.0,
											 498.0, 334.0, 434.0, 546.0, 379.0, 413.0, 375.0, 584.0, 497.0, 448.0,
											 466.0, 424.0, 480.0, 447.0, 331.0, 467.0, 266.0, 450.0, 415.0, 298.0,
											 459.0, 369.0, 353.0, 455.0, 516.0, 401.0, 429.0, 434.0, 429.0, 318.0,
											 349.0, 397.0, 406.0, 457.0, 640.0, 391.0, 474.0, 282.0, 389.0, 325.0,
											 300.0, 387.0, 409.0, 513.0, 432.0, 312.0, 485.0, 387.0, 363.0, 454.0,
											 287.0, 437.0, 401.0, 486.0, 428.0, 359.0, 355.0, 438.0, 431.0, 337.0,
											 509.0, 421.0, 522.0, 550.0, 498.0, 456.0, 446.0, 537.0, 460.0, 413.0,
											 448.0, 461.0, 541.0, 299.0, 482.0, 437.0, 455.0, 389.0, 379.0, 429.0,
											 440.0, 454.0, 478.0, 578.0, 410.0, 571.0, 424.0, 360.0, 485.0, 428.0,
											 392.0, 427.0, 442.0, 268.0, 318.0, 438.0, 397.0, 525.0, 300.0, 423.0,
											 393.0, 435.0, 439.0, 451.0, 419.0, 406.0, 379.0, 623.0, 426.0, 449.0,
											 473.0, 511.0, 406.0, 548.0, 344.0, 502.0, 286.0, 528.0, 523.0, 491.0,
											 482.0, 442.0, 483.0, 473.0, 363.0, 597.0, 483.0, 360.0, 570.0, 669.0,
											 308.0, 389.0, 384.0, 558.0, 284.0, 406.0, 501.0, 489.0, 453.0, 360.0,
											 448.0, 568.0, 466.0, 268.0, 449.0, 469.0, 353.0, 459.0, 280.0, 502.0,
											 355.0, 463.0, 468.0, 547.0, 378.0, 256.0, 412.0, 567.0, 450.0, 464.0,
											 494.0, 406.0, 442.0, 417.0, 392.0, 362.0, 416.0, 431.0, 299.0, 413.0,
											 393.0, 556.0, 527.0, 428.0, 459.0, 534.0, 348.0, 357.0, 504.0, 561.0,
											 304.0, 433.0, 560.0, 455.0, 539.0, 347.0, 420.0, 491.0, 318.0, 343.0,
											 387.0, 486.0, 409.0, 395.0, 498.0, 454.0, 492.0, 465.0, 360.0, 413.0,
											 508.0, 558.0, 448.0, 571.0, 399.0, 356.0, 584.0, 578.0, 467.0, 340.0,
											 471.0, 342.0, 419.0, 411.0, 419.0, 377.0, 600.0, 400.0, 307.0, 436.0,
											 427.0, 390.0, 486.0, 434.0, 402.0, 332.0, 342.0, 614.0, 527.0, 532.0]
	assert [float(x) for x in i_csv[500]] == [1269.0, 1649.0, 3297.0, 1152.0, 590.0, 4132.0, 5335.0, 22784.0, 1826.0,
											  568.0, 281.0, 214.0, 230.0, 540.0, 399.0, 2524.0, 4016.0, 11251.0, 1296.0,
											  1614.0, 277.0, 5897.0, 547.0, 3550.0, 1006.0, 5995.0, 487.0, 417.0, 238.0,
											  874.0, 3127.0, 8184.0, 618.0, 181.0, 325.0, 241.0, 224.0, 168.0, 295.0,
											  442.0, 759.0, 1846.0, 1790.0, 398.0, 207.0, 179.0, 248.0, 741.0, 1034.0,
											  178.0, 640.0, 361.0, 432.0, 269.0, 274.0, 300.0, 156.0, 264.0, 343.0,
											  359.0, 879.0, 174.0, 183.0, 2012.0, 232.0, 2191.0, 363.0, 359.0, 180.0,
											  190.0, 882.0, 176.0, 190.0, 923.0, 223.0, 456.0, 324.0, 270.0, 164.0,
											  166.0, 690.0, 198.0, 379.0, 167.0, 363.0, 162.0, 619.0, 220.0, 6843.0,
											  1986.0, 1314.0, 272.0]
	
	# Read .mz.csv and check values
	assert (outputdir / "andi_gcms_data.mz.csv").exists()
	i_csv = list(csv.reader((outputdir / "andi_gcms_data.mz.csv").open()))
	assert [float(x) for x in i_csv[5]] == [50.1, 51.1, 53.1, 54.1, 55.1, 56.1, 57.2, 58.2, 59.1, 60.2, 61.1,
											62.20, 63.10, 64.10, 65.10, 66.20, 67.20, 69.20, 70.10, 73.20, 74.10, 75.20,
											77.10, 78.10, 79.10, 81.10, 82.20, 84.10, 85.10, 86.10, 87.10, 88.10, 89.10,
											90.10, 91.10, 92.10, 93.10, 93.80, 95.30, 100.10, 101.0, 102.10, 103.0,
											104.10, 105.10, 105.90, 106.90, 108.10, 109.0, 110.10, 111.10, 113.20,
											113.90, 115.10, 116.10, 117.10, 118.10, 119.10, 120.10, 121.10, 122.10,
											123.0, 124.10, 125.10, 126.20, 127.10, 128.30, 131.10, 132.10, 133.10,
											134.10, 135.10, 136.10, 136.90, 137.90, 139.10, 140.0, 141.10, 143.10,
											144.20, 147.10, 148.10, 149.0, 150.10, 151.10, 151.90, 153.10, 154.20,
											155.0, 155.90, 156.80, 158.0, 159.20, 160.10, 161.20, 162.0, 163.10, 164.20,
											165.10, 168.10, 169.10, 170.0, 172.10, 173.10, 174.20, 176.10, 177.10,
											178.10, 179.0, 179.90, 181.10, 188.10, 189.10, 190.10, 192.10, 193.10,
											194.10, 195.0, 196.0, 197.10, 197.80, 198.30, 198.80, 199.30, 200.10,
											200.90, 201.70, 202.30, 203.0, 205.20, 206.20, 206.80, 208.0, 209.0, 209.80,
											211.10, 212.10, 213.0, 214.10, 215.0, 215.60, 217.10, 218.30, 218.90, 221.0,
											222.10, 223.0, 224.0, 225.10, 225.80, 226.60, 227.10, 228.0, 228.90, 229.40,
											230.0, 230.70, 232.20, 233.90, 238.10, 239.10, 240.20, 242.10, 243.0,
											244.10, 245.0, 245.90, 247.0, 248.10, 248.70, 249.70, 250.70, 251.60,
											252.30, 253.20, 257.10, 258.10, 259.10, 260.0, 261.0, 261.30, 261.90,
											262.40, 263.20, 264.20, 265.0, 265.40, 266.30, 267.10, 268.0, 268.60,
											269.50, 270.60, 271.10, 272.10, 272.90, 273.50, 274.50, 275.40, 275.80,
											276.30, 276.80, 277.30, 277.70, 278.40, 279.20, 279.70, 280.70, 281.70,
											282.90, 283.20, 283.60, 284.10, 285.10, 285.90, 286.50, 287.0, 288.0,
											288.60, 289.30, 289.90, 290.50, 291.20, 292.0, 292.60, 293.10, 293.80,
											294.80, 296.0, 297.20, 297.90, 299.10, 300.10, 301.10, 302.20, 303.20,
											304.0, 304.80, 305.30, 306.10, 307.0, 307.50, 308.20, 308.90, 309.50,
											310.70, 311.20, 311.70, 313.10, 314.0, 315.50, 316.20, 317.10, 317.90,
											318.70, 319.60, 320.30, 320.70, 321.50, 322.10, 322.60, 323.40, 324.20,
											325.10, 326.0, 327.10, 328.0, 329.10, 330.10, 331.40, 332.50, 333.70,
											334.30, 335.0, 335.80, 336.30, 336.90, 337.50, 338.30, 338.70, 339.60,
											340.30, 341.20, 341.80, 342.40, 342.70, 343.20, 343.70, 344.10, 344.70,
											345.30, 346.30, 347.10, 348.10, 348.60, 350.10, 350.70, 351.30, 351.90,
											352.40, 352.90, 353.50, 354.60, 355.10, 356.30, 356.80, 357.30, 358.50,
											359.40, 360.40, 361.30, 362.50, 363.20, 364.10, 364.40, 365.20, 365.90,
											366.80, 367.40, 368.70, 369.30, 370.70, 371.90, 372.50, 372.90, 373.50,
											374.10, 374.70, 375.40, 375.90, 376.20, 376.90, 377.50, 378.20, 379.0,
											379.60, 380.20, 380.90, 381.70, 382.80, 383.40, 384.10, 384.90, 385.50,
											386.0, 386.70, 387.60, 388.20, 388.80, 390.10, 390.90, 391.50, 392.60,
											393.20, 394.30, 394.80, 395.90, 396.70, 398.20, 398.80, 399.40, 399.90,
											400.70, 401.50, 402.30, 403.0, 403.40, 404.0, 404.60, 405.20, 406.10,
											406.40, 407.30, 408.80, 409.80, 410.60, 411.50, 413.30, 414.40, 415.50,
											415.90, 416.40, 417.0, 417.90, 418.80, 419.70, 420.50, 420.90, 421.80,
											422.50, 424.60, 425.30, 426.50, 427.40, 428.0, 428.50, 429.40, 430.20,
											430.70, 431.20, 432.20, 433.10, 434.20, 434.70, 435.30, 436.10, 436.50,
											437.50, 438.40, 439.10, 439.60, 440.0, 440.50, 441.80, 442.80, 443.80,
											444.40, 445.20, 445.90, 446.60, 447.20, 447.80, 449.0, 449.50, 450.30,
											451.20, 452.10, 452.70, 453.80, 454.50, 455.0, 455.80, 456.30, 457.50,
											458.20, 459.10, 459.90, 460.70, 461.20, 462.10, 462.70, 463.60, 464.30,
											465.40, 466.10, 467.10, 467.70, 468.10, 468.80, 469.60, 470.40, 470.90,
											471.40, 472.10, 473.10, 473.80, 474.30, 475.0, 475.50, 476.20, 477.40,
											478.40, 479.30, 480.10, 480.70, 481.20, 482.10, 482.70, 483.70, 484.20,
											485.20, 485.80, 486.60, 487.20, 488.0, 488.60, 489.30, 490.0, 490.90,
											491.50, 492.30, 492.70, 493.30, 493.90, 494.60, 495.40, 496.50, 497.10,
											497.90, 498.70, 499.70, 500.0, 500.60, 501.10, 501.70, 502.70, 503.30,
											504.20, 505.30, 505.90, 506.80, 507.60, 508.70, 509.30, 510.10, 511.10,
											512.0, 512.70, 513.20, 513.80, 514.70, 515.30, 515.90, 516.50, 517.50,
											518.0, 518.50, 519.0, 519.70, 520.50, 521.10, 521.60, 522.10, 523.10,
											523.70, 524.60, 525.30, 526.20, 526.80, 527.80, 528.30, 529.10, 530.0,
											530.90, 531.30, 532.0, 532.60, 533.40, 533.80, 534.90, 535.60, 535.90,
											536.80, 537.40, 537.90, 538.50, 539.10, 539.90, 540.60, 541.80, 542.0,
											543.0, 543.40, 543.90, 544.30, 544.90, 545.40, 546.0, 546.90, 547.70,
											548.40, 549.20, 549.80, 551.30, 551.70, 552.60, 553.10, 553.70, 554.40,
											554.80, 555.60, 556.20, 556.70, 558.0, 558.80, 559.50, 560.40, 561.0,
											561.80, 562.70, 563.20, 563.60, 564.40, 565.10, 565.70, 566.20, 566.90,
											567.30, 568.20, 568.80, 569.30, 569.90, 570.40, 571.0, 572.0, 572.70,
											573.20, 574.0, 575.30, 576.0, 577.0, 578.40, 578.90, 579.30, 580.20, 580.70,
											581.20, 581.80, 582.30, 583.30, 584.30, 584.90, 585.60, 586.10, 586.70,
											587.50, 588.10, 588.60, 589.20, 590.10, 590.70, 591.20, 591.90, 592.70,
											593.60, 594.20, 595.50, 596.30, 597.50, 598.90, 599.8000]
	assert [float(x) for x in i_csv[50]] == [50.10, 51.10, 53.10, 54.20, 55.10, 56.10, 57.20, 58.20, 59.10, 60.10,
											 61.20, 62.20, 63.10, 64.10, 65.20, 66.20, 67.20, 69.20, 70.10, 73.20,
											 74.10, 75.20, 77.10, 78.10, 79.10, 81.10, 82.20, 84.10, 85.20, 86.10,
											 87.10, 88.10, 89.10, 90.10, 91.10, 92.10, 93.0, 94.30, 95.20, 96.0, 100.10,
											 101.10, 102.10, 103.10, 104.10, 105.10, 106.10, 107.0, 107.90, 108.70,
											 109.20, 110.0, 111.20, 113.10, 113.80, 115.10, 116.10, 117.10, 118.0,
											 119.10, 120.10, 121.10, 122.10, 123.0, 124.10, 124.90, 126.20, 127.0,
											 128.20, 131.10, 132.20, 133.10, 134.10, 135.0, 136.0, 137.0, 138.10,
											 139.10, 140.0, 141.20, 143.10, 144.10, 147.10, 148.10, 149.10, 150.10,
											 151.0, 152.10, 153.0, 154.0, 155.0, 156.0, 156.90, 158.10, 159.10, 160.10,
											 160.90, 162.10, 163.10, 163.90, 165.10, 168.10, 169.10, 170.10, 172.10,
											 173.10, 174.20, 176.10, 177.10, 178.0, 178.90, 180.10, 180.70, 181.90,
											 188.10, 189.10, 190.10, 192.10, 193.10, 194.10, 194.90, 195.90, 197.10,
											 197.80, 198.30, 198.80, 199.40, 200.70, 201.50, 202.30, 202.90, 204.0,
											 205.10, 206.0, 207.0, 207.80, 208.0, 209.0, 209.70, 210.20, 210.70, 211.40,
											 212.10, 212.90, 213.90, 214.30, 214.90, 215.50, 216.0, 216.60, 217.20,
											 221.10, 222.0, 223.10, 224.0, 225.20, 225.70, 227.10, 228.0, 229.30,
											 230.10, 231.20, 232.50, 233.0, 233.50, 234.10, 238.10, 239.10, 240.20,
											 242.10, 243.10, 244.0, 244.90, 245.80, 246.90, 247.90, 248.80, 249.40,
											 250.80, 251.30, 252.0, 253.90, 257.10, 258.0, 259.10, 260.10, 260.90,
											 262.80, 263.50, 264.20, 264.90, 265.40, 266.30, 266.90, 267.60, 268.30,
											 269.80, 270.50, 271.10, 271.90, 272.40, 273.70, 274.40, 275.10, 275.90,
											 277.10, 278.0, 278.80, 279.90, 280.70, 281.20, 282.10, 283.20, 283.90,
											 284.50, 285.10, 285.60, 286.40, 287.10, 288.0, 288.50, 289.20, 290.0,
											 290.80, 291.30, 292.40, 293.30, 293.60, 294.30, 294.70, 295.70, 296.20,
											 296.80, 297.50, 298.30, 298.90, 299.80, 300.50, 301.40, 303.0, 303.60,
											 304.20, 305.40, 306.90, 307.50, 308.30, 309.10, 310.20, 310.60, 311.30,
											 311.70, 312.30, 313.80, 314.30, 315.30, 315.90, 316.40, 316.90, 317.50,
											 318.10, 318.70, 319.20, 319.80, 320.80, 321.60, 322.20, 322.90, 323.60,
											 323.80, 324.70, 325.30, 326.0, 327.0, 327.70, 328.50, 329.0, 329.80,
											 330.30, 331.90, 333.10, 333.70, 334.30, 334.80, 335.30, 336.60, 337.30,
											 337.90, 338.80, 339.50, 339.90, 340.50, 340.90, 341.70, 342.60, 343.40,
											 344.60, 345.50, 346.20, 347.10, 347.60, 348.60, 349.20, 350.0, 350.70,
											 351.50, 352.40, 353.40, 354.0, 354.60, 355.20, 355.60, 356.60, 357.20,
											 357.50, 358.90, 359.80, 360.20, 360.80, 362.10, 363.50, 364.40, 364.90,
											 365.30, 365.90, 366.90, 367.70, 369.40, 370.50, 370.90, 371.40, 372.0,
											 372.70, 373.60, 374.80, 376.0, 376.70, 377.40, 379.40, 380.20, 380.90,
											 381.60, 382.90, 383.20, 384.10, 384.90, 385.50, 386.30, 386.80, 387.20,
											 387.80, 388.30, 388.80, 389.70, 390.20, 390.80, 391.30, 392.30, 393.0,
											 393.50, 394.0, 394.60, 395.10, 395.90, 396.60, 397.30, 398.0, 398.50,
											 399.20, 399.70, 400.20, 400.80, 401.90, 403.10, 404.50, 405.90, 406.60,
											 407.60, 408.10, 408.80, 409.20, 409.90, 411.10, 411.60, 412.0, 412.90,
											 414.10, 414.70, 415.20, 415.70, 416.70, 417.20, 418.30, 419.10, 420.0,
											 420.60, 421.10, 422.0, 422.50, 423.40, 423.90, 424.80, 426.30, 426.80,
											 427.60, 428.70, 429.50, 430.10, 430.80, 431.40, 431.90, 432.70, 433.50,
											 434.40, 435.60, 436.20, 437.10, 437.50, 438.60, 439.30, 440.50, 441.50,
											 442.70, 443.30, 443.70, 444.60, 445.10, 445.70, 446.0, 446.80, 448.30,
											 448.80, 449.30, 450.10, 450.70, 451.40, 452.30, 453.20, 454.30, 455.0,
											 455.90, 456.70, 457.10, 457.70, 458.30, 458.80, 459.50, 460.40, 460.90,
											 461.70, 462.30, 463.10, 463.80, 465.0, 465.70, 466.50, 467.60, 468.60,
											 469.20, 469.80, 470.70, 471.80, 472.60, 473.10, 473.60, 474.70, 476.10,
											 476.80, 477.40, 478.40, 478.90, 479.70, 480.80, 481.40, 481.90, 482.30,
											 482.80, 483.30, 483.80, 484.70, 485.60, 486.70, 487.20, 488.10, 488.70,
											 489.20, 490.20, 491.10, 491.60, 492.40, 493.0, 493.70, 494.60, 495.80,
											 496.40, 497.90, 498.50, 499.40, 499.90, 500.40, 500.90, 501.90, 502.30,
											 502.90, 503.50, 504.30, 504.70, 505.20, 506.30, 506.70, 507.20, 508.20,
											 509.10, 509.90, 510.40, 511.0, 511.50, 512.40, 513.20, 513.90, 514.50,
											 515.30, 516.80, 517.60, 518.0, 518.50, 519.20, 519.60, 520.10, 520.80,
											 521.50, 522.70, 523.30, 524.80, 526.10, 526.80, 527.60, 528.20, 529.20,
											 530.20, 530.80, 531.50, 532.20, 532.80, 533.80, 534.70, 535.80, 536.60,
											 537.30, 538.0, 538.50, 539.50, 540.50, 541.20, 541.80, 542.70, 543.20,
											 543.80, 544.40, 544.90, 545.50, 546.50, 547.20, 547.50, 548.40, 548.90,
											 549.60, 550.0, 551.20, 551.70, 552.30, 553.10, 553.70, 554.50, 555.40,
											 556.10, 556.90, 557.80, 558.40, 559.30, 560.10, 561.0, 561.50, 562.0,
											 562.80, 563.70, 564.20, 565.0, 565.70, 566.20, 566.90, 567.80, 568.50,
											 568.90, 569.60, 570.30, 571.0, 571.40, 572.0, 572.90, 574.40, 575.10,
											 575.80, 576.50, 577.20, 577.80, 578.30, 578.90, 579.80, 580.30, 581.40,
											 582.30, 584.20, 585.10, 585.70, 586.60, 587.50, 588.30, 589.0, 589.80,
											 590.80, 591.70, 592.20, 592.70, 593.30, 594.30, 595.30, 595.90, 596.40,
											 597.0, 597.60, 598.40, 599.30, 599.9000]
	assert [float(x) for x in i_csv[500]] == [50.10, 51.10, 52.0, 53.10, 54.10, 55.0, 56.10, 57.10, 58.10, 59.0, 60.90,
											  63.0, 66.20, 67.10, 68.10, 69.0, 70.10, 71.10, 72.10, 73.0, 73.90, 75.10,
											  76.0, 77.10, 78.10, 79.10, 80.0, 81.10, 82.10, 83.10, 84.10, 85.10, 86.10,
											  86.80, 89.20, 89.90, 91.0, 92.10, 94.90, 96.20, 97.10, 98.10, 99.20,
											  100.0, 101.0, 103.10, 111.10, 112.0, 113.10, 114.0, 115.0, 116.0, 117.0,
											  117.90, 119.0, 119.80, 121.0, 125.0, 130.10, 130.90, 133.10, 134.0,
											  134.90, 142.20, 143.0, 147.10, 148.10, 148.90, 150.90, 154.10, 161.0,
											  161.90, 170.10, 177.10, 184.0, 188.0, 190.90, 192.90, 203.80, 205.0,
											  207.10, 208.0, 220.10, 220.90, 248.90, 250.0, 265.0, 266.90, 281.10,
											  282.0, 283.10, 284.0]


def test_write_intensities_stream(andi, outputdir):
	andi.write_intensities_stream(outputdir / "andi_intensity_stream.csv")
	
	# Errors
	for type in [test_list_strs, test_dict, test_list_ints, test_tuple, test_int, test_float]:
		with pytest.raises(TypeError):
			andi.write_intensities_stream(type)
	
	# Read and check values
	assert (outputdir / "andi_intensity_stream.csv").exists()
	intensity_stream = list((outputdir / "andi_intensity_stream.csv").open().readlines())
	assert intensity_stream[5] == "55416.0000\n"
	assert intensity_stream[50] == "7478.0000\n"
	assert intensity_stream[500] == "480.0000\n"


# Inherited Methods from pymsBaseClass

def test_dump(andi, outputdir):
	andi.dump(outputdir / "ANDI_dump.dat")
	
	# Errors
	for type in [test_list_strs, test_dict, test_list_ints, test_tuple, test_int, test_float]:
		with pytest.raises(TypeError):
			andi.dump(type)
	
	# Read and check values
	assert (outputdir / "ANDI_dump.dat").exists()
	loaded_data = pickle.load((outputdir / "ANDI_dump.dat").open("rb"))
	assert loaded_data == andi
	assert len(loaded_data) == len(andi)


# Inherited Methods from TimeListMixin

def test_time_list(andi):
	time = andi.time_list
	assert isinstance(time, list)
	# number of retention times
	assert len(time) == 9865
	# retention time of 1st scan:
	assert isinstance(time[0], float)
	assert time[0] == 305.582


def test_get_time_list(andi):
	with pytest.warns(DeprecationWarning):
		andi.get_time_list()


# Inherited Methods from MaxMinMassMixin

def test_get_max_mass(andi):
	with pytest.warns(DeprecationWarning):
		andi.get_max_mass()


def test_get_min_mass(andi):
	with pytest.warns(DeprecationWarning):
		andi.get_min_mass()


def test_max_mass(andi):
	# maximum mass found in all data
	assert isinstance(andi.max_mass, float)
	assert andi.max_mass == 599.9000244140625


def test_min_mass(andi):
	assert isinstance(andi.min_mass, float)
	# minimum mass found in all data
	assert andi.min_mass == 50.0


# Inherited Methods from GetIndexTimeMixin

def test_get_index_at_time(andi):
	# index of 400sec in time_list
	assert isinstance(andi.get_index_at_time(400.0), int)
	assert andi.get_index_at_time(400.0) == 252
	
	# Errors
	for type in [test_dict, test_list_ints, test_list_strs, test_string, test_tuple]:
		with pytest.raises(TypeError):
			andi.get_index_at_time(type)
	with pytest.raises(IndexError):
		andi.get_index_at_time(0)
	with pytest.raises(IndexError):
		andi.get_index_at_time(100000)


def test_get_time_at_index(andi):
	assert isinstance(andi.get_time_at_index(400), float)
	assert andi.get_time_at_index(400) == 455.71
	
	# Errors
	for type in [test_dict, test_list_ints, test_list_strs, test_string, test_tuple]:
		with pytest.raises(TypeError):
			andi.get_time_at_index(type)
	with pytest.raises(IndexError):
		andi.get_time_at_index(-1)
	with pytest.raises(IndexError):
		andi.get_time_at_index(1000000)

# Test GCMS.Function

# def test_diff(data):
# TODO


# def test_ic_window_points(data):
# todo
