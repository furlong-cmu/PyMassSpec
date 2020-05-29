
'''
A file for loading Agilent GCMS data.  This a python rewrite of code from MATLAB code developed by James Dillon (https://github.com/chemplexity/chromatography)


@author P. Michael Furlong
@date 11 October 2019

'''
import os.path
import struct
import pathlib
from datetime import datetime

import numpy as np

# Project files
from pyms.GCMS.Class import GCMS_data
from pyms.Spectrum import Scan



def uint8(data):
    return struct.unpack('>B', data)[0]

def uint16(data):
    return struct.unpack('>H', data)[0]

def uint32(data):
    return struct.unpack('>L', data)[0]

def int16(data):
    return struct.unpack('>h', data)[0]

def int32(data):
    return struct.unpack('>l', data)[0]

def Agilent_reader(file_name):

    if not isinstance(file_name, (str, pathlib.Path)):
            raise TypeError("'file_name' must be a string or a pathlib.Path object")

    if not isinstance(file_name, pathlib.Path):
            file_name = pathlib.Path(file_name)


    assert os.path.exists(file_name)
    if (file_name / 'DATA.MS').exists():
        d_file = open(file_name / 'DATA.MS', 'rb')
    elif (file_name / 'data.ms').exists():
        d_file = open(file_name / 'data.ms', 'rb')
    else:
        print(f'Error: {file_name} does not contain a data.ms file.')
        raise ValueError(f'Error: {file_name} does not contain a data.ms file.')

    data = AgilentGCMSData()
    options = Options()

    load_file_info(d_file, data, options)
    load_tic(d_file, data, options)
    load_xic(d_file, data, options)

    time_list = list(data.time)
    rows, _ = data.xic.shape
    scan_list = [Scan(data.mz,list(data.xic[r,:])) for r in range(rows)]

    return GCMS_data(time_list, scan_list)

def load_file_info(d_file, data, options):

    # TODO: Make this smarter.  it's all over the shop.
    # Sample Name
    d_file.seek(24)
    data.sample.name = d_file.read(uint8(d_file.read(1))).strip()

    # Sample description
    d_file.seek(86)
    data.sample.description = d_file.read(uint8(d_file.read(1))).strip()

    d_file.seek(252)
    data.sample.sequence = int16(d_file.read(2))
    data.sample.vial = int16(d_file.read(2))
    data.sample.replicate = int16(d_file.read(2))

    # Method Name
    d_file.seek(228)
    data.method.name = d_file.read(uint8(d_file.read(1))).strip()

    # Method Operator
    d_file.seek(148)
    data.method.operator = d_file.read(uint8(d_file.read(1))).strip()

    # Method date/time
    d_file.seek(178)
    date_str = d_file.read(uint8(d_file.read(1))).strip().decode('UTF-8')
    print(date_str)
    try:
        data.method.date = datetime.strptime(date_str, '%d %b %y   %I:%M %p')
    except ValueError:
        try:
            data.method.date = datetime.strptime(date_str, '%d/%b/%y   %I:%M:%S %p')
        except ValueError:
            try:
                data.method.date = datetime.strptime(date_str, '%d-%b-%y,   %H:%M:%S')
            except ValueError:
                data.method.date = date_str
            ### end try
        ### end try
    ### end try

    # Instrument name
    d_file.seek(208)
    data.instrument.name = d_file.read(uint8(d_file.read(1))).strip()

    # Instrument inlet
    d_file.seek(218)
    data.instrument.inlet = d_file.read(uint8(d_file.read(1))).strip()

    # Total scans
    d_file.seek(278)
    options.scans = uint32(d_file.read(4))

    # TIC Offset
    d_file.seek(260)
    options.offset_tic = 2 * int32(d_file.read(4)) - 2
    if options.offset_tic < 0:
        raise RuntimeError(f'The file {d_file.name} does not meet Agilent format specification')
    ### end if

    # XIC Offset
    d_file.seek(options.offset_tic)
    options.offset_xic = []
    for _ in range(options.scans):
        options.offset_xic.append(2 * int32(d_file.read(4)) - 2)
        d_file.seek(d_file.tell() + 8) # Skipping 8 bytes after each read

#     options.offset_xic = [2 * int32(d_file.read(4)) - 2 for x in range(options.scans)]

    # Normalization Offset
    dat =  d_file.read(4)
    if len(dat) == 4:
        options.offset_normalization = 2 * int32(dat) - 2
### end load_file_info

def load_tic(d_file, data, options):

    # Variables
    scans = options.scans
    offset = options.offset_tic
    time_scale = 60000.

    data.time = []
    data.tic = []

    # Time values
    d_file.seek(offset+4)
    for _ in range(scans):
        data.time.append(int32(d_file.read(4))/time_scale)
        d_file.seek(d_file.tell() + 8) # Skip 8 bytes after the last read
    ### end for

    # Total Intensity Values
    d_file.seek(offset+8)
    for _ in range(scans):
        data.tic.append(int32(d_file.read(4)))
        d_file.seek(d_file.tell() + 8) # Skip 8 bytes after the last read
### end load_tic

def load_xic(d_file, data, options):

    assert options.scans is not None and options.offset_xic is not None

    scans = options.scans
    offset = options.offset_xic

    mz = []
    xic = []
    ns = []
    for i in range(scans):
        # Scan size
        d_file.seek(offset[i])
        n_i = (int16(d_file.read(2)) -18) // 2 + 2

        # Mass values
        d_file.seek(offset[i] + 18)
        for _ in range(n_i):
            mz_val = uint16(d_file.read(2))
            mz_val = round((mz_val / 20.)*10**(options.precision))/options.precision
            mz.append(mz_val)

            xic_val = uint16(d_file.read(2))
            xic.append((xic_val & 16383) * (8**(xic_val >> 14)))
        ### end for _ in range(n_i)
        ns.append(n_i)
    ### end for

    data.mz =  np.unique(mz)
    index = np.zeros((scans,2),dtype=int)

    index[:,1] = np.cumsum(ns,dtype=int)
    index[:,0] = np.roll(index[:,1],1)
    index[0,0] = 0

    xic = np.array(xic)
    data.xic = np.zeros((len(data.time), len(data.mz)))
    # Index columns
    (in_data, cols) = ismember(mz, data.mz)

    for i in range(scans):
        data.xic[i,cols[index[i,0]:index[i,1]+1]] = xic[index[i,0]:index[i,1]+1]
    ### end for
### end load_xic # Don't need to return data, options because of side-effects.

def ismember(xs, ys):
    lia_b = np.array([(x in ys, np.min(np.where(x == ys))) for x in xs])
    return zip(*lia_b)


class SampleInfo:
    def __init__(self):
        self.name = None
        self.description = None
        self.sequence = None
        self.vial = None
        self.replicate = None
    def __str__(self):
        return 'Sample:\n\t Name: %s\n\t Desc: %s\n\t Seq: %s\n\t Vial: %s\n\t Replicate: %s\n'%(self.name,self.description,self.sequence, self.vial, self.replicate)
### end class SampleInfo

class MethodInfo:
    def __init__(self):
        self.name = None
        self.operator = None
        self.date = None

    def __str__(self):
        return 'Method:\n\t Name: %s\n\t Operator: %s\n\t Date: %s'%(self.name,self.operator,self.date)
### end class MethodInfo

class InstrumentInfo:
    def __init__(self):
        self.name = None
        self.inlet = None

    def __str__(self):
        return 'Instrument:\n\t Name: %s\n\t Inlet: %s'%(self.name,self.inlet)
### end class InstrumentInfo

class Options:
    def __init__(self):
        self.scans = None
        self.offset_tic = None
        self.offset_xic = None
        self.offset_normalization = None
        self.precision = 3 # default value
    def __str__(self):
        return 'Options:\n\t Scans: %d\n\t Offset TIC: %d\n\t Offset XIC: %s\n\t Offset Normalization: %d' % (self.scans, self.offset_tic, str(self.offset_xic[:10]), self.offset_normalization)
### end class Options

class AgilentGCMSData:
    def __init__(self):
        self.file = None
        self.sample = SampleInfo()
        self.method = MethodInfo()
        self.instrument = InstrumentInfo()
        self.time = None
        self.tic = None
        self.xic = None
        self.mz = None
    ### end __init__

    def __str__(self):
        return str(self.sample) + '\n' + str(self.method) + '\n Time: %s \n TIC: %s\n XIC: %s' %(self.time[:10], self.tic[:10], self.xic[:10])

