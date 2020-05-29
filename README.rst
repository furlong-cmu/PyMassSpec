************
PyMassSpec
************

.. start shields

.. image:: https://img.shields.io/travis/com/domdfcoding/PyMassSpec/master?logo=travis
    :target: https://travis-ci.com/domdfcoding/PyMassSpec
    :alt: Travis Build Status
.. image:: https://readthedocs.org/projects/pymassspec/badge/?version=latest
    :target: https://pymassspec.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
.. image:: https://img.shields.io/pypi/v/PyMassSpec.svg
    :target: https://pypi.org/project/PyMassSpec/
    :alt: PyPI
.. image:: https://img.shields.io/pypi/pyversions/PyMassSpec.svg
    :target: https://pypi.org/project/PyMassSpec/
    :alt: PyPI - Python Version
.. image:: https://img.shields.io/pypi/wheel/PyMassSpec
    :target: https://pypi.org/project/PyMassSpec/
    :alt: PyPI - Wheel
.. image:: https://img.shields.io/pypi/implementation/PyMassSpec
    :target: https://pypi.org/project/PyMassSpec/
    :alt: PyPI - Implementation
.. image:: https://img.shields.io/conda/v/domdfcoding/PyMassSpec
    :alt: Conda
    :target: https://anaconda.org/domdfcoding/PyMassSpec
.. image:: https://img.shields.io/conda/pn/domdfcoding/PyMassSpec?label=conda%7Cplatform
    :alt: Conda Platform
    :target: https://anaconda.org/domdfcoding/PyMassSpec
.. image:: https://coveralls.io/repos/github/domdfcoding/PyMassSpec/badge.svg?branch=master
    :target: https://coveralls.io/github/domdfcoding/PyMassSpec?branch=master
    :alt: Coverage
.. image:: https://img.shields.io/badge/license-GPLv2-orange
    :alt: License
    :target: https://github.com/domdfcoding/PyMassSpec/blob/master/LICENSE
.. image:: https://img.shields.io/github/languages/top/domdfcoding/PyMassSpec
    :alt: GitHub top language
.. image:: https://img.shields.io/github/commits-since/domdfcoding/PyMassSpec/v2.2.21
    :target: https://github.com/domdfcoding/PyMassSpec/pulse
    :alt: GitHub commits since tagged version
.. image:: https://img.shields.io/github/last-commit/domdfcoding/PyMassSpec
    :target: https://github.com/domdfcoding/PyMassSpec/commit/master
    :alt: GitHub last commit
.. image:: https://img.shields.io/maintenance/yes/2020
    :alt: Maintenance
.. image:: https://img.shields.io/codefactor/grade/github/domdfcoding/PyMassSpec
    :target: https://www.codefactor.io/repository/github/domdfcoding/PyMassSpec
    :alt: CodeFactor Grade

.. end shields

A Python toolkit for processing of chromatography--mass spectrometry data

PyMassSpec is a Python_ package for processing gas chromatography-mass spectrometry data.
PyMassSpec provides a framework and a set of components for rapid development and testing of methods for processing of chromatography--mass spectrometry data.
PyMassSpec can be used interactively through the Python shell, in a `Jupyter Notebook <https://jupyter.org/>`_, or the functions can be collected into scripts when it is preferable to perform data processing in the batch mode.

|

Forked from the original PyMS Repository: https://github.com/ma-bio21/pyms.
Originally by Andrew Isaac, Sean O'Callaghan and Vladimir Likić. The original publication can be found here: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-115

The original project seems to have been abandoned as there has been no activity since 2017.

|

.. contents:: Table of Contents
    :local:



The PyMassSpec project
=========================

The directory structure of PyMassSpec is as follows:

.. code-block:: text

    /
    ├── pyms:      The PyMassSpec code
    │
    ├── pyms-data: Example GC-MS data files
    │
    ├── pyms-demo: Examples of how to use PyMassSpec
    │
    ├── tests: pytest tests
    │
    └── UserGuide: Sphinx source for documentation

Features
=========

Installation
==============

PyMassSpec can be installed with the following command:

.. code-block:: bash

    $ pip --user install PyMassSpec

This will also install the following dependencies:

.. code-block:: bash

    numpy >= 1.16.2
    scipy >= 1.2.1
    pymzml >= 2.2.1
    matplotlib >= 3.0.2
    openpyxl >= 2.6.2
    netCDF4 >= 1.5.0
    biopython >= 1.74
    deprecation >= 2.0.6


PyMassSpec can also make use of 'mpi4py' if it is installed. See https://mpi4py.readthedocs.io/en/stable/ for further information.


Usage
=======

A tutorial illustrating various PyMassSpec features in detail is provided
in subsequent chapters of this User Guide. The commands executed
interactively are grouped together by example, and can be found
:ref:`here <pyms-demo>`.

.. If you are viewing this source, the examples can be found in the pyms-demo directory, and the data files in pyms-data

The data used in the PyMassSpec documentation and examples is available
:ref:`here <pyms-demo/data-files>`.

In the ":ref:`Demos and Examples <pyms-demo>`" section there
is a page corresponding to each example, coded with the chapter number
(ie. ":ref:`pyms-demo/20a/`" corresponds to the Example 20a, from Chapter 2).

Each example has a script named 'proc.py' which contains the commands given in the example.
These scripts can be run with the following command:

.. code-block:: bash

    $ python3 proc.py

Example processing GC-MS data
-------------------------------

Download the file ``gc01_0812_066.jdx`` and save it in the folder ``data``.
This file contains GC-MS data in the the JCAMP-DX format.

First the raw data is loaded:

    >>> from pyms.GCMS.IO.JCAMP import JCAMP_reader
    >>> jcamp_file = "data/gc01_0812_066.jdx"
    >>> data = JCAMP_reader(jcamp_file)
    -> Reading JCAMP file 'Data/gc01_0812_066.jdx'
    >>> data
    <pyms.GCMS.Class.GCMS_data at 0x7f3ec77da0b8>

The intensity matrix object is then built by binning the data:

    >>> from pyms.IntensityMatrix import build_intensity_matrix_i
    >>> im = build_intensity_matrix_i(data)

In this example, we show how to obtain the dimensions of the
newly created intensity matrix, then loop over all ion chromatograms,
and for each ion chromatogram apply Savitzky-Golay noise filter
and tophat baseline correction:

    >>> n_scan, n_mz = im.size
    >>> from pyms.Noise.SavitzkyGolay import savitzky_golay
    >>> from pyms.TopHat import tophat
    >>> for ii in range(n_mz):
    ...     print("working on IC", ii)
    ...     ic = im.get_ic_at_index(ii)
    ...     ic1 = savitzky_golay(ic)
    ...     ic_smooth = savitzky_golay(ic1)
    ...     ic_base = tophat(ic_smooth, struct="1.5m")
    ...     im.set_ic_at_index(ii, ic_base)

The resulting noise and baseline corrected ion chromatogram is saved back into the intensity matrix.

Further examples can be found in the `documentation`_

Contributing
==============

Contributions are very welcome. Tests can be run with `pytest`_. Please
ensure the coverage is at least .. image:: https://coveralls.io/repos/github/domdfcoding/pymassspec/badge.svg?branch=master
before you submit a pull request.

For further information see the section `Contributing to PyMassSpec`_

License
=========
PyMassSpec is Free and Open Source software released under the `GNU General Public License version 2 <GPL_>`__.


Issues
========

If you encounter any problems, please `file an issue`_ along with a
detailed description.


.. _`documentation`: https://pymassspec.readthedocs.io
.. _`Contributing to PyMassSpec`: https://pymassspec.readthedocs.io/en/master/Contributing/Contributing.html
.. _`pytest`: https://pytest.org/
.. _`file an issue`: https://github.com/domdfcoding/pymassspec/issues
.. _Python: https://www.python.org/
.. _GPL: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html