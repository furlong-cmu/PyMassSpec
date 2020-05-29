*********************
Tests
*********************

|pkgname| uses `pytest`_ for unit testing. Tests can be run with the following command:

.. code-block:: bash

    $ pytest tests/

pytest requires the following additional dependencies:

.. literalinclude:: ../../tests/requirements.txt

These can be installed with the following command:

.. code-block:: bash

    $ python3 -m pip install -r tests/requirements.txt


A series of reference images for ``test_Display.py`` are in the "tests/baseline" directory. If these files need to be regenerated, run the following command:

.. code-block:: bash

    $ pytest --mpl-generate-path="tests/baseline" tests/test_Display.py


|

Please ensure that the coverage percentage is at least the current value of ...