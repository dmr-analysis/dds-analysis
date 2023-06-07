
Installation
=================


.. contents::
    :local:


Download
--------

dds_analysis is written in Python. It can be installed and accessed from the command line and is available for both Linux and macOS operating systems. The package can be downloaded `here <https://github.com/dmr-analysis/dds-analysis/archive/refs/heads/master.zip>`_. Alternatively, you can use the following command:

.. code-block:: bash

   $ wget https://github.com/dmr-analysis/dds-analysis/archive/refs/heads/master.zip

Install
-------

It is highly recommended to create a separate virtual environment for the package to avoid any library conflicts. You can create a virtual environment using the following commands. We recommend installing and using Miniconda/Anaconda (`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_). You can find a tutorial on creating and updating virtual environments `here <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_.

If Miniconda is already installed, you can proceed with the following step-by-step installation. We have provided a quick installation setup file named ``quick_install.sh`` for your convenience. Running a simple bash command will automatically prepare the package for use:

.. code-block:: bash

   $ ./quick_install

If the quick_install.sh script is unsuccessful, you can follow the step-by-step details provided below:

1. Create a virtual environment:

.. code-block:: bash

   $ conda create --name dmr_env python==3.9.16

2. Activate the virtual environment:

.. code-block:: bash

   $ conda activate dmr_env

3. Install pip if not already installed:

.. code-block:: bash

   $ conda install pip


Please allow any other installations when prompted.

Requirements
------------

Before installing the package, the dependencies must be fulfilled. It is advised to install the dependencies using Miniconda. The list of dependencies is as follows:

- matplotlib==3.5.3
- numpy==1.21.5
- pandas==1.4.4
- scikit_learn==1.2.2
- scipy==1.9.1
- setuptools==65.6.3
- statsmodels==0.13.5
- bedtools==2.27.0

These dependencies can be installed one by one using the Conda package manager. For example:

.. code-block:: bash

   $ conda install numpy==1.21.5

A requirements.txt file is provided with the package. You can automatically install all the requirements using the following command:

.. code-block:: bash

   $ conda install --file requirements.txt

Alternatively, you can install the requirements using pip:

.. code-block:: bash

   $ pip install -r requirements.txt

Package Installation
--------------------

To install the package, navigate to the ``dds_analysis`` directory (the folder containing setup.py and pyproject.toml) and run the following command:

.. code-block:: bash

   $ pip install .

For more details, refer to the readme file in the package.

Package Contents
----------------

The package folder will contain the following:

- ``demo``: Contains function scripts.
- ``dds_analysis``: Contains the Python source code of the pipeline.
- ``readme.txt``: Instructions about the usage of the package.
- ``requirements.txt``: List of requirements that can be used for automatic installation using Miniconda or pip.
- ``setup.py``: Setup file for the package.
- ``project.toml``: Setup file for the package.
- ``data``: Contains input and output data for the secondary functions.
