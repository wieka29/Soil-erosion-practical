# Soil-erosion-practical
The aim of this practical is to have an idea of how coarse-resolution global models of soil loss and land surface processes can be used to study the soil organic carbon pool. We will make use of the CE-DYNAM model (Naipal et al., 2018, Biogeosciences), and we will only look at the carbon erosion rates triggered by soil removal on cropland. 
The Practical is accompanied by a general introduction into soil erosion and carbon modelling with land surface models.

The input data to execute the model is provided after contacting the author.

Have a look at the Jupyter notebook Erosion_and_carbon_training.ipynb to see how the code and results looks like, and the steps that need to be followed to run the model. 

# Installation of Python with Anaconda and necessary packages
1) Install Anaconda with python 2.7 (ideally): http://docs.continuum.io/anaconda/install/
2) Install fortran compilers to run f2py in Anaconda
For Windows:
I. Download and install MinGW-w64 through the anaconda prompt: conda install
mingw
II. Insert the path to the mingw bin folder to your system path. The path that
you want to add to your system variable should look something like C:\mingw\
mingw64\bin. Go to control panel > advanced system settings > environmental
variables > system variables > add new > variable is gfortran.exe and path is C:\
Users\adminuser\Anaconda2\MinGW\bin
III. In Anaconda create the environment in which you want f2py working, and
activate it:
conda create -n py36_test python=3.6
activate py36_test
VI. Configure your Anaconda environment to use MinGW when it compiles by
executing the command below from the Anaconda prompt: (echo [build] & echo compiler
= mingw32) > CONDA_PREFIX%\Lib\distutils\distutils.cfg
VIII. Install several necessary packages to your environment: conda install
numpy libpython m2w64-toolchain
For Linux :
Just make sure numpy is installed, which should be the case when you have installed
anaconda. And install gfortran: apt install gfortran
For Mac:
Just make sure numpy is installed, which should be the case when you have installed
anaconda. And install gcc: apt install gcc
3) Launch spyder from py27_test environment in which you want to work. Run the command: spyder
activate py27_test
You can also open the anaconda-navigator and change the environment from “base” to “py27_test”.
Then launch spyder
4) Launch spyder, then go to tools>preferences>Ipython console>graphics>change “backend” to
Automatic
5) Install Basemap: conda install -c conda-forge basemap
