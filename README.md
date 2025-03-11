chelsa_cmip6
-----------
This package contains functions to creates monthly high-resolution 
climatologies for min-, max-, and mean temperature, precipitation rate 
and bioclimatic variables from anomalies and using CHELSA V2.1 as 
baseline high resolution climatology. Only works for GCMs for
hich tas, tasmax, tasmin, and pr are available. It is part of the
CHELSA Project: (CHELSA, <https://www.chelsa-climate.org/>).


# INSTALLATION
------------
chelsa-cmip6 can be installed on your machine in different ways. It is compatible with Windows, MacOS, and Linux systems.
## On Linux

Linux usually come already with a version of python preinstalled. 
#### From Github
If you want the latest version its best if you install it directly from github:
Example for Ubuntu 20.04
if you have not installed Git on your linux machine yet, the apt package management tool is the easiest way to install Git.

1. To update the packages, launch a terminal window, and enter:
```bash
sudo apt-get update
```

2. To install from the default repositories, enter the following:
```bash
sudo apt-get install git
```

3. you can now install the latest version from github by typing the following command in the terminal:
```bash
python3 -m pip install git+https://gitlabext.wsl.ch/karger/chelsa_cmip6.git
```

#### From PyPI:
If you want the latest release version, you can also install it from PyPI using:

using pip:
```bash
python3 -m pip install chelsa-cmip6
```
## On Windows
#### Install python:
1. Check if python is already installed on your machine. To do so open the command line interface by:
Press: Win + R

2. This will open a run window:
Type: cmd.exe and press enter
This should open the command line interface of Windows

3. In the command line interface type:
```bash
python3
```
If python is already installed, you will automatically enter the command line interpreter of python. If
Python is not already installed, a Microsoft Store window will open asking you to install Python. Click 'Get'
and install Python. Make sure your version is at least 3.8.

4. Verify your installation by typing the following in the command line interface:
```bash
python3 --version
```

#### Install the chelsa-cmip6 package:
##### Using pip
If you want the latest release version, you can also install it from PyPI using:

```bash
python3 -m pip install chelsa-cmip6
```

##### Using Git
If you want the latest version, you can also install it from Github. If you dont have Git
installed on your machine, you can find the installers here: 
https://git-scm.com/download/win

Follow the steps in the installer and after completion open a new command line interface in windows and type:
```bash
python3 -m pip install git+https://gitlabext.wsl.ch/karger/chelsa_cmip6.git
```

##### GDAL and rasterio in chelsa-cmip6 V.1.0 & V1.1
chelsa-cmip6 up to version  1.1. depends on rasterio, which in turn depends on GDAL. Sometimes these two modules can create problems with the 
installation routine. If the installation hangs while installing the rasterio module, you can try installing both modules manually
before installing chelsa-cmip6.

GDAL can be quite complex to build and install, particularly on Windows and MacOS
If you run into problems follow the steps provided for your respective system here:

You can get information how to install GDAL on your system here: https://pypi.org/project/GDAL/

After the manual installation of GDAL you can install rasterio by:
```bash
python3 -m pip install rasterio
```

starting V.1.2 chelsa-cmip6 does not rely on rasterio anymore, but uses byte wise connections to remote CHELSA
data stored in netcdf files.

# HOW TO USE
----------
The chelsa_cmip6 module provides functions to create monthly climatologies from climate
simulation data from CMIP6 using climate observation data from CHELSA V.2.1
at a 0.0083333° grid resolution for a given area of choice.

The GetClim module contains classes and functions to connect to CMIP6 data
via the Google cloud storage and to read the data into xarrays. It also creates
monthly climatologies using the delta change anomally correction method for a given 
time period. 

The BioClim module contains classes calculating various bioclimatic parameters
from climatological data (see: https://chelsa-climate.org/bioclim).

The delta change method that is applied is relatively insensitive towards individual model 
bias of a GCM, as it only uses the difference (ratio) for a given variable between
a reference period and a future period. In case of temperature an additive delta change 
is applied. In case of precipitation a multiplicative delta change is applied by 
adding a constant of 0.0000001 kg m^-2 s^⁻1 to both the reference and the future data
to avoid division by zero. 

The code only runs for CMIP6 models for which all needed variables tas, tasmax, tasmin, pr,
are available for both the reference and the future period.

The standard reference period is 1981-01-01 - 2010-12-31. If another reference period is 
chosen, the code conducts a delta change for this period as well. Best practice would be to 
choose the standard reference period.





EXAMPLES: 
------------
You can use the package in two ways 1. by importing the module in python, and 2. by using
the run_chelsa_cmip6.py wrapper function in the terminal (Linux, MAC) or command prompt (Windows).

Figure 1 gives an visual example of the different model parameters for Example 1. 


![explanation of the delta change method](/figs/Fig1-1.png)*Example of the delta change method applied on the model MPI-ESM1-2-LR for ssp585, the reference period 1981-2010, and the future period 2041-2070 (Example 1). The MPI-ESM1-2-LP gives the low resolution reference period temperature tas_low^ref  and the low resolution future period temperature tas_low^fut. The difference between these two temperatures is interpolated to 30 arcsec resolution using a qubic-spline CS(Δtas_fut^ref) and then added to the high resolution reference temperature tas_high^ref to get tas_high^fut. The parameters that need to be set to achieve the shown delta change downscaling are shown in the upper left corner. Both tas_high^ref and tas_high^fut are given as output for a specific region that is specified using the arguments: --xmin 5.3 --xmax 10.4 --ymin 46.0 --ymax 47.5 .*



EXAMPLE 1: 
------------

To create future climate data within python, first import the main function by:

## Using python
------------
```python
from chelsa_cmip6.GetClim import chelsa_cmip6
```

Creating long term climatological normals and the related bioclimatic variables that are commonly used in species distribution modeling 
is controlled via the fefps and fefpe parameters of the chelsa_cmip6 function. You can use function by running the following command 
in python a python prompt. Open a python prompt by either typing python 
in your terminal in Linux, or a command prompt in Windows.

You can then set the parameters of the chelsa_cmip6 function to create the climate data for the CMIP6 model you want.
If we want to create climatologies and bioclimatic variables for the model MPI-ESM1-2-LR and ssp585 for the years
2041-2070 for the region between 5.3° - 10.4° longitude, and 46.0° - 47.5° latitude and save them in your home 
directory (~/, on a linux system), we need to set the parameters of the function as follows: 
```python
chelsa_cmip6(activity_id='ScenarioMIP', 
             table_id='Amon', 
             experiment_id='ssp585', 
             institution_id='MPI-M', 
             source_id='MPI-ESM1-2-LR', 
             member_id='r1i1p1f1', 
             refps='1981-01-15', 
             refpe='2010-12-15', 
             fefps='2041-01-15', 
             fefpe='2070-12-15', 
             xmin=5.3, 
             xmax=10.4,
             ymin=46.0, 
             ymax=47.5,
             output='~/',
             use_esgf=False) 
```
#### on Windows:
If you are on a windows system the 'output' parameter should be in the form windows requires it, e.g.

```python
output='C:/Users/your_user_name/' # the directory you want the output to be saved in 
```

Important: Plese be aware that depending on the computer you use the function will only run if you have access to the internet, and it might take a considerable amount of time to finish as the data transfer needed is relativly large.


EXAMPLE 2: 
------------

Creating a monthly timeseries for the same model requires only an adaptation of the fefps, and fefpe parameter 
of the function. Here we show an example using a simple loop in python. The output will be a netCDF files for each 
month from 2016, 2100 for tas, tasmax, tasmin, and pr, and an annual timeseries for the bioclimatic variables.
Notice the change in dates in fefps and fefpe, which need to be the end and start of the year. If the model you choose uses a 360 day calender (e.g. UKESM1-0-LL), the last day of the year is the 30th. 

```python
for year in range(2016,2101):
    chelsa_cmip6(activity_id='ScenarioMIP', 
                 table_id='Amon', 
                 experiment_id='ssp585', 
                 institution_id='MPI-M', 
                 source_id='MPI-ESM1-2-LR', 
                 member_id='r1i1p1f1', 
                 refps='1981-01-15', 
                 refpe='2010-12-15', 
                 fefps=str(year) + '-01-01', 
                 fefpe=str(year) + '-12-31', 
                 xmin=5.3, 
                 xmax=10.4,
                 ymin=46.0, 
                 ymax=47.5,
                 output='~/')
```


## Use chelsa_cmip6 without using python directly
------------
You can also use the function directly from the terminal if you like. The chelsa_cmip6 package comes with a wrapper function,
run_chelsa_cmip6.py that allows you to simply use it via the terminal without any python knowledge required. It means however that you have to add directory where the wrapper function run_chelsa_cmip6.py is located to your PATH variable, so that your system knows where to look for it. 

#### on Linux
1. Find the directory in which chelsa_cmip6 is located by typing the following in the terminal:
```bash
python3 -c "import chelsa_cmip6; print(chelsa_cmip6.__file__);"
```
The printed path without the '\__init__.py' at the end is the path you will need to add to your PATH variable.

2. You can now add it to your PATH variable by using:
```bash
export PATH="your_path:$PATH"
```

3. Restart your terminal. You should now be able to use chelsa_cmip6 by running e.g. the following in the terminal:
```bash 
run_chelsa_cmip6.py --activity_id "ScenarioMIP" --table_id "Amon" \
--experiment_id "ssp585" --institution_id "MPI-M" --source_id "MPI-ESM1-2-LR" \
--member_id "r1i1p1f1" --refps "1981-01-15" --refpe "2010-12-15" \
--fefps "2041-01-15" --fefpe "2070-12-15" --xmin 5.3 --xmax 10.4 \
--ymin 46.0 --ymax 47.5 --output "~/"
```

#### on Windows
On Windows you need to first find out where the python package is installed. You can do so by typing the following in the command line interface:

```bash
python3 -c "import chelsa_cmip6; print(chelsa_cmip6.__file__);"
```

The printed path without the '\__init__.py' at the end is the path you will need to add to your path variable. On my machine it 
looks like this:

```bash 
C:\Users\Administrator\AppData\Local\Packages\PythonSoftwareFoundation.Python.3.10_qbz5n2kfra8p0\LocalCache\local-packages\Python310\site-packages\chelsa_cmip6\__init__.py
```

You need to add this path to your PATH environment variable. To do so you need to open a command prompt as Administrator. To do so, you need to:
1. Press: Win + R
2. Type: cmd.exe
3. Use Ctrl + Shift + Click/Tap on the OK button

4. Type the following in the command prompt, replacing C:\your\path\here with the your local path
```bash 
setx /M path "%path%;C:\your\path\here"
```
5. Close the command prompt and open a new one. You now should be able to run the the function now by e.g. typing:
```bash 
run_chelsa_cmip6.py --activity_id "ScenarioMIP" --table_id "Amon" \
--experiment_id "ssp585" --institution_id "MPI-M" --source_id "MPI-ESM1-2-LR" \
--member_id "r1i1p1f1" --refps "1981-01-15" --refpe "2010-12-15" \
--fefps "2041-01-15" --fefpe "2070-12-15" --xmin 5.3 --xmax 10.4 \
--ymin 46.0 --ymax 47.5 --output "C:/Users/<your_user_name>/"
```

```bash
run_chelsa_cmip6.py --activity_id "ScenarioMIP" --table_id "Amon" \
--experiment_id "ssp585" --institution_id "MPI-M" --source_id "MPI-ESM1-2-LR" \
--member_id "r1i1p1f1" --refps "1981-01-15" --refpe "2010-12-15" \
--fefps "2041-01-15" --fefpe "2070-12-15" --xmin 5.3 --xmax 10.4 \
--ymin 46.0 --ymax 47.5 --output "~/"
```

important is that the combination of activity_id 'ScenarioMIP' and e.g. experiment_id 'ssp585' is set to a combination that exists.
You can also get historical data but in that case, activity_ID, experiment_id, and fefps and fefps need to be changed. E.g. 

```bash
run_chelsa_cmip6.py --activity_id "CMIP" --table_id "Amon" \
--experiment_id "historical" --institution_id "MPI-M" \
--source_id "MPI-ESM1-2-LR" --member_id "r1i1p1f1" \
--refps "1981-01-15" --refpe "2010-12-15" \
--fefps "1851-01-15" --fefpe "1880-12-15" \
--xmin 5.3 --xmax 10.4 --ymin 46.0 --ymax 47.5 --output '~/'
```

it is also important that your fefps and fefpe are covered by the experiment_id and activity_id.

These reference periods are possible for example:

'ScenarioMIP' - 2016-01-01 - 2100-12-31

'CMIP' - 1850-01-01 - 2015-12-31

refps and refpe need to be in the range 1850-01-01 - 2015-12-31.

## chelsa_cmip6 in R via the reticulate package
------------
You can also use the chelsa_cmip6 package in R if you want to integrate it into your R workflow. To do so open an R console or R Studio and follow the steps below. Important: You still need to have python and the chelsa_cmip6 package installed (see instructions above).
Tested with R 4.2.1 in Windows 10

Install and load the reticulate package
```R
install.packages("reticulate",dependencies = TRUE)
library(reticulate)
```

Install and load the chelsa-cmip6 package
```R
py_install("chelsa-cmip6", pip=T)
```

The import() function enables you to import the chelsa_cmip6 module and call it’s functions directly from R.
```R
chelsa_cmip6 <- import('chelsa_cmip6')
```

You can then use the chelsa_cmip6 function in R the same you would in python. Be aware that the function might run for a while before it creates any output.
```R
chelsa_cmip6$GetClim$chelsa_cmip6(activity_id='ScenarioMIP', 
                                  table_id='Amon', 
                                  experiment_id='ssp585', 
                                  institution_id='MPI-M', 
                                  source_id='MPI-ESM1-2-LR', 
                                  member_id='r1i1p1f1', 
                                  refps='1981-01-15', 
                                  refpe='2010-12-15', 
                                  fefps='2041-01-15', 
                                  fefpe='2070-12-15', 
                                  xmin=5.3, 
                                  xmax=10.4,
                                  ymin=46.0, 
                                  ymax=47.5,
                                  output='~/')
```


# REQUIREMENTS
------------
chelsa_cmip6 is written in Python 3. It has been tested to run well with the
following Python release and package versions. The dependencies will be installed automatically.
- python 3.8.10
- xarray 0.16.2
- requests 2.25.1
- numpy 1.19.5
- rasterio 1.2.1 (not needed after version 1.1)
- pandas 1.1.5
- zarr 2.6.1
- gcsfs 0.7.2
- datetime 3.9.2
- scipy 0.19.1
- fsspec
- dask
- netcdf4
- h5netcdf
- pyesgf



## SINGULARITY
------------
All dependencies are also resolved in the singularity container '/singularity/chelsa_cmip6.sif'. Singularity needs to be installed on the respective linux system you are using. 
An installation guide can be found here: https://sylabs.io/guides/3.3/user-guide/quick_start.html#quick-installation-steps

The singularity container is available only here: https://gitlabext.wsl.ch/karger/chelsa_cmip6/-/blob/master/singularity/chelsa_cmip6.sif

If you use chelsa_cmip6 together with singularity the command should be slightly modified:

```bash
singularity exec /singularity/chelsa_cmip6.sif python3 chelsa_cmip6.py \
--activity_id 'CMIP' --table_id 'Amon' --experiment_id 'historical' \
--institution_id 'MPI-M' --source_id 'MPI-ESM1-2-LR' --member_id 'r1i1p1f1' \
--refps '1981-01-15' --refpe '2010-12-15' --fefps '1851-01-15' \
--fefpe '1880-12-15' --xmin 5.3 --xmax 10.4 --ymin 46.0 \
--ymax 47.5 --output '/home/karger/scratch/'
```

tested with singularity version 3.3.0-809.g78ec427cc
but newer versions usually work as well.


## CHECKING IF ALL NEEDED INPUT IS AVAILABLE
------------
Not all models and activities provide all the necessary input needed for chelsa_cmip6.py.
chelsa_cmip6.py will only work for GCMs that are both available for the historical period
and the respective scenario of interest. You can check this by using the CMIP6 data search
interface on e.g. https://esgf-node.llnl.gov/search/cmip6/ 
There you can filter for the different parameters (e.g. experiment_id) and see if a dataset
exists. E.g. by using the parameters given in the example. To check if also the historical
data exists for the model, just change the activity_id to 'CMIP' and the experiment_id to 'historical'.
Make sure the four variables needed do exist both for the scenario and the historical period:

These variables are needed:
- pr
- tas
- tasmax
- tasmin

Important: As chelsa_cmip6 is using by default the Pangeo archive, there might be GCMs which are not available. In 
that case the 'use_esgf' flag in che chelsa_cmip6() function can be set to 'True', to use the ESGF data archive instead.


## OUTPUT
------------
The output consist of netCDF4 files. There will be different files for each variable and seperatly for
the reference (refps - refpe) and the future period (fefps - fefpe). 
Additionally, there will be netCDF4 files for the 
different bioclimatic variables each for both the reference (refps - refpe) and the future period (fefps - fefpe). 


##  CAVEATS
------------
Important: 
1. Please be aware that depending on the computer you use the function will only run if you have access to the internet. 
2. It might take a considerable amount of time to finish as the data transfer needed is relativly large.
3. You will need a considerable amount of RAM depending on the geographical size you choose. The example given runs with 16 GB of RAM. Below that we cannot garantee that the function will succeed.


# COPYRIGHT
---------
(C) 2022 Dirk Nikolaus Karger


# LICENSE
-------
chelsa_cmip6 is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the
Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

chelsa_cmip6 is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with chelsa_cmip6. If not, see <http://www.gnu.org/licenses/>.


## CITATION:
------------
If you need a citation for the output, please refer to the article describing the high
resolution climatologies:

Karger, D.N., Conrad, O., Böhner, J., Kawohl, T., Kreft, H., Soria-Auza, R.W., Zimmermann, N.E., Linder, P., Kessler, M. (2017). Climatologies at high resolution for the Earth land surface areas. Scientific Data. 4 170122. https://doi.org/10.1038/sdata.2017.122


## CONTACT
-------
<dirk.karger@wsl.ch>


## AUTHOR
------
Dirk Nikolaus Karger
Swiss Federal Research Institute WSL
Zürcherstrasse 111
8903 Birmensdorf 
Switzerland
