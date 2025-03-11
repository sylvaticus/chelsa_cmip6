from distutils.core import setup
setup(
  name = 'chelsa_cmip6',         
  packages = ['chelsa_cmip6'],  
  version = '1.3',
  license='GNU',
  description = "This package contains function to create monthly high-resolution climatologies for a selected geographic area for min-, max-, and mean temperature, precipitation rate and bioclimatic variables from anomalies and using CHELSA V2.1 as baseline high resolution climatology. Only works for GCMs for which tas, tasmax, tasmin, and pr are available.", 
  author = 'Dirk Nikolaus Karger',                  
  author_email = 'dirk.karger@wsl.ch',     
  url = 'https://gitlabext.wsl.ch/karger/chelsa_cmip6.git',
  download_url = 'https://gitlabext.wsl.ch/karger/chelsa_cmip6/-/archive/v1.0/chelsa_cmip6-v1.0.tar.gz', 
  keywords = ['CMIP6', 'climate', 'delta-change', 'CHELSA', 'bioclimate', 'growing degree days', 'gdd', 'bio'],  
  install_requires=[
          'numpy',
          'xarray',
          'pandas',
          'zarr',
          'gcsfs',
          'scipy',
          'fsspec',
          'dask',
          'netcdf4',
          'h5netcdf',
          'esgf-pyclient',
          'siphon'
      ],
  classifiers=[
    'Development Status :: 5 - Production/Stable',      
    'Intended Audience :: Science/Research',      
    'Topic :: Scientific/Engineering :: GIS'
  ],
)
