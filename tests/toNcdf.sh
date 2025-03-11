for var in 'hurs' 'clt' 'sfcWind'
do
for month in {1..12}
do
infile='/storage/brunp/Data/CHELSA_V.2.1/'${var}'/Climatology/CHELSA_'${var}'_'$(printf "%02d" ${month})'_1981-2010_V.2.1.tif'
outfile='/storage/karger/chelsa_V2/scaled/CHELSA_'${var}'_'$(printf "%02d" ${month})'_1981-2010_V.2.1.nc'
tempfile='/home/karger/scratch/tmp.nc'
gdal_translate ${infile} ${tempfile}
ncks -4 -L 6 --cnk_dmn lat,500 --cnk_dmn lon,500 ${tempfile} /home/karger/scratch/out.nc
cp /home/karger/scratch/out.nc ${outfile}
rm /home/karger/scratch/out.nc
rm ${tempfile}
done
done