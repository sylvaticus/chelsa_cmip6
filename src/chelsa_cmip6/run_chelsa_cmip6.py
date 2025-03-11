#!/usr/bin/env python

#This file is part of chelsa_cmip6.
#
#chelsa_cmip6 is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#chelsa_cmip6 is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with chelsa_cmip6.  If not, see <https://www.gnu.org/licenses/>.

from GetClim import chelsa_cmip6
import argparse

ap = argparse.ArgumentParser(
    description='''# This script creates monthly high-resolution 
    for min-, max-, and mean temperature, precipitation rate 
    and bioclimatic variables from anomalies and using CHELSA V2.1 as 
    baseline high resolution climatology. Only works for GCMs for
    which tas, tasmax, tasmin, and pr are available.
    ''',
    epilog='''author: Dirk N. Karger, dirk.karger@wsl.ch'''
)
ap.add_argument('-s', '--source_id', type=str,
                help="Source model (GCM), e.g. MPI-ESM1-2-LR, string")
ap.add_argument('-i', '--institution_id', type=str,
                help="Institution ID, e.g. MPI-M, string")
ap.add_argument('-t', '--table_id', type=str,
                help="table id, e.g. Amon, string")
ap.add_argument('-a', '--activity_id', type=str,
                help="activity id, e.g. ScenarioMIP, string")
ap.add_argument('-e', '--experiment_id', type=str,
                help="experiment id, e.g. ssp585, string")
ap.add_argument('-m', '--member_id', type=str,
                help="ensemble member, e.g. r1i1p1f1, string")
ap.add_argument('-rs', '--refps', type=str,
                help="reference period start, e.g. 1981-01-01, date format YYYY-MM-DD, string")
ap.add_argument('-re', '--refpe', type=str,
                help="reference period end, e.g. 2010-12-31, date format YYYY-MM-DD, string")
ap.add_argument('-fs', '--fefps', type=str,
                help="anomaly period start, e.g. 2041-01-01, date format YYYY-MM-DD, string")
ap.add_argument('-fe', '--fefpe', type=str,
                help="anomaly period end, e.g. 2070-01-01, date format YYYY-MM-DD, string")
ap.add_argument('-xn', '--xmin', type=float,
                help="minimum longitudinal extent of the boundary box, float")
ap.add_argument('-xm', '--xmax', type=float,
                help="maximum longitudinal extent of the boundary box, float")
ap.add_argument('-yn', '--ymin', type=float,
                help="minimum latitudinal extent of the boundary box, float")
ap.add_argument('-ym', '--ymax', type=float,
                help="maximum latitudinal extent of the boundary box, float")
ap.add_argument('-o', '--output', type=str,
                help="output directory, needs to exist, string")
ap.add_argument('-esfg', '--esgf', type=bool,
                help="use esgf instead of pangeo, default=False, bool")
ap.add_argument('-b', '--bio', type=bool,
                help="compute bio variables, default=True, bool")
                
                

args = ap.parse_args()
print("Downscaling:")
print(args)


if __name__ == '__main__':
    chelsa_cmip6(source_id=args.source_id,
                 institution_id=args.institution_id,
                 table_id=args.table_id,
                 activity_id=args.activity_id,
                 experiment_id=args.experiment_id,
                 member_id=args.member_id,
                 refps=args.refps,
                 refpe=args.refpe,
                 fefps=args.fefps,
                 fefpe=args.fefpe,
                 xmin=args.xmin,
                 xmax=args.xmax,
                 ymin=args.ymin,
                 ymax=args.ymax,
                 output=args.output,
                 use_esgf=args.esgf,
                 bio=args.bio)
