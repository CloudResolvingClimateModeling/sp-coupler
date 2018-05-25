#!/usr/bin/env python

# Superparameteriation coupling code for OpenIFS <--> Dales
#
# Fredrik Jansson, Gijs van den Oord 
# 2017-2018
# 

from __future__ import print_function
import argparse
import logging

import os
import shapely.geometry
import sys
import json

from splib import splib, modfac

logging.basicConfig(level=logging.DEBUG)

# Logger
log = logging.getLogger(__name__)


# Checks whether the argument is a readable directory, used for our argparse logic
def readable_dir(dirname):
    if not os.path.isdir(dirname):
        raise argparse.ArgumentTypeError("Input path {0} is not a directory".format(dirname))
    if not os.access(dirname, os.R_OK):
        raise argparse.ArgumentTypeError("Input path {0} is not readable".format(dirname))
    return dirname

# Splits the input into latitude/longitude pairs
def parse_lat_lons(coordinate_list):
    n = len(coordinate_list)
    if n % 2:
        log.info("Odd number of point components encountered... omitting the last latitude")
        coordinate_list = coordinate_list[0:n - 1]
    return [(float(coordinate_list[2 * i + 1]) % 360, float(coordinate_list[2 * i])) for i in range(n / 2)]
    #      [(lon1, lat1), (lon2, lat2), ...]
    # % 360 on longitudes, to map negative longitudes to 180...360. Still doesn't handle polygons
    # over the 0-meridian


# read geoJSON, convert to a shapely polygon
# TODO: Currently returns the first polygon, if many are defined
def read_poly_file(polyfile):
    try:
        with open(polyfile) as f:
            print ('Reading polygon from file', polyfile)
            js = json.load(f)
            for feature in js['features']:
                polygon = shapely.geometry.shape(feature['geometry'])
                print('  Found polygon', polygon)
                return polygon            
    except Exception as e:
        print('Failed to read or parse the polygon file:', polyfile, e)
        sys.exit(1)
    

# Main function
def main():
    les_types = [modfac.dales_type, modfac.dummy_type, modfac.ncbased_type]
    gcm_types = [modfac.oifs_type, modfac.dummy_type, modfac.ncbased_type]

    parser = argparse.ArgumentParser(description="GCM-LES superparametrization run script", fromfile_prefix_chars="@")

    parser.add_argument("--steps", dest="gcm_steps",
                        metavar="N",
                        type=int,
                        default=splib.gcm_num_steps,
                        help="Nr. of (GCM) time steps")

    parser.add_argument("--conf", dest="conf",
                        metavar="FILE.json",
                        type=str,
                        default=None,
                        help="Configuration file")

    parser.add_argument("--lesdir", dest="les_input_dir",
                        metavar="DIR",
                        type=readable_dir,
                        default=splib.les_input_dir,
                        help="LES input directory")

    parser.add_argument("--lestype", dest="les_type",
                        metavar="TYPE",
                        choices=les_types,
                        type=str,
                        default=splib.les_type,
                        help="LES model type")

    parser.add_argument("--lesprocs", dest="les_num_procs",
                        metavar="N",
                        type=int,
                        default=splib.les_num_procs,
                        help="Nr. of MPI tasks per LES")

    parser.add_argument("--les_dt", dest="les_dt",
                        metavar="dt",
                        type=float,
                        default=60,
                        help="Time step (s) between saving LES statistics e.g. LWP fields")

    parser.add_argument("--spinup", dest="les_spinup",
                        metavar="T",
                        type=int,
                        default=0,
                        help="Time for initial spinup of LES models to initial profile")

    parser.add_argument("--spinup_steps", dest="les_spinup_steps",
                        metavar="N",
                        type=int,
                        default=1,
                        help="Number of iterations for spinup nudging")

    parser.add_argument("--spinup_forcing", dest="les_spinup_forcing_factor",
                        metavar="f",
                        type=float,
                        default=1.,
                        help="Forcing strength during les spinup")

    parser.add_argument("--gcmdir", dest="gcm_input_dir",
                        metavar="DIR",
                        type=readable_dir,
                        default=splib.gcm_input_dir,
                        help="GCM input directory")

    parser.add_argument("--gcmtype", dest="gcm_type",
                        metavar="TYPE",
                        choices=gcm_types,
                        type=str,
                        default=splib.gcm_type,
                        help="GCM model type")

    parser.add_argument("--gcmprocs", dest="gcm_num_procs",
                        metavar="N",
                        type=int,
                        default=splib.gcm_num_procs,
                        help="Nr. of MPI tasks for gcm")

    parser.add_argument("--gcmexp", dest="gcm_exp_name",
                        metavar="NAME",
                        type=str,
                        default=splib.gcm_exp_name,
                        help="GCM experiment name")

    parser.add_argument("--odir", dest="output_dir",
                        metavar="DIR",
                        type=str,
                        default=splib.output_dir,
                        help="Output directory")

    parser.add_argument("--dryrun", action="store_true",
                        default=False,
                        help="Only initialize the GCM, save the grid point coordinates to a file.")

    parser.add_argument("--points", metavar="lat1 lon1 ... latn lonn",
                        nargs="+",
                        default="",
                        help="Space-separated list of lat/lon pairs. To fill all boxes, use \"all\"")

    parser.add_argument("--poly", metavar="lat1 lon1 ... latn lonn",
                        nargs="+",
                        default="",
                        help="Space-separated list of polygon corner lat/lon pairs.")

    parser.add_argument("--polyfile", metavar="filename",
                        default=None,
                        help="geoJSON file containing a polygon for superparameterization.")

    parser.add_argument("--output_poly", metavar="lat1 lon1 ... latn lonn",
                        nargs="+",
                        default="",
                        help="Designate non-superparametrized columns for inclusion in netCDF output. Space-separated "
                             "list of polygon corner lat/lon pairs.")

    parser.add_argument("--output_polyfile", metavar="filename",
                        default=None,
                        help="geoJSON file containing a polygon for statistics output")

    parser.add_argument("-a", "--all", action="store_true",
                        default=False,
                        help="Superparametrize all IFS grid columns")

    parser.add_argument("--numles", dest="max_num_les",
                        metavar="N",
                        type=int,
                        default=-1,
                        help="Nr. of LES instances to run. If a single point is selected,"
                             "use this option to select a given number of closest gridpoints")

    parser.add_argument("--queue", dest="les_queue_threads",
                        metavar="N",
                        type=int,
                        default=splib.les_queue_threads,
                        help="Nr. of LES models to run concurrently. 1 denotes serial execution. Default: fully "
                             "parallel")

    parser.add_argument("--channel", dest="channel_type",
                        metavar="TYPE",
                        choices=["mpi", "sockets", "nospawn"],
                        type=str,
                        default=splib.channel_type,
                        help="Amuse communication type")

    parser.add_argument("--restart", action="store_true",
                        default=False,
                        help="Restart an old run")

    parser.add_argument("--cplsurf", action="store_true",
                        default=False,
                        help="Couple surface fluxes and roughness lengths")

    parser.add_argument("--qt_forcing", dest="qt_forcing",
                        metavar="TYPE",
                        choices=["sp", "variance", "local"],
                        type=str,
                        default="sp",
                        help="qt forcing type on LES (stimulate ql alignment with variance or local forcing)")


    args = parser.parse_args()

    geometries = []
    for p in parse_lat_lons(args.points):
        geometries.append(shapely.geometry.Point(p))
    if any(parse_lat_lons(args.poly)):
        geometries.append(shapely.geometry.Polygon(parse_lat_lons(args.poly)))
    if args.all:
        geometries = [shapely.geometry.box(-float("inf"), -float("inf"), float("inf"), float("inf"))]
    # Read geoJSON from file, convert to a shapely object
    if args.polyfile:
        p = read_poly_file(args.polyfile)        
        geometries.append(p)
        
    output_geometries = []
    if any(parse_lat_lons(args.output_poly)):
        output_geometries.append(shapely.geometry.Polygon(parse_lat_lons(args.output_poly)))
    if args.output_polyfile:
        p = read_poly_file(args.output_polyfile)        
        output_geometries.append(p)

        
    splib.read_config(args.conf)

    splib.initialize(args.__dict__, geometries, output_geometries)

    splib.run(args.gcm_steps)

    splib.finalize()


if __name__ == "__main__":
    print("-- Spmaster starting --")
    main()
