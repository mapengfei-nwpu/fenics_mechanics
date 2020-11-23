from fenics import mpi_comm_world
from fenicshotools.io import save_geometry
from fenicshotools.gmsh import geo2dolfin

import argparse
import os.path

def get_args() :
    """ Parse command line arguments"""
    descr = 'Convert a .geo GMSH geometry file to a FEniCS .h5 file.'
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('ifile', metavar='myfile.msh',
            type=str, 
            help='GMSH .geo filename.')
    parser.add_argument('-o', metavar='myfile.h5',
            type=str,
            help='FEniCS .h5 filename.')
    # topological dimension must be provided since in general
    # GMSH is not able to deduce it from the geometry file.
    parser.add_argument('-t', '--topological_dimension',
            type=int, default=3,
            help='Topological dimension of the mesh.')
    parser.add_argument('-g', '--geometric_dimension',
            type=int, default=None,
            help='Geometric embedding dimension of the mesh.')

    return parser

def main_func():

    comm = mpi_comm_world()
    args = get_args().parse_args()
    ifile = args.ifile
    ofile = args.o or '{}.h5'.format(os.path.splitext(args.ifile)[0])

    # 1. parsing
    with open(ifile, 'r') as f:
        code = f.read()

    # 2. convert
    mesh,phi,markers = geo2dolfin(code,args.topological_dimension,
                                  args.geometric_dimension,comm)

    # 3. writing
    save_geometry(comm,mesh,phi,ofile,'',markers=markers,overwrite_file=True)

