from fenicshotools.io import save_geometry
from fenicshotools.gmsh import gmsh2dolfin
from fenics import mpi_comm_world

import argparse
import os.path

def get_args():
    descr = 'Convert a .msh GMSH mesh file to a FEniCS .h5 file.'
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('ifile', metavar='myfile.msh',
            type=str, 
            help='GMSH .msh filename.')
    parser.add_argument('-o', metavar='myfile.h5',
            type=str,
            help='FEniCS .h5 filename.')
    parser.add_argument('-g', '--geometric_dimension',
            type=int, default=None,
            help='Geometric embedding dimension of the mesh.')

    return parser


def main_func():

    # 1. parsing
    args = get_args().parse_args()
    ifile = args.ifile
    ofile = args.o or '{}.h5'.format(os.path.splitext(args.ifile)[0])

    # 2. conversion
    comm = mpi_comm_world()
    gdim = args.geometric_dimension
    mesh,phi,markers = gmsh2dolfin(ifile,comm=comm,geometric_dimension=gdim)

    # 3. writing
    save_geometry(comm,mesh,phi,ofile,'',markers=markers,overwrite_file=True)

