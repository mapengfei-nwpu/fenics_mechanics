from fenicshotools.gmsh import GmshFile
from fenicshotools.gmsh import gmsh2dolfin
from fenicshotools.io import XdmfHoFile
import argparse
from dolfin import *

def gmshfiletype(filename) :
    try :
        gmsh = GmshFile(filename)
    except Exception as e :
        raise argparse.ArgumentTypeError(e)
    return gmsh

def get_args() :
    parser = argparse.ArgumentParser(description='Parse gmsh file.')
    parser.add_argument('ifile', metavar='file.msh',
            type=gmshfiletype, 
            help='Gmsh input filename.')
    parser.add_argument('--plot', action='store_true', help='Plot the mesh.')
    return parser.parse_args()

if __name__ == "__main__" :

    args = get_args()
    gmsh = args.ifile
    #gmsh.print_info()

    mesh,phi,markers = gmsh2dolfin(gmsh, use_coords=True)

    if phi is not None:
        order = phi.function_space().ufl_element().degree()
        #domain = create_mesh(phi)
        #V = FunctionSpace(domain,"CG",order)
        #u = interpolate(Expression("x[0]*x[0]-x[1]*x[1]"),V)

        # export
        with XdmfHoFile("test.xdmf") as ofile:
            ofile.add("mesh",mesh)
            #ofile.add("field", u)
            #ofile.add("field", u)
