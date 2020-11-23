from fenicshotools.gmsh import GmshFile
from fenicshotools.gmsh import gmsh2dolfin
import argparse
from dolfin import *

def gmshfiletype(filename):
    try:
        gmsh = GmshFile(filename)
    except Exception as e:
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

    parameters["allow_extrapolation"] = True

    args = get_args()
    gmsh = args.ifile

    gmsh.print_info()

    mesh,phi,markers = gmsh2dolfin(gmsh,use_coords=True)

    if args.plot :
        if phi is not None :
            meshr = refine(refine(mesh))
            Vr = VectorFunctionSpace(meshr,"CG",1)
            phir = interpolate(phi,Vr)
            ALE.move(meshr,phir)
        else :
            meshr = mesh
        plot(meshr,interactive=True)


