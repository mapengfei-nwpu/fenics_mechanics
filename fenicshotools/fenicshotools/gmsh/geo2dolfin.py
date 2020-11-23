"FIXME"

from mpi4py import MPI as mpi
from fenics import *
from ..io import save_geometry, load_geometry
from .gmsh2dolfin import gmsh2dolfin
from .gmshfile import GmshFile
from .inline_backend import gmsh_cpp_geo2msh

import os
import tempfile
import shutil

__all__ = ["geo2dolfin"]

def geo2dolfin(code,topological_dimension=3,
               geometric_dimension=None,
               comm=None,marker_ids=None):
    """
    FIXME
    """

    # We generate the msh file from a geo file, and then feed into gmsh2dolfin
    comm = comm if comm is not None else mpi_comm_world()

    if not MPI.is_receiver(comm):
        # create a temporary directory
        tmpdir = tempfile.mkdtemp()

        # generates a new geo file
        geoname = os.path.join(tmpdir,'mesh.geo')
        with open(geoname, 'w') as f:
            f.write(code)

        # generates the mesh
        info("--- Generating .msh file from .geo (may take a while)")
        mshname = os.path.join(tmpdir,'mesh.msh')
        logname = os.path.join(tmpdir,'mesh.log')
        curdir = os.getcwd()
        os.chdir(tmpdir)
        gmsh_cpp_geo2msh(geoname,topological_dimension,mshname,logname)
        os.chdir(curdir)

        # communicate the filename
        comm.tompi4py().bcast(mshname,root=0)
    else:
        # receive the filename
        mshname = comm.tompi4py().bcast(None,root=0)

    # import the mesh
    tmpdir = os.path.dirname(mshname)
    info("--- Importing from .msh")
    mesh,phi,markers = gmsh2dolfin(mshname,comm=comm,tmpdir=tmpdir,
        geometric_dimension=geometric_dimension,marker_ids=marker_ids)

    # Clean-up the temporary directory
    if not MPI.is_receiver(comm):
        shutil.rmtree(tmpdir)

    # In case of parametrized meshes, we return the backbone mesh and the
    # displacement, rather than the parametrized mesh obtained with
    # create_mesh utility.  Indeed, I don't know how to retrieve the backbone
    # mesh from the parametrized one in a trivial way.
    return mesh,phi,markers

