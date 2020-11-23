"""
Module for saving and loading a geometry to and from a dolfin.HDF5File.
"""
try:
    import h5py
    has_h5py = True
except:
    has_h5py = False

from hashlib import sha1
import os
from fenics import *
import fenics
import re

__all__ = ["save_geometry", "load_geometry", "save_function", "load_function"]

# If dolfin compiled with petsc4py support we expect different types for the mpi_comm.
mpi_comm_type = type(fenics.mpi_comm_world())
str_repr_mpi_comm_type = mpi_comm_type.__name__
rst_repr_mpi_comm_type = "petsc4py.PETSc.Comm" if fenics.has_petsc4py() \
                         else "dolfin.MPI_Comm"

def save_geometry(comm,mesh,phi,h5name,h5group='',markers=None,\
                  overwrite_file=False,overwrite_group=True):

    if not isinstance(comm, mpi_comm_type):
        raise TypeError("expected  a '{}' for the first argument".format(\
            str_repr_mpi_comm_type))

    if not isinstance(mesh, Mesh):
        raise TypeError("expected  a 'ufl.Domain' for the second argument")

    if not isinstance(h5name, str):
        raise TypeError("expected a 'str' as the third argument")

    if not(len(h5name)>3 and h5name[-3:]==".h5"):
        raise ValueError("expected 'h5name' argument to end with '.h5'")

    if not isinstance(h5group, str):
        raise TypeError("expected the fourth argument to be a 'str'")

    markers = markers or {}

    if not (isinstance(markers, dict), ):
        raise TypeError("expected the 'markers' argument to be a 'dict' "\
                        "with positive int as keys and a size 2 tuple with str "\
                        "and a positive int as value.")

    file_mode = "a" if os.path.isfile(h5name) and not overwrite_file else "w"

    # IF we should append the file but overwrite the group we need to
    # check that the group does not exist. If so we need to open it in
    # h5py and delete it.
    if file_mode == "a" and overwrite_group and h5group!="":
        if has_h5py:
            with h5py.File(h5name) as h5file:
                if h5group in h5file:
                    debug("Deleting existing group: '{}'".format(h5group))
                    del h5file[h5group]
        else:
            raise RuntimeError("Cannot overwrite group '{}'. "\
                               "Need h5py for that.".format(h5group))
    # Open the file for writing
    with HDF5File(comm, h5name, file_mode) as h5file:
        ggroup = '{}/geometry'.format(h5group)

        # Save the mesh
        mgroup = '{}/mesh'.format(ggroup)
        h5file.write(mesh, mgroup)
        h5file.attributes(ggroup)['label'] = mesh.label()
        h5file.attributes(mgroup)['label'] = mesh.label()
        h5file.attributes(mgroup)['name']  = mesh.name()

        # Save high-order information if available
        if phi is not None:
            order  = phi.ufl_element().degree()
            family = phi.ufl_element().family()
            fspace = '{}_{}'.format(family, order)
            pgroup = '{}/coordinates'.format(ggroup)
            h5file.write(phi, pgroup)
            h5file.attributes(pgroup)['function_space'] = fspace
            h5file.attributes(pgroup)['label'] = phi.label()
            h5file.attributes(pgroup)['name'] = phi.name()

        # Save the boundary markers as MeshValueCollection, that is suitable
        # for reading in parallel. MeshDomains is not supported by HDF5File
        # since DOLFIN 2016.1.0.
        for dim,(mvc,mm) in markers.items():
            # Skip if entities are not marked
            if mvc is None: continue
            # Save it similarly to MeshDomain
            dgroup = '{}/mesh/domain_{}'.format(ggroup,dim)
            h5file.write(mvc,dgroup)
            # We save the labels of the markers as HDF5 attributes
            if mm is None: continue
            for marker,lbl in mm.items():
                aname = 'marker_name_{}'.format(marker)
                h5file.attributes(dgroup)[aname] = str(lbl)


save_geometry.func_doc = """
    Save a geometry to a h5 file

    *Arguments*
      comm (:py:class:`{}`)
        An MPI communicator used to partioned the domain.
      domain (:py:class:`ufl.Domain`)
        The domain which should be saved.
      h5name (str)
        The name of the h5 file that the domain should be saved to.
      h5group (str)
        A h5 group name the domain optionally should be stored to
      markers (dict)
        A dict with marker information. Expects positiv int as keys and a size
        2 tuple with str and a positive int as value.
        
    """.format(rst_repr_mpi_comm_type)


def load_geometry(comm,h5name,h5group=''):

    if not isinstance(comm, mpi_comm_type):
        raise TypeError("expected  a '{}' for the first argument".format(\
            str_repr_mpi_comm_type))

    if not isinstance(h5name, str):
        raise TypeError("expected a 'str' as the second argument")

    if not os.path.isfile(h5name):
        raise IOError("expected 'h5name' to be a file.")

    if not isinstance(h5group, str):
        raise TypeError("expected the third argument to be a 'str'")

    # Open the file
    with HDF5File(comm, h5name, 'r') as h5file:
        if h5group and not h5file.has_dataset(h5group):
            raise ValueError("The '{}' file does not have the '{}' dataset".\
                             format(h5name, h5group))

        ggroup = '{}/geometry'.format(h5group)
        if not h5file.has_dataset(ggroup):
            raise ValueError("Expected the '{}' file to have the '{}' dataset".\
                             format(h5name, ggroup))

        mgroup = '{}/mesh'.format(ggroup)
        if not h5file.has_dataset(mgroup):
            raise ValueError("Expected the '{}' file to have the '{}' dataset".\
                             format(h5name, mgroup))

        # Read attributes
        ggroup_attrs = h5file.attributes(ggroup)
        mgroup_attrs = h5file.attributes(mgroup)
        if "label" not in ggroup_attrs:
            raise ValueError("Expected 'label' to be an attribute in '{}'".format(ggroup))
        if "label" not in mgroup_attrs:
            raise ValueError("Expected 'label' to be an attribute in '{}'".format(mgroup))
        if "name" not in mgroup_attrs:
            raise ValueError("Expected 'name' to be an attribute in '{}'".format(mgroup))

        dlabel = ggroup_attrs['label']
        mlabel = mgroup_attrs['label']
        mname  = mgroup_attrs['name']

        # Read the mesh
        mesh = Mesh(comm)
        h5file.read(mesh, mgroup, False)
        mesh.rename(mname, mlabel)

        # Load possible high-order information
        pgroup = '{}/geometry/coordinates'.format(h5group)
        if h5file.has_dataset(pgroup):

            # Read string labels
            pgroup_attrs = h5file.attributes(pgroup)
            for attr in ["function_space", "name", "label"]:
                if attr not in pgroup_attrs:
                    raise ValueError("Expected '{}' to be an attribute in '{}'".\
                                     format(attr, pgroup))

            fspace = pgroup_attrs['function_space']
            pname  = pgroup_attrs['name']
            plabel = pgroup_attrs['label']

            # Recontruct function space
            family, order = fspace.split('_')
            gdim = mesh.geometry().dim()
            V = VectorFunctionSpace(mesh, family, int(order), gdim)

            # Read the coordinates
            phi = Function(V)
            h5file.read(phi, pgroup)
            phi.rename(pname, plabel)
        else:
            phi = None

        # Load the boundary markers
        markers = {}
        for dim in range(mesh.ufl_domain().topological_dimension()+1):
            dgroup = '{}/mesh/domain_{}'.format(ggroup, dim)

            # If dataset is not present
            if not h5file.has_dataset(dgroup):
                markers[dim] = (None,None)
                continue

            # Read MeshValueCollection from the file
            mvc = MeshValueCollection('size_t')
            mvc.init(mesh,dim)
            h5file.read(mvc,"{}/mesh/domain_{}".format(ggroup,dim))

            # Read the attributes
            labels = {}
            for aname in h5file.attributes(dgroup).str().strip().split(' '):
                m = re.match('marker_name_(\d+)', aname)
                if not m: continue
                marker = m.group(1)
                name = h5file.attributes(dgroup)[m.group(0)]
                if name == 'None': name = None
                labels[int(marker)] = name
            if len(labels)==0:
                labels = None

            markers[dim] = (mvc,labels)

    return mesh,phi,markers

load_geometry.func_doc = """
    Load a dolfin mesh stored in a h5 file

    *Arguments*
      comm (:py:class:`{}`)
        An MPI communicator used to partioned the domain.
      h5name (str)
        The name of the h5 file that the domain should be load from.
      h5group (str)
        A h5 group name the domain optionally should be load from.
    """.format(rst_repr_mpi_comm_type)

def save_function(comm, func, h5name, h5group='', index=0,\
                  overwrite_file=False, overwrite_group=True):

    if not isinstance(comm, mpi_comm_type):
        raise TypeError("expected  a '{}' for the first argument".format(\
            str_repr_mpi_comm_type))

    if not isinstance(func, Function):
        raise TypeError("expected  a 'dolfin.Function' for the second argument")

    if not isinstance(h5name, str):
        raise TypeError("expected a 'str' as the third argument")

    if not(len(h5name)>3 and h5name[-3:]==".h5"):
        raise ValueError("expected 'h5name' argument to end with '.h5'")

    if not isinstance(h5group, str):
        raise TypeError("expected the fourth argument to be a 'str'")

    if not isinstance(index, int) or index < 0:
        raise TypeError("expected a positive 'int' for the fifth argument")

    # Open a file
    file_mode = "a" if os.path.isfile(h5name) and not overwrite_file else "w"

    data_name = "{}_{0:03d}".format(h5group, index) \
                if h5group else "data_{0:03d}".format(index)

    # If no h5py we just save the function directly
    if not has_h5py:

        # Write file without link
        with HDF5File(comm, h5name, file_mode) as h5file:
            debug("Write function to '{}' dataset.".format(data_name))
            h5file.write(func, data_name)

    # If h5py we create a link to function data that does not change
    else:

        # IF we should append the file but overwrite the group we need to
        # check that the group does not exist. If so we need to open it in
        # h5py and delete it.
        if file_mode == "a" and overwrite_group and data_name!="":
            with h5py.File(h5name) as h5file:
                if data_name in h5file:
                    del h5file[data_name]

        with HDF5File(comm, h5name, file_mode) as h5file:

            # Create "good enough" hash. This is done to avoid data
            # corruption when restarted from different number of
            # processes, different distribution or different function
            # space
            local_hash = sha1()
            local_hash.update(str(func.function_space().mesh().num_cells()))
            local_hash.update(str(func.function_space().ufl_element()))
            local_hash.update(str(func.function_space().dim()))
            local_hash.update(str(MPI.size(comm)))

            # Global hash (same on all processes), 10 digits long
            global_hash = MPI.sum(comm, int(local_hash.hexdigest(), 16))
            global_hash = "hash_" + str(int(global_hash%1e10)).zfill(10)

            # Create hashed named for writing the function the first time
            if not h5file.has_dataset(global_hash):
                debug("Write function to hashed dataset.")
                h5file.write(func, global_hash)

            else:

                # If hashed data set already excist just write the vector
                if h5file.has_dataset(data_name):
                    raise IOError("Cannot write dataset '{}'. It "\
                                  "already excist.".format(data_name))
                debug("Write vector to '{}' dataset.".format(data_name))
                h5file.write(func.vector(), "{}/vector_0".format(data_name))

        # Reopen file with h5py to create link, only on rank 0 process
        if MPI.rank(comm) == 0:
            with h5py.File(h5name, "a") as h5file:
                assert global_hash in h5file

                # Get dataset to link to
                link_to = h5file[global_hash]
                
                # If data set created for the first time
                if data_name not in h5file:
                    link_from = h5file.create_group(data_name)

                    # Create a hard link to the vector created in the
                    # global hash
                    link_from["vector_0"] = link_to["vector_0"]
                else:
                    link_from = h5file[data_name]

                debug("Link '{}' datasets in h5 file.".format(data_name))
                for dname in ['x_cell_dofs', 'cell_dofs', 'cells']:

                    # Create a hard link to the Function data
                    link_from[dname] = link_to[dname]
                
        MPI.barrier(comm)

save_function.func_doc = """
    Save a dolfin Function to an h5 file

    *Arguments*
      comm (:py:class:`{}`)
        An MPI communicator used to partioned the domain.
      func (:py:class:`dofin.Function`)
        The function that will be saved
      h5name (str)
        The name of the h5 file that the function should be saved to.
      h5group (str)
        A h5 group name the function optionally should be saved from.
      index (int)
        A saving index
    """.format(rst_repr_mpi_comm_type)

def load_function(comm, func, h5name, h5group='', index=0):
    if not isinstance(comm, mpi_comm_type):
        raise TypeError("expected  a '{}' for the first argument".format(\
            str_repr_mpi_comm_type))
    
    if not isinstance(func, Function):
        raise TypeError("expected  a 'dolfin.Function' for the second argument")

    if not isinstance(h5name, str):
        raise TypeError("expected a 'str' as the third argument")

    if not(len(h5name)>3 and h5name[-3:]==".h5"):
        raise ValueError("expected 'h5name' argument to end with '.h5'")

    if not isinstance(h5group, str):
        raise TypeError("expected the fourth argument to be a 'str'")

    if not isinstance(index, int) or index < 0:
        raise TypeError("expected a positive 'int' for the fifth argument")
    
    data_name = "{}_{0:03d}".format(h5group, index) \
                if h5group else "data_{0:03d}".format(index)

    
    # Write file without link
    with HDF5File(comm, h5name, "r") as h5file:
        debug("Load function from '{}' dataset.".format(data_name))
        h5file.read(func, data_name)
    
load_function.func_doc = """
    Load a dolfin Function to a h5 file

    *Arguments*
      comm (:py:class:`{}`)
        An MPI communicator used to partioned the domain.
      func (:py:class:`dofin.Function`)
        The function that will be loaded
      h5name (str)
        The name of the h5 file that the function should be loaded from.
      h5group (str)
        A h5 group name the function optionally should be loaded from.
      index (int)
        A saving index
    """.format(rst_repr_mpi_comm_type)
