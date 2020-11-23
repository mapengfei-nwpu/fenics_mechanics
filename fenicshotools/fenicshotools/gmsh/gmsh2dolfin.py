"FIXME"
from dolfin import *
import numpy as np
from six.moves import zip
from .generate_points import *
from .make_affine_mapping import *
from .reference_topology import *
from .gmshfile import GmshFile
from ..io import save_geometry,load_geometry

from scipy.spatial import cKDTree

import tempfile,os,shutil

__all__ = ["gmsh2dolfin"]

def _tabulate_all_gmsh_dofs(gmsh):

    # Gmsh provides the position of the dofs in the deformed state
    # and not in the reference.
    # Here we construct the list of all reference dofs in the same
    # order of the nodes.

    # the type of the cell
    cell_shape = gmsh.cell_shape()
    vert_per_cell = gmsh_shape[cell_shape].num_vertices()

    # cell dofs in the reference element
    cell_ref_dofs = np.array(generate_points(cell_shape, gmsh.mesh_order))
    cell_ref_vert = cell_ref_dofs[:vert_per_cell]

    # cell with dofs
    cells, cells_dofs = gmsh.cells(), gmsh.cells(True)
    nodes = gmsh.dof_coords()
    global_dofs = np.zeros(nodes.shape)

    for cell_vert,cell_dofs in zip(nodes[cells],cells_dofs):
        # create the affine mapping
        A, b = make_affine_mapping(cell_ref_vert, cell_vert)
        # maps the dofs
        dofs = cell_ref_dofs.dot(A.T) + np.array([b,]*cell_ref_dofs.shape[0])
        global_dofs[cell_dofs,:] = dofs

    return global_dofs


def gmsh2dolfin(ifile,use_coords=True,comm=None,
                geometric_dimension=None,tmpdir=None,
                marker_ids=None):

    # Use provided communicator, or world otherwise
    comm = comm if comm is not None else mpi_comm_world()

    # NOTE: markers are computed for the global mesh, with global indices.
    # In parallel, we could use build_distributed_value_collection from
    # MeshPartitioning.h, but this is not exposed to Python.
    # MeshDomains is not supported anymore for distributed meshes.
    # A dirty workaround is to save the MeshValueCollection to hdf5 file, then
    # read it again.
    # We do the same for the coordinate function and the mesh, without using
    # the build_distributed_mesh function.

    # NOTE: with coordinate map, we need to initialize a Function with a
    # MPI_COMM_SELF communicator. This will fail if PETSc has not been
    # initialized externally. Indeed, by default PETSc used MPI_COMM_WORLD,
    # but only MPI_COMM_SELF process is doing MPI calls, so we get an error
    # somewhere. To avoid this, we force PETSc initialization here.
    SubSystemsManager.init_petsc()

    # This is equivalent to execute this portion of code only on the
    # master node (broadcaster), that on DOLFIN is rank=0.

    if not MPI.is_receiver(comm):
        # Read the file into a Gmsh object
        gmsh = GmshFile(ifile,geometric_dimension,False,marker_ids)
        # Generate the dolfin mesh
        mesh,phi,markers = _generate_mesh(gmsh,mpi_comm_self(),use_coords)
        # Create the temporary hdf5 file
        loctmpdir = tmpdir or tempfile.mkdtemp()
        h5name = os.path.join(loctmpdir,"mesh.h5")
        # save to hdf5
        save_geometry(mpi_comm_self(),mesh,phi,h5name,markers=markers)
    else:
        h5name = None

    # communicate the filename (we have a barrier here)
    h5name = comm.tompi4py().bcast(h5name,root=0)

    # Now we read (in parallel) the hdf5 file, to obtain the distributed mesh
    mesh,phi,markers = load_geometry(comm,h5name)

    # Clean-up the temporary directory
    if not MPI.is_receiver(comm):
        if tmpdir is None:
            shutil.rmtree(loctmpdir)

    return mesh,phi,markers


def _generate_mesh(gmsh,comm,use_coords=True):
    # ---------------------------
    # Construction of DOLFIN mesh
    # ---------------------------
    mesh = Mesh(comm)
    editor = MeshEditor()

    # The embedding space dimension in always 3
    # in gmsh, but not in DOLFIN.
    mdim  = gmsh.topo_dim
    gdim  = gmsh.geom_dim
    shape = gmsh.cell_shape()
    editor.open(mesh, shape, mdim, gdim)

    # vertices are not indexed in a specific order,
    # so a map is necessary to build the mesh
    msh2dolfin = {}
    vert = gmsh.vertices()
    dofs = gmsh.dof_coords()
    editor.init_vertices(len(vert))
    for dolfin_id, msh_id in enumerate(vert):
        editor.add_vertex(dolfin_id, dofs[msh_id,:])
        msh2dolfin[msh_id] = dolfin_id

    # cells
    # translate gmsh ids to dolfin ids for each cell
    cells = np.vectorize(msh2dolfin.get, otypes=[np.uintp])(gmsh.cells())
    # cell vertices needs to be sorted (ufc assumption)
    if shape not in [ "quadrilateral", "hexahedron" ]:
        cells.sort(axis=1)
    editor.init_cells(len(cells))
    for idx, vert in enumerate(cells):
        editor.add_cell(idx, vert)

    # Finalise the mesh
    editor.close(False)

    # -----------------
    # Physical entities
    # -----------------
    # Since DOLFIN 2016.2.0, MeshDomains are no longer saved to HDF5 files.
    # It is not possible to easily convert MeshDomains to MeshValueCollection
    # and viceversa, which is now the (hopefully) default and stable way to
    # store mesh informations to file.  Therefore we switch to this and we
    # keep this separate from the mesh.
    mvc = {i:(None,None) for i in range(gmsh.topo_dim+1)}
    markers = dict(gmsh.marker_name)

    # Cells markers can be directly stored since we have one-to-one mapping
    # of Gmsh indices and DOLFIN indices. Local entity id of the cell is 0.
    em = gmsh.entity_markers(0)
    mc = MeshValueCollection('size_t')
    mc.init(mesh,gmsh.topo_dim)
    for c in dolfin.cells(mesh):
        idx = c.index()
        mc.set_value(idx,0,int(em[idx]))
    mark = { m:l for m,(l,d) in markers.items() if d == gmsh.topo_dim }
    mvc[gmsh.topo_dim] = (mc,mark)


    # other entities
    for codim in range(1, gmsh.topo_dim + 1):
        entities = gmsh.entities(codim)
        if entities is None: continue
        # map the numbering of the vertices
        ent = []
        for e in entities:
            try:
                ent.append([ msh2dolfin[c] for c in e ])
            except KeyError:
                # we skip marked entities which are not belonging to
                # the mesh (e.g. external points)
                continue
        ent = np.array(ent, dtype = np.uintp)
        ent_dim = mdim - codim
        # to look for the entity we need the vertex to entity connectivity
        # map and the entity to cell map
        mesh.init(0, ent_dim)
        mesh.init(ent_dim, mdim)
        vert_to_enti = mesh.topology()(0, ent_dim)
        enti_to_cell = mesh.topology()(ent_dim, mdim)

        # the entity id is obtained by intersecating all the entities
        # connected to the given vertices.
        # The result should always be one entity.
        from functools import reduce
        eids = np.concatenate([ reduce(np.intersect1d, map(vert_to_enti, e)) \
                for e in ent ])

        mc = MeshValueCollection('size_t')
        mc.init(mesh,ent_dim)
        for e, marker in zip(eids, gmsh.entity_markers(codim)):
            mc.set_value(e,int(marker))
        mark = { m:l for m,(l,d) in markers.items() if d == ent_dim }
        mvc[ent_dim] = (mc,mark)

    if gmsh.mesh_order == 1:
        return mesh,None,mvc

    # -----------------------
    # Construction of the map
    # -----------------------
    # we need to compare the list of dofs generated by gmsh and
    # the list from ufc.
    fe = VectorElement("P",mesh.ufl_cell(),gmsh.mesh_order,gdim)
    V = FunctionSpace(mesh,fe)

    # first I need to reorganize the ufc dofs
    idx = np.column_stack([ V.sub(i).dofmap().dofs()
                for i in range(0, gdim) ])
    ufc_dofs = V.tabulate_dof_coordinates().reshape(-1,gdim)
    ufc_dofs = ufc_dofs[idx[:,0],:]

    # gmsh dofs are generated from topology
    msh_dofs = _tabulate_all_gmsh_dofs(gmsh)

    # now we compute the permutation
    tree = cKDTree(msh_dofs)
    msh2ufc = tree.query(ufc_dofs)[1]

    phi = Function(V)
    if use_coords:
        displ = gmsh.dof_coords()
    else:
        displ = gmsh.dof_coords() - msh_dofs

    for i in range(0,gdim):
        phi.vector()[idx[:,i].copy()] = (displ[msh2ufc,i]).astype(np.float)

    return mesh,phi,mvc

