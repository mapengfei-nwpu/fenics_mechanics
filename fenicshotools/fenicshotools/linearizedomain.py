"""
FIXME
"""
from fenics import *
from .gmsh import geo2dolfin
import numpy as np

__all__ = [ 'linearize_domain_and_fields' ]

def _get_refined_connectivity(domain):
    """
    Overcomplicate procedure to obtain a linear mesh of a reference
    finite element using all the dof (of a Lagrangian space).
    """
    V = domain.coordinates().function_space()
    gdim = domain.geometric_dimension()
    tdim = domain.topological_dimension()
    if tdim == 1:
        meshref = UnitIntervalMesh(1)
    elif tdim == 2:
        meshref = UnitTriangleMesh()
    elif tdim == 3:
        meshref = UnitTetrahedronMesh()
    else:
        raise NotImplementedError

    # reference element with coordinates in the
    # right order
    degree = V.ufl_element().degree()
    Vref = FunctionSpace(meshref, 'P', degree)
    coords = Vref.dofmap().tabulate_all_coordinates(meshref).reshape(-1,gdim)
    cell_dofs = Vref.dofmap().cell_dofs(0)
    idx = range(len(cell_dofs))
    idx.sort(key=lambda i: cell_dofs[i])
    coords = coords[idx,:]

    if tdim == 1:
        coords_i = range(degree+1)
        conn_i = [ (i,i+1) for i in range(degree) ]

        perm = { coords_i.index(int(degree*v)): i
                 for i, v in enumerate(coords) }
        conn = [ map(coords_i.index, c) for c in conn_i ]
        conn = np.vectorize(perm.get, otypes=[np.uintp])(conn)

    elif tdim == 2:
        coords_ij = [ (i,j) for i in range(degree+1)
                            for j in range(i+1) ]

        conn_ij  = [ [(i,j), (i+1,j), (i+1,j+1)]
                     for i in range(degree)
                     for j in range(i+1) ]
        conn_ij += [ [(i,j), (i,j+1), (i+1,j+1)]
                     for i in range(1,degree)
                     for j in range(i) ]

        A = np.array([[1, 1], [0, 1]])
        perm = { coords_ij.index(tuple(np.rint(degree*A.dot(v)).astype(int))): i
                 for i, v in enumerate(coords) }
        conn = [ map(coords_ij.index, c) for c in conn_ij ]
        conn = np.vectorize(perm.get, otypes=[np.uintp])(conn)

    elif tdim == 3:
        # coordinate (i,j,k) indices
        coords_ijk = [ (i,j,k) for i in range(degree+1)
                               for j in range(i+1)
                               for k in range(j+1) ]
        # connectivity
        conn_ijk  = [ [(i,j,k), (i+1,j,k), (i+1,j+1,k), (i+1,j+1,k+1)]
                      for i in range(degree)
                      for j in range(i+1) for k in range(j+1) ]

        conn_ijk += [ [(i,j,k), (i,j+1,k+1), (i+1,j+1,k+1), (i,j+1,k)]
                      for i in range(1,degree)
                      for j in range(i) for k in range(j+1) ]
        conn_ijk += [ [(i,j,k), (i,j+1,k), (i+1,j+1,k+1), (i+1,j+1,k)]
                      for i in range(1,degree)
                      for j in range(i) for k in range(j+1) ]
        conn_ijk += [ [(i,j+1,k), (i,j+1,k+1), (i+1,j+2,k+1), (i+1,j+1,k+1)]
                      for i in range(1,degree)
                      for j in range(i) for k in range(j+1) ]
        conn_ijk += [ [(i+1,j+1,k), (i+1,j+1,k+1), (i+1,j+2,k+1), (i,j+1,k)]
                      for i in range(1,degree)
                      for j in range(i) for k in range(j+1) ]

        conn_ijk += [ [(i,j,k), (i,j,k+1), (i,j+1,k+1), (i+1,j+1,k+1)]
                      for i in range(2,degree)
                      for j in range(i) for k in range(j) ]

        A = np.array([[1, 1, 1], [0, 1, 1], [0, 0, 1]])
        perm = { coords_ijk.index(tuple(np.rint(degree*A.dot(v)).astype(int))): i
                 for i, v in enumerate(coords) }
        conn = [ map(coords_ijk.index, c) for c in conn_ijk ]
        conn = np.vectorize(perm.get, otypes=[np.uintp])(conn)

    else:
        raise NotImplementedError

    return conn

def _mesh_from_conn(coords, conn):
    """
    Generate a DOLFIN mesh from an array
    of connectivity and coordinate of the
    vertices.
    """
    mesh = Mesh()
    editor = MeshEditor()
    editor.open(mesh, coords.shape[1], conn.shape[1]-1)
    editor.init_vertices_global(coords.shape[0], coords.shape[0])
    for i, v in enumerate(coords):
        editor.add_vertex(i, v)
    editor.init_cells_global(conn.shape[0], conn.shape[0])
    for i, c in enumerate(conn):
        editor.add_cell(i, c)
    editor.close()
    return mesh

def linearize_domain_and_fields(domain, *args, **kwargs):
    """
    Linearize a domain

    Given a (generally non-linear) domain, and an arbitrary number of
    functions, this method returns a linearized domain obtained by
    first meshing the dofs and then further refining it by splitting
    the elements (ndiv times).  The functions are then interpolated as
    linear function on the new mesh.

    *Arguments*
      domain (FIXME)
        FIXME
    """
    # 0. parameters
    ndiv = kwargs['ndiv'] if 'ndiv' in kwargs else 0

    # 1. build a mesh with dofs as vertices
    gdim = domain.geometric_dimension()
    tdim = domain.topological_dimension()
    mesh = domain.data()
    if domain.coordinates():
        V = domain.coordinates().function_space()
        # coordinates of the dofs
        # we do not use actual coordinates provided
        # by the coordinate function, since we wish
        # the refine the domain even further and then
        # interpolate everything.
        coords = V.dofmap().tabulate_all_coordinates(mesh).reshape(-1,gdim)
        idx = np.hstack([ V.sub(i).dofmap().dofs()[:,np.newaxis]
                      for i in range(V.num_sub_spaces()) ])
        coords = coords[idx[:,0],:]
        # connectivity
        cell_dofs = np.vstack([ V.sub(0).dofmap().cell_dofs(i)
                                for i in range(mesh.num_cells()) ])
        perm = { v: k for k, v in enumerate(idx[:, 0]) }
        cell_dofs = np.vectorize(perm.get, otypes=[np.uintp])(cell_dofs)
        # apply the refinement rule
        refrule = _get_refined_connectivity(domain)
        conn = []
        for k, cell in enumerate(cell_dofs):
            conn += [ cell[refrule] ]
        conn = np.vstack(conn)
        # construct the mesh from connectivity
        refmesh = _mesh_from_conn(coords, conn)
    else:
        # no refinement is required
        refmesh = mesh

    # 2. additional refinement (if requested)
    refmesh = reduce(lambda m, _: refine(m), range(ndiv), refmesh)

    # 3. interpolation of all the fields
    fields  = list(args)
    fields += [ domain.coordinates() ] if domain.coordinates() else []
    I = LagrangeInterpolator()
    ifields = []
    for f0 in fields:
        # construction of linear space
        elm = f0.function_space().ufl_element()
        shape = elm.value_shape()
        if len(shape) == 0:
            Vref = FunctionSpace(refmesh, 'P', 1)
        elif len(shape) == 1:
            Vref = VectorFunctionSpace(refmesh, 'P', 1, shape[0])
        else:
            # FIXME symmetry?
            Vref = TensorFunctionSpace(refmesh, 'P', 1, shape=shape)
        # interpolation
        fun = Function(Vref)
        I.interpolate(fun, f0)
        ifields.append(fun)

    # 4. remap the coordinates if necessary
    if domain.coordinates():
        icoords = ifields.pop()
        vvals = icoords.compute_vertex_values(refmesh).reshape((gdim,-1))
        refmesh.coordinates()[:] = vvals.transpose()

    # 5. we return the mesh and the fields
    if len(fields) == 0:
        return refmesh
    else:
        return tuple([refmesh] + ifields)

