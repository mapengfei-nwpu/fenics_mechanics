# vtk is required only for high-order plot
try :
    from vtk import *
    from tvtk.array_handler import *
    hasvtk = True
    import numpy as np
    # check for old vtk
    from distutils.version import StrictVersion
    vtkversion = vtk.vtkVersion.GetVTKVersion()
    hasvtk6 = StrictVersion(vtkversion) >= StrictVersion('6.0.0')
except ImportError :
    hasvtk = False

__all__ = [ "dolfin2vtk", "vtk_add_field", "vtk_write_file", "plot_vtk" ]

def dolfin2vtk(domain, V=None) :
    """
    FIXME:
    """
    
    if not hasvtk:
        raise RuntimeError("need vtk to use this function")

    mesh = domain.data()
    gdim = domain.geometric_dimension()
    mdim = domain.topological_dimension()

    # elaborate point coordinates
    if domain.coordinates() :
        phi = domain.coordinates()
        V = phi.function_space()
        order = V.ufl_element().degree()
        # first we need to reorganize the entries of phi as matrix
        # whose column are components.
        # we do this by tabulating dofs of each subspace of V
        idx = np.column_stack([ V.sub(i).dofmap().dofs()
                    for i in xrange(0, gdim) ])
        coords = phi.vector().array()[idx]
        # compute the connectivity from cell dofs
        reorder = np.array([[ 0, 1, 2 ],
                            [ 0, 1, 2, 5, 3, 4 ],
                            [ 0, 1, 2, 3, 9, 6, 8, 7, 5, 4 ]][mdim-1])
        if order == 1 :
            reorder = reorder[0:mdim+1]
        conn = np.array([ V.sub(0).dofmap().cell_dofs(c)[reorder]
                    for c in xrange(0, mesh.num_cells()) ])
        # now remap the ids with the new order
        perm = { v: k for k, v in enumerate(idx[:, 0]) }
        conn = np.vectorize(perm.get, otypes=[np.uintp])(conn)
    elif V is not None :
        # linear mesh, but dof coordinates from V
        order = V.ufl_element().degree()
        idx = V.sub(0).dofmap().dofs()
        coords = V.dofmap().tabulate_all_coordinates(mesh).reshape(-1, gdim)
        coords = coords[idx,:]
        # compute the connectivity from cell dofs
        reorder = np.array([[ 0, 1, 2 ],
                            [ 0, 1, 2, 5, 3, 4 ],
                            [ 0, 1, 2, 3, 9, 6, 8, 7, 5, 4 ]][mdim-1])
        if order == 1 :
            reorder = reorder[0:mdim+1]
        conn = np.array([ V.sub(0).dofmap().cell_dofs(c)[reorder]
                    for c in xrange(0, mesh.num_cells()) ])
        # now remap the ids with the new order
        perm = { v: k for k, v in enumerate(idx) }
        conn = np.vectorize(perm.get, otypes=[np.uintp])(conn)
    else :        
        order = 1
        # coordinates of the mesh
        coords = mesh.coordinates().copy()
        # connectivity
        conn = mesh.cells()

    coords = np.hstack([coords, np.zeros((coords.shape[0], 3-gdim))])

    # only these are supported by dolfin
    vtk_shape = { 1 : { 1 : VTK_LINE,
                        2 : VTK_TRIANGLE,
                        3 : VTK_TETRA },
                  2 : { 1 : VTK_QUADRATIC_EDGE,
                        2 : VTK_QUADRATIC_TRIANGLE,
                        3 : VTK_QUADRATIC_TETRA } }[order][mdim]

    # create the grid
    grid = vtkUnstructuredGrid()
    grid.SetPoints(array2vtkPoints(coords))
    grid.SetCells(vtk_shape, array2vtkCellArray(conn))

    return grid

def vtk_add_field(grid, fun) :
    if not hasvtk:
        raise RuntimeError("need vtk to use this function")

    V = fun.function_space()
    family = V.ufl_element().family()
    degree = V.ufl_element().degree()

    if fun.value_rank() > 0 :
        idx = np.column_stack([ V.sub(i).dofmap().dofs()
                    for i in xrange(0, V.num_sub_spaces()) ])
        fval = fun.vector().array()[idx]
    else :
        fval = fun.vector().array()

    if fun.name() == 'displacement' :
        # add zero columns if necessary
        gdim = V.num_sub_spaces()
        fval = np.hstack([fval, np.zeros((fval.shape[0], 3-gdim))])

    funvtk = array2vtk(fval)
    funvtk.SetName(fun.name())
    if family == 'Discontinuous Lagrange' and degree == 0 :
        grid.GetCellData().AddArray(funvtk)
    else :
        grid.GetPointData().AddArray(funvtk)

def vtk_write_file(grid, fname) :
    if not hasvtk:
        raise RuntimeError("need vtk to use this function")
    writer = vtkXMLUnstructuredGridWriter()
    writer.SetInputData(grid)
    writer.SetFileName(fname)
    writer.Write()

def plot_vtk(grid, ndiv = 0) :
    if not hasvtk:
        raise RuntimeError("need vtk to use this function")

    ugg = vtkUnstructuredGridGeometryFilter()
    ugg.SetInputData(grid)

    lingeo = vtkDataSetSurfaceFilter()
    lingeo.SetNonlinearSubdivisionLevel(ndiv)
    lingeo.SetInputConnection(ugg.GetOutputPort())
    lingeo.Update()

    mapper = vtkPolyDataMapper()
    mapper.SetInputConnection(lingeo.GetOutputPort())

    actor = vtkActor()
    actor.SetMapper(mapper)
    #actor.GetProperty().EdgeVisibilityOn()
    actor.GetProperty().SetRepresentationToWireframe()
    actor.GetProperty().SetColor(0.8,0.8,0.8)
    actor.GetProperty().SetEdgeColor(0.2,0.2,0.2)

    renderer = vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(1,1,1)
    renderer.SetUseDepthPeeling(1)
    renderer.SetMaximumNumberOfPeels(6)

    renWin = vtkRenderWindow()
    renWin.AddRenderer(renderer)
    renWin.SetAlphaBitPlanes(1)
    renWin.SetMultiSamples(0)

    renWinInteractor = vtkRenderWindowInteractor()
    renWinInteractor.SetRenderWindow(renWin)

    style = vtk.vtkInteractorStyleTrackballCamera()
    renWinInteractor.SetInteractorStyle(style)

    renWin.Render()
    renWin.SetSize(800, 800)
    renWinInteractor.Start()

