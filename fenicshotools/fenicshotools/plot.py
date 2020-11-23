"""Module for ploting potentially higher order domains using vtk tools
or fenics"""
from fenicshotools.vtkutils import *

__all__ = ["plot_geometry"]

def plot_geometry(domain):
    """
    Plot a potentially higher order Domain using either vtk or dolfin

    *Arguments*
      domain (:py:class:`ufl.Domain`)
        The domain being plotted
    """
    if hasvtk :
        grid = dolfin2vtk(domain)
        plot_vtk(grid, 4)
    else :
        plot(domain.data(), interactive=True)

