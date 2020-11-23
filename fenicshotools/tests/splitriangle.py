from fenics import *
from fenicshotools.gmsh import geo2dolfin
from fenicshotools.linearizedomain import *
import numpy as np

parameters['allow_extrapolation'] = True

if __name__ == "__main__" :
    code = """\
    Point(1) = { 0.0, 0.0, 0.0, 1.0 };
    Point(2) = { 1.0, 0.0, 0.0, 1.0 };
    Point(3) = { 0.0, 1.0, 0.0, 1.0 };
    Circle(1) = { 2, 1, 3 };
    Line(2) = { 3, 1 };
    Line(3) = { 1, 2 };
    Line Loop(1) = {1 , 2, 3 };
    Plane Surface(1) = { 1 };
    //Extrude {0.0, 0.0, 0.5} { Surface{1}; }
    Mesh.ElementOrder = 5;
    """
    domain, _ = geo2dolfin(code)

    if domain.coordinates() :
        degree = domain.coordinates().ufl_element().degree()
    else :
        degree = 1

    # a function
    V = FunctionSpace(domain, 'P', degree)
    u = interpolate(Expression("x[0]*x[0]+x[1]*x[1]", element=V.ufl_element()), V)

    refmesh, refu = linearize_domain_and_fields(domain, u, ndiv=0)
    plot(refu, interactive=True)
    #plot(refmesh, interactive=True)

