#!/usr/bin/env python
# -*- coding: utf-8 -*-

# .. _generate_from_geo:
#
# .. py:currentmodule:: fenicshotools
#

# Mesh generation from geo file strings
# =====================================
#
# This demo shows how one can generate dolfin meshes on the fly using
# gmsh geo file.

# Implementation
# **************

# We begin with importing the fenics and geo2dolfin function

from fenics import *
from fenicshotools import geo2dolfin
from textwrap import dedent

def square_example():
    "Create a square mesh"

    code = dedent(\
    """\
    Point(1) = { 0.0, 0.0, 0.0, 0.1 };
    Point(2) = { 1.0, 0.0, 0.0, 0.1 };
    Point(3) = { 1.0, 1.0, 0.0, 0.1 };
    Point(4) = { 0.0, 1.0, 0.0, 0.1 };
    Line(1) = { 1, 2 };
    Line(2) = { 2, 3 };
    Line(3) = { 3, 4 };
    Line(4) = { 4, 1 };

    Line Loop(1) = { 1, 2, 3, 4 };
    Plane Surface(1) = { 1 };

    Physical Line("bottom") = { 1 };
    Physical Point("origin") = { 1 };
    Physical Surface("tissue") = { 1 };

    Mesh.ElementOrder = 1;
    """)
    return code

def disk_example():
    "Create a disk mesh"
    code = dedent(\
    """\
    Rin   = 0.5;
    Rout  = 2.0;
    psize = 1.5;

    Geometry.CopyMeshingMethod = 1;
    Mesh.ElementOrder = 2;
    Mesh.Optimize = 1;
    Mesh.OptimizeNetgen = 1;
    Mesh.HighOrderOptimize = 1;

    Point(1) = { Rin,  0.0, 0.0, psize };
    Point(2) = { Rout, 0.0, 0.0, psize };

    Line(1) = { 1, 2 };

    line_id = 1;
    For num In {0:2}
        out[] = Extrude
            { { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 0.0 }, Pi/2 }
            { Line{line_id}; };
        line_id = out[0];
        surf[num]    = out[1];
        bc_epi[num]  = out[2];
        bc_endo[num] = out[3];
    EndFor
    out[] = Extrude
        { { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 0.0 }, Pi/2 }
        { Line{line_id}; };
    surf[num]    = out[1];
    bc_epi[num]  = out[2];
    bc_endo[num] = out[3];

    Physical Surface("Myocardium") = { surf[] };
    Physical Line("Epicardium")    = { bc_epi[] };
    Physical Line("Endocardium")   = { bc_endo[] };
    """)
    return code

def prolate_example() :
    "Create a prolate mesh"
    code = dedent(\
    """
    Geometry.CopyMeshingMethod = 1;
    Mesh.ElementOrder = 2;
    Mesh.Optimize = 1;
    Mesh.OptimizeNetgen = 1;
    Mesh.HighOrderOptimize = 1;

    d_focal = 3.7;
    l_epi   = 0.7;
    l_endo  = 0.4;
    mu_base = 120.0 / 180.0 * Pi;

    Function EllipsoidPoint
        Point(id) = { r1 * Sin(mu) * Cos(theta),
                      r1 * Sin(mu) * Sin(theta),
                      r2 * Cos(mu), 0.5 };
    Return

    center = newp; Point(center) = { 0.0, 0.0, 0.0 };

    theta = 0.0;

    r1 = d_focal * Sinh(l_endo);
    r2 = d_focal * Cosh(l_endo);
    mu = 0.0;
    apex_endo = newp; id = apex_endo; Call EllipsoidPoint;
    mu = mu_base;
    base_endo = newp; id = base_endo; Call EllipsoidPoint;

    r1 = d_focal * Sinh(l_epi);
    r2 = d_focal * Cosh(l_epi);
    mu = 0.0;
    apex_epi = newp; id = apex_epi; Call EllipsoidPoint;
    mu = Acos(Cosh(l_endo) / Cosh(l_epi) * Cos(mu_base));
    base_epi = newp; id = base_epi; Call EllipsoidPoint;

    apex = newl; Line(apex) = { apex_endo, apex_epi };
    base = newl; Line(base) = { base_endo, base_epi };
    endo = newl; Ellipse(endo) = { apex_endo, center, apex_endo, base_endo };
    epi  = newl; Ellipse(epi) = { apex_epi, center, apex_epi, base_epi };

    ll1 = newll; Line Loop(ll1) = { apex, epi, -base, -endo };
    s1 = news; Plane Surface(s1) = { ll1 };

    out[] = Extrude { { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 0.0 }, Pi/2 }
                    { Surface{s1}; };
    out[] = Rotate { { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 0.0 }, Pi/2 }
                   { Duplicata{Volume{1};} };
    out[] = Rotate { { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 0.0 }, Pi/2 }
                   { Duplicata{Volume{out[0]};} };
    out[] = Rotate { { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 0.0 }, Pi/2 }
                   { Duplicata{Volume{out[0]};} };
    """)
    return code

if __name__ == "__main__" :

    # linear mesh
    mesh,phi,markers = geo2dolfin(square_example())

    # quadratic mesh
    mesh,phi,markers = geo2dolfin(disk_example())

    # 3d quadratic mesh
    mesh,phi,markers = geo2dolfin(prolate_example())
