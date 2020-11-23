"""FIXME:"""
# Gmsh topology of reference element
import six

__all__ = ['gmsh_shape', 'gmsh_shape_entity', 'msh2ufc_shape']

class GmshEntity(object) :
    def __init__(self, topology) :
        self._topology = topology
        self._dim = max(topology.keys())
    def dimension(self) :
        return self._dim
    def num_vertices(self) :
        return len(self._topology[self._dim][0])
    def topology(self) :
        return self._topology

# FIXME Check quad and hexa facets!
gmsh_shape = {
    'point' :
        GmshEntity({ 0 : [ (0, ) ] }),
    'interval' :
        GmshEntity({ 1 : [ (0, 1) ] }),
    'triangle' :
        GmshEntity({ 1 : [ (0, 1), (1, 2), (2, 0) ],
                     2 : [ (0, 1, 2) ] }),
    'quadrilateral' :
        GmshEntity({ 1 : [ (0, 1), (1, 2), (2, 3), (3, 0) ],
                     2 : [ (0, 1, 2, 3) ] }),
    'tetrahedron' :
        GmshEntity({ 1 : [ (0, 1), (1, 2), (0, 2), (0, 3), (2, 3), (1, 3) ],
                     2 : [ (0, 1, 3), (0, 2, 1), (0, 3, 2), (1, 2, 3) ],
                     3 : [ (0, 1, 2, 3) ] }),
    'hexahedron' :
        GmshEntity({ 1 : [ (0, 1), (0, 3), (0, 4), (1, 2), (1, 5), (2, 3),
                           (2, 6), (3, 7), (4, 5), (4, 7), (5, 6), (6, 7) ],
                     2 : [ (0, 1, 2, 3), (0, 1, 5, 4), (0, 3, 7, 4),
                           (1, 2, 6, 5), (3, 2, 6, 7), (4, 7, 6, 5) ],
                     3 : [ (0, 1, 2, 3, 4, 5, 6, 7) ] }) }

# UFC shape to gmsh shapes. The position in the list
# reveals the order of the element.
ufc2msh_shape = {
    # 0d entities
    'point' : [ 15 ],
    # 1d entities
    'interval' : [ 1, 8, 26, 27, 28, 62 ],
    # 2d entities
    'triangle' : [ 2, 9, 21, 23, 25, 42 ],
    'quadrilateral' : [ 3, 10 ],
    # 3d entities
    'tetrahedron' : [ 4, 11, 29, 30, 31 ],
    'hexahedron'  : [ 5, 12, 92, 93 ] }

# gmsh to (UFC shape, order), computed from the previous one.
msh2ufc_shape = {gtype:(shape,order)
        for shape,gtypes in six.iteritems(ufc2msh_shape)
        for order,gtype in enumerate(gtypes,1) }

# entity to facet
gmsh_shape_entity = {
    "point" : [],
    "interval" : [ "point" ],
    "triangle" : [ "point", "interval" ],
    "quadrilateral" : [ "point", "interval" ],
    "tetrahedron" : [ "point", "interval", "triangle" ],
    "hexahedron" : [ "point", "interval", "quadrilateral" ] }

