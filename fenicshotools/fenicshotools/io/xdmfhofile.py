from textwrap import dedent
from dolfin import *
import numpy as np
from six import StringIO

__all__ = ["XdmfHoFile"]

class XdmfHoFile(object) :
    """
    File(ish) object for saving higher order geometry to XDMF, such
    that it can be visualized with for example ParaView.
    """
    def __init__(self, ofilename):
        """
        Create a File like object that stores Domains and Functions to a XDMF file.
        
        *Arguments*
          ofilename (str)
            The name of the xdmf file
        """
        self._fname = ofilename
        self._grids = []

    def add(self, grid_name, data, component=None):
        """
        Add a field or mesh to be written to the XDMF file.

        *Arguments*
          grid_name (str)
            The name of the grid stored in the XDMF file
          data (:py:class:`ufl.Domain` or :py:class:`ufl.Function`)
            The higher order data that should be saved to file
          component (int)
            Optional argument for saving a sub-component in the given data.
        """

        if isinstance(data, Domain):
            dom = data
            saving_function = False
            
        elif isinstance(data, Function):
            dom = data.domain()
            saving_function = True
            if component is not None and not (\
                isinstance(component, int) and 0 <= component < \
                data.function_space().num_sub_spaces()):
                raise ValueError("expected component to be an int in range "\
                                 "[0,{})".format(data.function_space().num_sub_spaces()))
        else:
            raise TypeError("expected either Domain or Function as the second argument")
        
        phi = dom.coordinates()
        if phi is None:
            raise ValueError("Got linear domain, expected a higher order.")

        Vphi = phi.function_space()
        mesh = dom.data()
        mdim = dom.topological_dimension()
        gdim = dom.geometric_dimension()
        morder = Vphi.ufl_element().degree()

        # first we need to reorganize the entries of phi as matrix
        # whose column are components.
        # we do this by tabulating dofs of each subspace of Vphi
        idx = np.column_stack([ Vphi.sub(i).dofmap().dofs() 
                    for i in xrange(0, gdim) ])
        coords = phi.vector().array()[idx]

        # always 3d coordinates
        coords = np.column_stack([coords,
                                  np.zeros((coords.shape[0], 3-gdim))])

        # compute the connectivity from cell dofs
        reorder = np.array([[ 0, 1, 2 ],
                            [ 0, 1, 2, 5, 3, 4 ],
                            [ 0, 1, 2, 3, 9, 6, 8, 7, 5, 4 ]][mdim-1])
        if morder == 1 :
            reorder = reorder[0:mdim+1]
        conn = np.array([ Vphi.sub(0).dofmap().cell_dofs(c)[reorder]
                    for c in xrange(0, mesh.num_cells()) ])

        # now remap the ids with the new order
        perm = { v: k for k, v in enumerate(idx[:, 0]) }
        conn = np.vectorize(perm.get, otypes=[np.uintp])(conn)

        topo = self._generate_topology(conn, mdim, morder)
        geom = self._generate_geometry(coords)

        # If saving Function as an attribut
        attr = ""
        if saving_function:
            if component:
                Vval = data.function_space().sub(component)
                data = data.sub(component, True)
            else :
                Vval = data.function_space()
                
            vals = data.vector().array()
            rank = data.rank()

            if rank > 0 :
                val_idx = np.column_stack([ Vval.sub(i).dofmap().dofs() 
                            for i in xrange(0, mdim) ])
                vals = vals[idx]
                # always 3d coordinates
                vals = np.column_stack([vals,
                                        np.zeros((vals.shape[0], 3-gdim))])

            attr = self._generate_attribute(data.name(), vals, data.rank())

        self._grids.append({ "grid_name" : grid_name,
                             "grid_content" : topo + geom + attr })

    def _dataitem_dump(self, data, fmt) :

        if isinstance(data, str) :
            # assuming HDF5 reference
            return "HDF", data
        elif isinstance(data, np.ndarray) :
            # dumping data into str
            data_str = StringIO.StringIO()
            np.savetxt(data_str, data, fmt = fmt)
            return "XML", data_str.getvalue()
        else :
            raise ValueError

    def _generate_geometry(self, coords) :

        geo_datatype, geo_data = self._dataitem_dump(coords, '%f')

        geometry = dedent(\
        """
        <Geometry GeometryType="XYZ">
          <DataItem Format="{geo_datatype}" DataType="Float" Dimensions="{num_points} 3">
{coords}
          </DataItem>
        </Geometry>""").format(geo_datatype = geo_datatype,
                    num_points = coords.shape[0],
                    coords = geo_data)

        return geometry

    def _generate_topology(self, conn, mdim, order) :

        vtk_types = [[ "Polyline", "Triangle", "Tetrahedron" ],
                     [ "Edge_3", "Tri_6", "Tet_10" ]]
        cell_dofs = [[ 1, 3, 4 ],
                     [ 3, 6, 10 ]]

        conn_datatype, conn_data = self._dataitem_dump(conn, '%d')

        topology = dedent(\
        """
        <Topology NumberOfElements="{num_cells}" TopologyType="{cell_type}">
          <DataItem Format="{data_type}" DataType="Int" Dimensions="{num_cells} {cell_dofs}">
{connectivity}
          </DataItem>
        </Topology>""").format(num_cells = conn.shape[0],
                    cell_type = vtk_types[order-1][mdim-1],
                    data_type = conn_datatype,
                    cell_dofs = cell_dofs[order-1][mdim-1],
                    connectivity = conn_data)

        return topology

    def _generate_attribute(self, name, values,
                            rank, center = "Node", symmetry = False) :

        vtk_shapes  = [ "Scalar", "Vector", "Tensor" ]

        comp = pow(3, rank)
        
        if symmetry and rank == 2:
            comp /= 2
            vtk_shapes[2] += "6"

        if center not in [ "Node", "Cell" ] :
            raise ValueError

        attr_datatype, attr_data = self._dataitem_dump(values, '%f')
        
        attribute = dedent(\
        """
        <Attribute Name=\"{name}\" AttributeType=\"{shape}\" Center=\"{center}\">
          <DataItem Format=\"{datatype}\" DataType=\"Float\" Dimensions=\"{num_vals} {comp}\">
{attribute_data}</DataItem>
        </Attribute>""").format(name = name,
                    shape = vtk_shapes[rank],
                    center = center,
                    datatype = attr_datatype,
                    num_vals = values.shape[0],
                    comp = comp,
                    attribute_data = attr_data)

        return attribute

    def write_to_file(self) :

        grid_form = """
      <Grid Name=\"{grid_name}\" GridType="Uniform">{grid_content}
      </Grid>"""
        time_form = """
      <Time TimeType=\"List\">
        <DataItem Format=\"XML\" Dimensions=\"{}\"> {} </DataItem>
      </Time>"""

        body = """\
<?xml version=\"1.0\"?>
<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>
<Xdmf Version=\"2.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\">
  <Domain>
    <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">{}{}
    </Grid>
  </Domain>
</Xdmf>
""".format(time_form.format(len(self._grids), " ".join(\
      str(i+1) for i in range(len(self._grids)))),
      "".join(grid_form.format(**grid) for grid in self._grids))
        
        with open(self._fname, 'w') as f :
            f.write(body)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.write_to_file()
        
