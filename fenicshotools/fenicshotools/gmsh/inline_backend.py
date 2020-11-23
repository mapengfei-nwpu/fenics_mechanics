"""FIXME:"""
from instant import inline
from textwrap import dedent
from os.path import splitext, expanduser, join
from os import getenv

GMSH_INCLUDE_DIRS = [ '/usr/include/gmsh', '/usr/local/include/gmsh',
                      join(expanduser('~'), 'local', 'include', 'gmsh')]
GMSH_LIBRARY_DIRS = [ '/usr/lib', '/usr/local/lib',
                      join(expanduser('~'), 'local', 'lib')]

if getenv('GMSH_DIR') :
    GMSH_INCLUDE_DIRS += [ join(getenv('GMSH_DIR'), 'include', 'gmsh') ]
    GMSH_LIBRARY_DIRS += [ join(getenv('GMSH_DIR'), 'lib') ]

def gmsh_cpp_inline(code,
                    extra_include_dirs=None,
                    extra_library_dirs=None,
                    extra_libraries=None,
                    extra_system_headers=None) :

    p = { 'system_headers' : [ 'Gmsh.h' ] + (extra_system_headers or []),
          'include_dirs' : GMSH_INCLUDE_DIRS + (extra_include_dirs or []),
          'library_dirs' : GMSH_LIBRARY_DIRS + (extra_library_dirs or []),
          'libraries' : [ 'Gmsh' ] + (extra_libraries or []) }

    return inline(code, **p)

def gmsh_cpp_geo2msh(geoname, dim=3, mshname=None, logname=None) :
    """
    Converts a GEO to a MSH directly using GMSH library.

    Note: GmshInitialize() possibly initializes MPI and PETSc if not
          already initialized, so be sure to do it beforehand.
    """

    mshname = mshname or splitext(geoname)[0] + '.msh'
    logname = logname or splitext(geoname)[0] + '.log'

    code = dedent(\
    """\
    void func(char geoname[], int dim, char mshname[], char logname[])
    {{
        fflush(stdout);
        int old = dup(1);
        int out = open(logname, O_CREAT | O_TRUNC | O_WRONLY, 0644);
        dup2(out, 1);
        close(out);
        GmshInitialize();
        GModel* model(new GModel);
        model->readGEO(geoname);
        model->mesh(dim);
        model->writeMSH(mshname);
        delete model;
        GmshFinalize();
        fflush(stdout);
        dup2(old, 1);
        close(old);
    }}\
    """)

    headers = [ 'GModel.h', 'fcntl.h', 'cstdio' ]
    func = gmsh_cpp_inline(code, extra_system_headers=headers)
    func(geoname, dim, mshname, logname)

