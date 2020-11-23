"""
FIXME: Not used anymore...
"""
from hashlib import sha1

def githash(data, chunk_size=8192):
    """
    Compute a hash of a binary file:
    FIXME:...
    """
    s = sha1()
    if isinstance(data, file) :
        fsize = os.fstat(data.fileno()).st_size
        s.update("blob {}\0".format(fsize))
        for chunk in iter(lambda : data.read(chunk_size), b'') :
            s.update(chunk)
    else :
        s.update("blob {}\0".format(long(len(data))))
        s.update(data)
    return s.hexdigest()


