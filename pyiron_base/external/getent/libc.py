"""Libc functions for the getent package."""

# pylint: disable=invalid-name

from ctypes import CDLL, POINTER
from ctypes import c_int, c_void_p, c_uint
from ctypes.util import find_library

from pyironbase.external.getent.constants import ctypes_c_char_p, gid_t, uid_t

from pyironbase.external.getent import headers

__all__ = (
    'endgrent',
    'endhostent',
    'endnetent',
    'endprotoent',
    'endpwent',
    'endrpcent',
    'endservent',
    'endspent',
    'getgrent',
    'getgrgid',
    'getgrnam',
    'gethostbyaddr',
    'gethostbyname2',
    'gethostent',
    'getnetbyname',
    'getnetent',
    'getprotobyname',
    'getprotobynumber',
    'getprotoent',
    'getpwent',
    'getpwnam',
    'getpwuid',
    'getrpcbyname',
    'getrpcbynumber',
    'getrpcent',
    'getservbyname',
    'getservbyport',
    'getservent',
    'getspent',
    'inet_pton',
    'libc',
    'setgrent',
    'sethostent',
    'setnetent',
    'setprotoent',
    'setpwent',
    'setrpcent',
    'setservent',
    'setspent',
)

#: libc object
libc = CDLL(find_library("c"))

# Map libc function calls

if hasattr(libc, 'getaliasent'):
    endaliasent = libc.endaliasent
    getaliasent = libc.getaliasent
    getaliasent.restype = POINTER(headers.AliasStruct)
    setaliasent = libc.setaliasent
else:
    endaliasent = None
    getaliasent = None
    setaliasent = None

if hasattr(libc, 'gethostent'):
    endhostent = libc.endhostent
    gethostent = libc.gethostent
    gethostent.restype = POINTER(headers.HostStruct)
    sethostent = libc.sethostent
else:
    endhostent = None
    gethostent = None
    sethostent = None

if hasattr(libc, 'getnetent'):
    endnetent = libc.endnetent
    getnetent = libc.getnetent
    getnetent.restype = POINTER(headers.NetworkStruct)
    setnetent = libc.setnetent
else:
    endnetent = None
    getnetent = None
    setnetent = None

if hasattr(libc, 'getprotoent'):
    endprotoent = libc.endprotoent
    getprotoent = libc.getprotoent
    getprotoent.restype = POINTER(headers.ProtoStruct)
    setprotoent = libc.setprotoent
else:
    endprotoent = None
    getprotoent = None
    setprotoent = None

if hasattr(libc, 'getrpcent'):
    endrpcent = libc.endrpcent
    getrpcent = libc.getrpcent
    getrpcent.restype = POINTER(headers.RPCStruct)
    setrpcent = libc.setrpcent
else:
    endrpcent = None
    getrpcent = None
    setrpcent = None

if hasattr(libc, 'getservent'):
    endservent = libc.endservent
    getservent = libc.getservent
    getservent.restype = POINTER(headers.ServiceStruct)
    setservent = libc.setservent
else:
    endservent = None
    getservent = None
    setservent = None

if hasattr(libc, 'getgrent'):
    endgrent = libc.endgrent
    getgrent = libc.getgrent
    getgrent.restype = POINTER(headers.GroupStruct)
    setgrent = libc.setgrent
else:
    endgrent = None
    getgrent = None
    setgrent = None

if hasattr(libc, 'getpwent'):
    endpwent = libc.endpwent
    getpwent = libc.getpwent
    getpwent.restype = POINTER(headers.PasswdStruct)
    setpwent = libc.setpwent
else:
    endpwent = None
    getpwent = None
    setpwent = None

if hasattr(libc, 'getspent'):
    endspent = libc.endspent
    getspent = libc.getspent
    getspent.restype = POINTER(headers.ShadowStruct)
    setspent = libc.setspent
else:
    endspent = None
    getspent = None
    setspent = None

#: inet_pton(family, src, dst)
inet_pton = libc.inet_pton
inet_pton.argtypes = (c_int, ctypes_c_char_p, c_void_p)
inet_pton.restype = c_int

#: gethostbyaddr(addr)
gethostbyaddr = libc.gethostbyaddr
gethostbyaddr.argtypes = (c_void_p, )
gethostbyaddr.restype = POINTER(headers.HostStruct)

#: gethostbyname2(name, family)
gethostbyname2 = libc.gethostbyname2
gethostbyname2.argtypes = (ctypes_c_char_p, c_uint)
gethostbyname2.restype = POINTER(headers.HostStruct)

#: getnetbyname(name)
getnetbyname = libc.getnetbyname
getnetbyname.argtypes = (ctypes_c_char_p,)
getnetbyname.restype = POINTER(headers.NetworkStruct)

#: getgrnam(name)
getgrnam = libc.getgrnam
getgrnam.argtypes = (ctypes_c_char_p,)
getgrnam.restype = POINTER(headers.GroupStruct)

#: getgrgid(gid)
getgrgid = libc.getgrgid
getgrgid.argtypes = (gid_t,)
getgrgid.restype = POINTER(headers.GroupStruct)

#: getpwnam(name)
getpwnam = libc.getpwnam
getpwnam.argtypes = (ctypes_c_char_p,)
getpwnam.restype = POINTER(headers.PasswdStruct)

#: getpwuid(uid)
getpwuid = libc.getpwuid
getpwuid.argtypes = (uid_t,)
getpwuid.restype = POINTER(headers.PasswdStruct)

#: getprotobyname(name)
getprotobyname = libc.getprotobyname
getprotobyname.argtypes = (ctypes_c_char_p,)
getprotobyname.restype = POINTER(headers.ProtoStruct)

#: getprotobynumber(n)
getprotobynumber = libc.getprotobynumber
getprotobynumber.argtypes = (c_int,)
getprotobynumber.restype = POINTER(headers.ProtoStruct)

#: getrpcbyname(name)
getrpcbyname = libc.getrpcbyname
getrpcbyname.argtypes = (ctypes_c_char_p,)
getrpcbyname.restype = POINTER(headers.RPCStruct)

#: getrpcbynumber(n)
getrpcbynumber = libc.getrpcbynumber
getrpcbynumber.argtypes = (c_int,)
getrpcbynumber.restype = POINTER(headers.RPCStruct)

#: getservbyname(name)
getservbyname = libc.getservbyname
getservbyname.argtypes = (ctypes_c_char_p, ctypes_c_char_p)
getservbyname.restype = POINTER(headers.ServiceStruct)

#: getservbyport(port)
getservbyport = libc.getservbyport
getservbyport.argtypes = (c_int, ctypes_c_char_p)
getservbyport.restype = POINTER(headers.ServiceStruct)

# Not supported on all platforms
try:
    getspnam = libc.getspnam
    getspnam.argtypes = (ctypes_c_char_p,)
    getspnam.restype = POINTER(headers.ShadowStruct)
except AttributeError:
    getspnam = None
