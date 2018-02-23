"""Headers for the getent package."""

# pylint: disable=too-few-public-methods

from ctypes import ARRAY, POINTER
from ctypes import Structure, Union
from ctypes import c_int, c_long, c_ulong, c_void_p

from pyironbase.external.getent.constants import ctypes_c_char_p
from pyironbase.external.getent.constants import size_t, uint8_t, uint16_t, uint32_t

__all__ = (
    'AliasStruct',
    'GroupStruct',
    'HostStruct',
    'InAddr6Struct',
    'InAddr6Union',
    'InAddrStruct',
    'NetworkStruct',
    'PasswdStruct',
    'ProtoStruct',
    'RPCStruct',
    'ServiceStruct',
    'ShadowStruct',
)


class AliasStruct(Structure):

    """Struct `aliasent` from `<aliasent.h>`."""

    _fields_ = [
        ('name', ctypes_c_char_p),
        ('members_len', size_t),
        ('members', POINTER(ctypes_c_char_p)),
        ('local', c_int),
    ]


class InAddrStruct(Structure):

    """Struct `in_addr` from `<netinet/in.h>`."""

    _fields_ = [
        ('s_addr', uint32_t),
    ]


class InAddr6Union(Union):

    """Struct `in_addr6` from `<netinet6/in6.h>`."""

    _fields_ = [
        ('u6_addr8', ARRAY(uint8_t, 16)),
        ('u6_addr16', ARRAY(uint16_t, 8)),
        ('u6_addr32', ARRAY(uint32_t, 4)),
    ]


class InAddr6Struct(Structure):

    """Struct `in_addr6` from `<netinet6/in6.h>`."""

    _anonymous_ = ('in6_u',)
    _fields_ = [
        ('in6_u', InAddr6Union),
    ]


class HostStruct(Structure):

    """Struct `hostent` from `<netdb.h>`."""

    _fields_ = [
        ('name', ctypes_c_char_p),
        ('aliases', POINTER(ctypes_c_char_p)),
        ('addrtype', c_int),
        ('addr_list_len', c_int),
        ('addr_list', POINTER(c_void_p)),
    ]


class NetworkStruct(Structure):

    """Struct `netent` from `<netdb.h>`."""

    _fields_ = [
        ('name', ctypes_c_char_p),              # official network name
        ('aliases', POINTER(ctypes_c_char_p)),  # alias list
        ('addrtype', c_int),            # net address type
        ('net', uint32_t),              # network number
    ]


class ProtoStruct(Structure):

    """Struct `protoent` from `<netdb.h>`."""

    _fields_ = [
        ('name', ctypes_c_char_p),              # official protocol name
        ('aliases', POINTER(ctypes_c_char_p)),  # alias list
        ('proto', c_int),               # protocol number
    ]


class RPCStruct(Structure):

    """Struct `rpcent` from `<netdb.h>`."""

    _fields_ = [
        # Name of server for RPC program
        ('name', ctypes_c_char_p),
        # Alias list
        ('aliases', POINTER(ctypes_c_char_p)),
        # RPC program number
        ('number', c_long),
    ]


class ServiceStruct(Structure):

    """Struct `servent` from `<netdb.h>`."""

    _fields_ = [
        ('name', ctypes_c_char_p),              # official service name
        ('aliases', POINTER(ctypes_c_char_p)),  # alias list
        ('port', c_int),                # port number
        ('proto', ctypes_c_char_p),            # protocol to use
    ]


class GroupStruct(Structure):

    """Struct `group` from `<grp.h>`."""

    _fields_ = [
        ("name", ctypes_c_char_p),
        ("password", ctypes_c_char_p),
        ("gid", c_int),
        ("members", POINTER(ctypes_c_char_p)),
    ]


class PasswdStruct(Structure):

    """Struct `passwd` from `<pwd.h>`."""

    _fields_ = [
        ('name', ctypes_c_char_p),
        ('password', ctypes_c_char_p),
        ('uid', c_int),
        ('gid', c_int),
        ('gecos', ctypes_c_char_p),
        ('dir', ctypes_c_char_p),
        ('shell', ctypes_c_char_p),
    ]


class ShadowStruct(Structure):

    """Struct `spwd` from `<shadow.h>`."""

    _fields_ = [
        ('name', ctypes_c_char_p),
        ('password', ctypes_c_char_p),
        ('change', c_long),
        ('min', c_long),
        ('max', c_long),
        ('warn', c_long),
        ('inact', c_long),
        ('expire', c_long),
        ('flag', c_ulong),
    ]
