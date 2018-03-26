"""Constants and types for the getent package."""

# pylint: disable=invalid-name

import sys

from ctypes import c_uint, c_int
from ctypes import c_ubyte, c_ushort
from ctypes import c_char_p as ctypes_c_char_p
import socket

# Python 2 and 3 merry go round.
if sys.version_info[0] < 3:
    c_char_p = ctypes_c_char_p

else:
    def c_char_p(value):
        """C char pointer."""
        return value.encode('utf-8')

# Socket
AF_INET = socket.AF_INET
AF_INET6 = socket.AF_INET6
INADDRSZ = 4
IN6ADDRSZ = 16

# Types
uint8_t = c_ubyte
uint16_t = c_ushort
uint32_t = c_uint
uid_t = c_uint
gid_t = c_uint
size_t = c_int
sa_family_t = c_ushort
