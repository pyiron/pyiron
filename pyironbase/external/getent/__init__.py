"""Getent offers libc type lookup functions."""

# pylint: disable=too-few-public-methods

import socket
import struct
import sys
from ctypes import POINTER, cast, create_string_buffer, pointer
from datetime import datetime

from pyironbase.external.getent.constants import (
    AF_INET,
    AF_INET6,
    IN6ADDRSZ,
    INADDRSZ,
    c_char_p,
    gid_t,
    uid_t
)
from pyironbase.external.getent.libc import (
    getspnam,
    getpwnam,
    getpwuid,
    endaliasent,
    endgrent,
    endhostent,
    endnetent,
    endprotoent,
    endpwent,
    endrpcent,
    endservent,
    endspent,
    getaliasent,
    getgrent,
    getgrgid,
    getgrnam,
    getrpcbyname,
    getrpcbynumber,
    gethostbyaddr,
    gethostbyname2,
    gethostent,
    getnetent,
    getprotobyname,
    getprotobynumber,
    getprotoent,
    getpwent,
    getrpcent,
    getservent,
    getspent,
    inet_pton,
    setaliasent,
    setgrent,
    sethostent,
    setnetent,
    setprotoent,
    setpwent,
    setrpcent,
    setservent,
    setspent,
    getnetbyname,
    getservbyname,
    getservbyport
)

from pyironbase.external.getent import headers

__all__ = (
    'alias', 'group', 'host', 'network', 'passwd', 'proto', 'rpc', 'service',
    'shadow'
)

if sys.version_info[0] < 3:
    def convert23(value):
        """Python 2/3 compatibility function for emitting encoded strings."""
        return value

else:
    def convert23(value):
        """Python 2/3 compatibility function for emitting encoded strings."""
        return value.decode('utf-8') if isinstance(value, bytes) else value


class StructMap(object):

    """Base class for mapped C structs."""

    def __init__(self, p):
        """Generic struct constructor.

        Will also define getter and setter attributes on the instance for all
        struct members.
        """
        self.p = p

        for attr in dir(self.p.contents):

            if attr.startswith('_'):
                continue

            elif not hasattr(self, attr):
                value = getattr(self.p.contents, attr)
                setattr(self, attr, convert23(value))

    def __iter__(self):
        """Iterate over the mapped struct members.

        Yields `(key, value)` pairs.
        """
        for attr in dir(self.p.contents):
            if attr.startswith('_'):
                continue

            else:
                yield (attr, getattr(self, attr))

    def _map(self, attr):
        i = 0
        obj = getattr(self.p.contents, attr)

        while obj[i]:
            yield convert23(obj[i])
            i += 1


def _resolve(addrtype, addr):
    if addrtype == AF_INET:
        p = cast(addr, POINTER(headers.InAddrStruct))
        # Here we pack to little-endian to reverse back to
        # network byte order
        packed = struct.pack('<I', p.contents.s_addr)
        return socket.inet_ntop(addrtype, packed)
    elif addrtype == AF_INET6:
        p = cast(addr, POINTER(headers.InAddr6Struct))
        packed = ''.join([
            struct.pack('<L', bit)
            for bit in p.contents.in6_u.u6_addr32
        ])
        return socket.inet_ntop(addrtype, packed)


# pep257: disable=D102

class Host(StructMap):

    """Struct ``hostent`` from ``<netdb.h>``.

    .. py:attribute:: name

       Official name of the host.

    .. py:attribute:: aliases

       List of alternative names for the host.

    .. py:attribute:: addresses

       List of network addresses for the host.

    """

    addrtype = socket.AF_INET

    def __init__(self, p):
        super(Host, self).__init__(p)
        self.aliases = list(self._map('aliases'))
        self.addresses = [_resolve(self.addrtype, addr)
                          for addr in self._map('addr_list')]


class Proto(StructMap):

    """Struct ``protoent`` from ``<netdb.h>``.

    .. py:attribute:: name

       The official name of the protocol.

    .. py:attribute:: aliases

       List of alternative names for the protocol.

    .. py:attribute:: proto

       The protocol number.
    """

    def __init__(self, p):
        super(Proto, self).__init__(p)
        self.aliases = list(self._map('aliases'))


class RPC(StructMap):

    """Struct ``rpcent`` from ``<netdb.h>``.

    .. py:attribute:: name

       The name of the server for this RPC program.

    .. py:attribute:: aliases

       List of alternate names for the RPC program.

    .. py:attribute:: number

       The RPC program number for this service.
    """

    def __init__(self, p):
        super(RPC, self).__init__(p)
        self.aliases = list(self._map('aliases'))


class Service(StructMap):

    """Struct ``servent`` from ``<netdb.h>``.

    .. py:attribute:: name

       The official name of the service.

    .. py:attribute:: aliases

       List of alternative names for the service.

    .. py:attribute:: port

       The port number of the service.

    .. py:attribute:: proto

       The name of the protocol to use with this service.
    """

    def __init__(self, p):
        super(Service, self).__init__(p)
        self.aliases = list(self._map('aliases'))


class Network(StructMap):

    """Struct ``netent`` from ``<netdb.h>``.

    .. py:attribute:: name

       The official name of the network.

    .. py:attribute:: aliases

       List of alternative names for the network.

    .. py:attribute:: addrtype

       Type of the network number; always ``AF_INET``.

    .. py:attribute:: net

       The network number.
    """

    def __init__(self, p):
        super(Network, self).__init__(p)
        self.aliases = list(self._map('aliases'))


class Alias(StructMap):

    """Struct ``aliasent`` from ``<aliases.h>``.

    .. py:attribute:: name

       The name of the alias.

    .. py:attribute:: members

       List of members for the alias.

    .. py:attribute:: local

       Local flag.
    """

    def __init__(self, p):
        super(Alias, self).__init__(p)
        self.members = list(self._map('members'))


class Group(StructMap):

    """Struct ``grp`` from ``<grp.h>``.

    .. py:attribute:: name

       The name of the group.

    .. py:attribute:: passwd

       The password of the group.

    .. py:attribute:: gid

       The group ID number.

    .. py:attribute:: members

       List of group members.
    """

    def __init__(self, p):
        super(Group, self).__init__(p)
        self.members = list(self._map('members'))


class Passwd(StructMap):

    """Struct ``pwd`` from ``<pwd.h>``.

    .. py:attribute:: name

       The username.

    .. py:attribute:: passwd

       The password.

    .. py:attribute:: uid

       The user ID number.

    .. py:attribute:: gid

       The primary group ID number.

    .. py:attribute:: gecos

       User information.

    .. py:attribute:: dir

       Home directory.

    .. py:attribute:: shell

       Shell program.
    """

    pass


class Shadow(StructMap):

    """Struct ``spwd`` from ``<shadow.h>``.

    .. py:attribute:: name

       Login name.

    .. py:attribute:: pwd

       Encrypted password.

    .. py:attribute:: change

       Date of last change.

    .. py:attribute:: expire

       Date when account expires.

    .. py:attribute:: min

       Minimum number of days between changes.

    .. py:attribute:: max

       Maximum number of days between changes.

    .. py:attribute:: warn

       Number of days before password expires to warn the user to change it.

    .. py:attribute:: inact

       Number of days after password expires until account is disabled.
    """

    def __init__(self, p):
        super(Shadow, self).__init__(p)
        self.change = datetime.fromtimestamp(p.contents.change)
        self.expire = datetime.fromtimestamp(p.contents.expire)


def alias(search=None):
    """Perform a (mail) alias lookup.

    To iterate over all alias entries::

        >>> for item in alias():
        ...

    To lookup a single alias::

        >>> mail = alias('postmaster')
        >>> print mail.members
        root
    """
    if not callable(setaliasent):
        raise NotImplementedError

    if search is None:
        setaliasent()
        p = True
        r = []
        while p:
            p = getaliasent()
            if p:
                r.append(Alias(p))
        endaliasent()
        return r


def host(search=None):
    """Perform a host lookup.

    To iterate over all host entries::

        >>> for item in host():
        ...

    To lookup a single host by name::

        >>> import socket
        >>> server = host('localhost')
        >>> server.addrtype in [socket.AF_INET, socket.AF_INET6]
        True

    """
    if not callable(sethostent):
        raise NotImplementedError

    if search is None:
        sethostent()
        p = True
        r = []
        while p:
            p = gethostent()
            if p:
                r.append(Host(p))
        endhostent()
        return r

    else:
        def lookup():
            addr = create_string_buffer(IN6ADDRSZ)

            # Test if input is an IPv6 address
            if inet_pton(AF_INET6, c_char_p(search), pointer(addr)) > 0:
                host = gethostbyaddr(addr, IN6ADDRSZ, AF_INET6)
                return host

            # Test if input is an IPv4 address
            if inet_pton(AF_INET, c_char_p(search), pointer(addr)) > 0:
                host = gethostbyaddr(addr, INADDRSZ, AF_INET)
                return host

            # Test if input is a hostname with an IPv6 address
            host = gethostbyname2(c_char_p(search), socket.AF_INET6)
            if host:
                return host

            # Test if input is a hostname with an IPv4 address
            host = gethostbyname2(c_char_p(search), socket.AF_INET)
            if host:
                return host

        host = lookup()
        if bool(host):
            return Host(host)


def proto(search=None):
    """Perform a protocol lookup.

    To lookup all protocols::

        >>> for item in proto():
        ...

    To lookup a single protocol number::

        >>> tcp = proto('tcp')
        >>> print tcp.proto
        6
    """
    if not callable(setprotoent):
        raise NotImplementedError

    if search is None:
        setprotoent()
        prt = True
        res = []
        while prt:
            prt = getprotoent()
            if prt:
                res.append(Proto(prt))
        endprotoent()
        return res

    else:
        search = str(search)
        if search.isdigit():
            prt = getprotobynumber(uid_t(int(search)))
        else:
            prt = getprotobyname(c_char_p(search))

        if bool(prt):
            return Proto(prt)


def rpc(search=None):
    """Perform a remote procedure call lookup.

    To lookup all rpc services::

        >>> for item in rpc():
        ...

    To lookup one rpc service by name::

        >>> nfs = rpc('nfs')
        >>> print nfs.number
        100003
        >>> print nfs.aliases
        ['portmap', 'sunrpc']

    """
    if not callable(setrpcent):
        raise NotImplementedError

    if search is None:
        setrpcent()
        ent = True
        res = []
        while ent:
            ent = getrpcent()
            if ent:
                res.append(RPC(ent))
        endrpcent()
        return res

    else:
        search = str(search)
        if search.isdigit():
            ent = getrpcbynumber(uid_t(int(search)))
        else:
            ent = getrpcbyname(c_char_p(search))

        if bool(ent):
            return RPC(ent)


def service(search=None, protocol=None):
    """Perform a service lookup.

    To lookup all services::

        >>> for item in service():
        ...

    To lookup one service by port number::

        >>> http = service(0, 'tcp')
        >>> print http.port
        80

    Or by service name::

        >>> smtp = service('smtp', 'tcp')
        >>> print smtp.port
        25

    Or by short notation::

        >>> snmp = service('udp/snmp')
        >>> print snmp.port
        161

    """
    if search is None:
        setservent()
        srv = True
        res = []
        while srv:
            srv = getservent()
            if srv:
                res.append(Service(srv))
        endservent()
        return res

    else:
        search = str(search)
        if not protocol and '/' in search:
            protocol, search = search.split('/')
        if proto not in ['tcp', 'udp']:
            raise ValueError('Unsupported protocol "%s"' % (str(protocol),))
        if search.isdigit():
            srv = getservbyport(uid_t(int(search)), c_char_p(protocol))
        else:
            srv = getservbyname(c_char_p(search), c_char_p(protocol))

        if bool(srv):
            return Service(srv)


def network(search=None):
    """Perform a network lookup.

    To lookup all services::

        >>> for item in network():
        ...

    To lookup one network by name::

        >>> net = network('link-local')

    """
    if search is None:
        setnetent()
        net = True
        res = []
        while net:
            net = getnetent()
            if net:
                res.append(Network(net))
        endnetent()
        return res

    else:
        net = getnetbyname(c_char_p(search))
        if bool(net):
            return Network(net)


def group(search=None):
    """Perform a group lookup.

    To lookup all groups::

        >>> for item in group():
        ...

    To lookup one group by group id (gid)::

        >>> root = group(0)
        >>> print root.name
        'root'

    To lookup one group by name::

        >>> root = group('root')
        >>> print root.gid
        0

    """
    # Iterate over all group entries
    if search is None:
        setgrent()
        grp = True
        res = []
        while grp:
            grp = getgrent()
            if grp:
                res.append(Group(grp))
        endgrent()
        return res

    else:
        search = str(search)
        if search.isdigit():
            grp = getgrgid(gid_t(int(search)))
        else:
            grp = getgrnam(c_char_p(search))

        if bool(grp):
            return Group(grp)


def passwd(search=None):
    """Perform a passwd lookup.

    To lookup all passwd entries::

        >>> for item in passwd():
        ...

    To lookup one user by user id (uid)::

        >>> root = passwd(0)
        >>> print root.name
        'root'

    To lookup one user by name::

        >>> root = passwd('root')
        >>> print root.uid
        0

    """
    if not callable(setpwent):
        raise NotImplementedError

    # Iterate over all passwd entries
    if search is None:
        setpwent()
        pwd = True
        res = []
        while pwd:
            pwd = getpwent()
            if pwd:
                res.append(Passwd(pwd))
        endpwent()
        return res

    else:
        search = str(search)
        if search.isdigit():
            pwd = getpwuid(uid_t(int(search)))
        else:
            pwd = getpwnam(c_char_p(search))

        if bool(pwd):
            return Passwd(pwd)


def shadow(search=None):
    """Perform a shadow lookup.

    To lookup all shadow entries::

        >>> for item in shadow():
        ...

    To lookup one user by name::

        >>> root = shadow('root')
        >>> print root.warn # doctest: +SKIP
        99999

    """
    # Iterate over all shadow entries
    if search is None:
        if not callable(setspent):
            raise NotImplementedError('Not available on %s' % (sys.platform,))
        setspent()
        spe = True
        res = []
        while spe:
            spe = getspent()
            if spe:
                res.append(Shadow(spe))
        endspent()
        return res

    else:
        if not callable(getspnam):
            raise NotImplementedError

        spe = getspnam(search)
        if bool(spe):
            return Shadow(spe)


if __name__ == '__main__':
    # pylint: disable=superfluous-parens
    print(dict(host('127.0.0.1')))
    print(dict(host('localhost')))

    for h in host():
        print(dict(h))

    print(dict(network('link-local') or {}))

    for p in proto():
        print(p.name, p.aliases)

    for r in rpc():
        print(r.name, r.aliases)

    for s in service():
        print(s.name, s.port, s.proto)

    for n in network():
        print(n.name, n.aliases)

    for g in group():
        print(g.name, g.members, dict(g))

    for p in passwd():
        print(p.name, p.shell, p.dir, dict(p))

    try:
        for s in shadow():
            print(s.name, s.change, s.expire, dict(s))
    except NotImplementedError as error:
        print('shadow() failed:', error)
