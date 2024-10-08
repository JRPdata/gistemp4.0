#!/usr/local/bin/python3.4
#
# fetch.py
#
# David Jones and Nick Barnes, Climate Code Foundation.
# Copyright (C) 2008-2010 Ravenbrook Limited.
# Copyright (C) 2011 Climate Code Foundation.
# Avi Persin, Revision 2017-07-28

"""
fetch.py [--help] [--list] [--force] [--store <dir>] [--config <file>] [pattern] ...

Script to fetch (download from the internet) the inputs required for
the cccgistemp program.

The groups, bundles, bundle members, and individual files fetchable
are defined by the configuration file specified by --config <file>
(default 'config/sources').

Everything fetched ends up in the directory specified by --store <dir>
(default 'input').

If no such arguments are given, the default is to fetch those files in
the default group.  In the config file provided with cccgistemp, that
means the files required for normal cccgistemp operation.

Any arguments are treated as regular expressions and matched against
groups, bundles, files, or bundle members (in that order).  The first
matching item for each argument is fetched.

Unless --force is set, no file that already exists is created.

--list lists all things that can be fetched.

The config file syntax is as follows:

    Comments begin '#' and run to the end of the line.
    Every line begins with a keyword, followed by a colon.
    Keywords are not case-sensitive.

    A fetchable item is either:

        file: <url> [<local filename>]

    which denotes a fetchable item which is also a source dataset, or

        bundle: <url> [<local filename>]

    which denotes a 'bundle': a file which can be unpacked into
    a number of files, one or more of which may be source datasets,
    identified thusly:

        member: <pattern> [<local filename>]

    <pattern> is a regular expression matching the tail of a pathname
    within the most-recently described bundle.

    The system works out for itself how to unpack a bundle.  These may
    be based on the bundle's name or contents: you shouldn't have to
    worry about it.

    <local filename> in each of the above is an optional name to give
    the fetched item or extracted member.  If absent, the system uses
    a filename derived from the fetched item or extracted member.

    <url> may be any ftp:// or http:// URL.  It may also be of this
    form:

       ftpmatch://<site>/<path>/<pattern>

    In which case the directory <path> on the FTP site <site> is
    searched for filenames matching <pattern> and the last such file
    is fetched.  This 'feature' was developed for USHCN version 2
    datasets.

    All the contents of this file may be divided into disjoint
    'groups'.  Each group is named.  Groups are introduced with group
    lines:

      group: <group name>

    The default group name is the empty string.
"""

# http://www.python.org/doc/2.4.4/lib/module-getopt.html
import getopt
# http://docs.python.org/release/2.4.4/lib/module-os.html
# http://www.python.org/doc/2.4.4/lib/module-sys.html
import sys
# https://docs.python.org/2.6/library/urllib2.html
import urllib.request
import ssl

import itertools
import re

from settings import *

import tarfile
import zipfile


class Fetcher(object):
    def __init__(self, **kwargs):

        self.force = kwargs.pop('force', False)
        self.output = kwargs.pop('output', sys.stdout)
        self.output.write("Fetching Input Files:\n")
        self.prefix = kwargs.pop('prefix', INPUT_DIR)
        self.config_file = kwargs.pop('config_file', SOURCES_DIR + 'sources.txt')
        self.requests = kwargs.pop('requests', None)

    def fetch(self):
        (bundles, files) = self.find_requests(self.requests)
        for url, local in files:
            # first, check if ghcn file exists
            if "ghcnm.tavg.qcf.dat" in url:
                if not os.path.exists(INPUT_DIR + "ghcnm.tavg.qcf.dat"):
                    self.get_ghcn_file(url)
            else:
                self.fetch_one(url, local)
        for ((url, local), members) in bundles.items():
            self.fetch_one(url, local, members=members)
        sys.stdout.flush()

    # Get most recent GHCN data file
    def get_ghcn_file(self, url):
        public_dir = url.replace("ghcnm.tavg.qcf.dat", "")
        import datetime
        import urllib.request as urllib
        response = urllib.urlopen(public_dir)
        html = response.read().decode()
        filenames = re.findall(r'href=[\'"]?([^\'" >]+)', html)
        ghcn_filenames = [x for x in filenames if x[-4:] == ".dat" and "v4" in x]
        last_modifieds = []
        for file in ghcn_filenames:
            conn = urllib.urlopen(public_dir + file)
            last_modified = conn.headers["last-modified"]
            last_modified = datetime.datetime.strptime(last_modified[5:-4], '%d %b %Y %H:%M:%S')
            last_modifieds.append((file, last_modified))
        last_modifieds.sort(key=lambda x: x[1])
        recent_ghcn = [x[0] for x in last_modifieds[-2:]]
        if "qcf" in recent_ghcn[0]:
            qcf_file = recent_ghcn[0]
        else:
            qcf_file = recent_ghcn[1]
        print(qcf_file)
        urllib.urlretrieve(public_dir + qcf_file, INPUT_DIR + "ghcnm.tavg.qcf.dat")

    def make_prefix(self):
        try:
            os.makedirs(self.prefix)
        except OSError:
            # Expected if the directories already exist.
            pass

    def key_lines(self):
        comment_re = re.compile(r'((.*?[^\\])??)#')
        key_re = re.compile(r'^([a-zA-Z_]+)\s*:\s*(.*)$')
        for (no, l) in zip(itertools.count(1), open(self.config_file)):
            m = comment_re.match(l)
            if m:
                bare = m.group(1)
            else:
                bare = l
            bare = bare.strip()
            # ignore blank lines
            if len(bare) == 0:
                continue
            m = key_re.match(bare)
            if m:
                yield (no, m.groups())
            else:
                raise Error("%s:%d: malformed line '%s'" % (self.config_file, no, l.strip()))

    def read_config(self):
        valid_keys = dict(group=re.compile(r'^\s*(.*?)\s*$'),
                          file=re.compile(r'^([^\s]+)(\s+.*)?\s*$'),
                          bundle=re.compile(r'^([^\s]+)(\s+.*)?\s*$'),
                          member=re.compile(r'^([^\s]+)(\s+.*)?\s*$'))
        group = ''
        config = {'': dict(files=[], bundles={})}
        for (no, (k, v)) in self.key_lines():
            k = k.lower()
            if k not in valid_keys:
                raise Error("%s:%d: unknown key '%s'" % (self.config_file, no, k))
            m = valid_keys[k].match(v)
            if not m:
                raise Error("%s:%d: malformed '%s' line" % (self.config_file, no, k))

            # 'bundle' only persists over 'member' lines.
            if k != 'member':
                bundle = None

            if k == 'group':
                group = m.group(1)
                config[group] = dict(files=[], bundles={})
            elif k == 'file':
                config[group]['files'].append(m.groups())
                pattern = m.group(1)
                local = m.group(2)
            elif k == 'bundle':
                bundle = m.groups()
                members = []
                config[group]['bundles'][bundle] = members
                pattern = m.group(1)
                local = m.group(2)
            elif k == 'member':
                if bundle is None:
                    raise Error("%s:%d: 'member' line with no bundle." % (self.config_file, no))
                config[group]['bundles'][bundle].append(m.groups())
                pattern = m.group(1)
                local = m.group(2)

        return config

    def list_things(self):
        """List the things that we know how to fetch."""

        config = self.read_config()
        group_names = config.keys()
        group_names = sorted(group_names)
        for g in group_names:
            if g == '':
                self.output.write("Default group: \n")
            else:
                self.output.write("Group '%s':\n" % g)
            bs = config[g]['bundles'].items()
            bs = sorted(bs)
            for ((pattern, local), members) in bs:
                self.output.write("  bundle '%s':\n" % pattern)
                if local:
                    self.output.write("   (read to '%s')\n" % local)
                for (p, l) in members:
                    self.output.write("    member '%s'\n" % p)
                    if l:
                        self.output.write("    (read to '%s')\n" % l)
            fs = config[g]['files']
            fs = sorted(fs)
            for (pattern, local) in fs:
                self.output.write("  file '%s'\n" % pattern)
                if local:
                    self.output.write("   (read to '%s')\n" % local)

    def find_requests(self, requests):
        config = self.read_config()

        bundles = {}
        files = []

        def add(fs, bs):
            for f in fs:
                files.append(f)
            for (b, ms) in bs.items():
                bundles[b] = bundles.get(b, []) + ms

        if not requests:
            requests = ['']
        for request in list(requests):
            if request in config:
                add(config[request]['files'], config[request]['bundles'])
                requests.remove(request)
        for request in list(requests):
            for group_name in config.keys():
                if re.search(request, group_name):
                    self.output.write("No group named '%s', using '%s' instead.\n"
                                      % (request, group_name))
                    add(config[group_name]['files'], config[group_name]['bundles'])
                    try:
                        requests.remove(request)
                    except ValueError:
                        # Happens when request matches several groups.
                        pass
        for request in list(requests):
            for dict in config.values():
                for (b, ms) in dict['bundles'].items():
                    (pattern, local) = b
                    if re.search(request, pattern) or (local is not None and re.search(request, local)):
                        self.output.write("No group matching '%s',\n"
                                          "    using bundle '%s:%s' instead.\n"
                                          % (request, pattern, local))
                        add([], {(pattern, local): ms})
                        try:
                            requests.remove(request)
                        except ValueError:
                            # Happens when request matches several bundles.
                            pass
        for request in list(requests):
            for dict in config.values():
                for (pattern, local) in dict['files']:
                    if re.search(request, pattern) or (local is not None and re.search(request, local)):
                        self.output.write("No group or bundle matching '%s',\n"
                                          "    using file '%s:%s' instead.\n"
                                          % (request, pattern, local))
                        add([(pattern, local)], {})
                        try:
                            requests.remove(request)
                        except ValueError:
                            # Happens when request matches several files.
                            pass
        for request in list(requests):
            for dict in config.values():
                for (b, ms) in dict['bundles'].items():
                    for (pattern, local) in ms:
                        if re.search(request, pattern) or (local is not None and re.search(request, local)):
                            self.output.write("No group or bundle matching '%s',\n"
                                              "    using member '%s:%s'\n"
                                              "    of bundle '%s:%s' instead.\n"
                                              % (request, pattern, local, b[0], b[1]))
                            add([], {b: [(pattern, local)]})
                            try:
                                requests.remove(request)
                            except ValueError:
                                # Happens when request matches several members.
                                pass
        if requests:
            raise Error("Don't know how to fetch these items: %s" % requests)
        return bundles, files

    def fetch_one(self, url, local, members=None):
        if members is None:
            members = []
        m = re.match('([a-z]+)://([^/]+)/(.*/)([^/]+)$', url)
        if m is None:
            raise Error("Malformed URL '%s'" % url)
        protocol = m.group(1)
        if protocol in 'https http ftp'.split():
            self.fetch_url(url, local, members)
        elif protocol == 'ftpmatch':
            host = m.group(2)
            path = m.group(3)
            pattern = m.group(4)
            self.ftpmatch(host, path, pattern, local, members)
        else:
            raise Error("Unknown protocol '%s' in URL '%s'" % (protocol, url))

    def fetch_url(self, url, local, members=None):
        import os

        if local is None:
            local = url.split('/')[-1]
        name = os.path.join(self.prefix, local.strip())

        if os.path.exists(name) and os.path.getsize(name) == 0:
            self.output.write("%s is empty; removing it.\n" % name)
            os.remove(name)
        if os.path.exists(name) and not self.force:
            self.output.write("%s already exists.\n" % name)
        else:
            self.make_prefix()
            self.output.write("Fetching %s to %s\n" % (url, name))
            # We have to set a User-Agent header in order to
            # fetch GISTEMPv4_sources.tar.gz (the web server
            # rejects the HTTP request otherwise). urllib2 does
            # this, but urllib does not.
            remote = urllib.request.urlopen(url)
            # Check getcode(), but only for HTTP.
            if remote.getcode() and remote.getcode() != 200:
                raise Error(
                    "Fetching %s to %s failed (status code %s)." %
                    (url, name, remote.getcode()))
            with open(name, 'wb') as out:
                copy_progress(remote, out, self.output)
        if os.path.getsize(name) == 0:
            raise Error("%s is empty." % name)
        if members:
            self.extract(name, members)

    def ftpmatch(self, host, path, pattern, local, members):
        regexp = re.compile(pattern)
        # http://www.python.org/doc/2.4.4/lib/module-ftplib.html
        import ftplib

        remote = ftplib.FTP(host, 'ftp', 'info@climatecode.org')
        remote.cwd(path)
        dir = remote.nlst()
        good = filter(regexp.match, dir)
        good = sorted(good)

        if not good:
            raise Error("Could not find any file matching '%s' at ftp://%s/%s" % (pattern, host, path))
        remotename = good[-1]
        path = path.strip('/')
        self.fetch_url('ftp://%s/%s/%s' % (host, path, remotename), local, members)

    def extract(self, name, members):
        exts = name.split('.')
        if exts[-1] in 'gz bz bz2'.split():
            exts = exts[:-1]
        if exts[-1] in 'tar tgz tbz tbz2'.split():
            self.extract_tar(name, members)
        elif exts[-1] in 'zip'.split():
            self.extract_zip(name, members)
        elif name.endswith('.gz'):
            self.extract_gzip(name, members)
        else:
            raise Error("Can't extract members from this type of file: %r", name)

    def extract_tar(self, archive, members):
        """
        `archive` is the name of a tar file (possibly
        compressed). `members` is a list of members to extract
        from it.
        """

        # Could figure out compression type here, and pass it in to
        # tarfile.open, but apparently these days the tarfile module
        # does it for us.

        # Watch out for bugs in some version of Python.
        # See http://code.google.com/p/cccgistemp/issues/detail?id=26
        tar = tarfile.open(name=archive, mode='r')
        for info in tar:
            # would like to use 'any', but that requires Python 2.5
            matches = [member for member in members if re.search(member[0] + '$', info.name)]
            if matches:
                if len(matches) > 1:
                    raise Error("Multiple patterns match '%s': %s" % (info.name, matches))
                members.remove(matches[0])
                local = matches[0][1]
                if local is None:
                    local = info.name.split('/')[-1]
                local = os.path.join(self.prefix, local.strip())
                if os.path.exists(local) and not self.force:
                    self.output.write("  ... %s already exists.\n" % local)
                else:
                    self.make_prefix()
                    out = open(local, 'wb')
                    self.output.write("  ... %s from %s.\n" % (local, info.name))
                    # The following used to be simply
                    # ``out.writelines(tar.extractfile(info))``, but the Python2.4
                    # tarfile.py does not provide iteration support.
                    member = tar.extractfile(info)
                    while True:
                        buf = member.read(4096)
                        if not buf:
                            break
                        out.write(buf)
        if members:
            raise Error("Couldn't find these members in '%s': %s" % (archive, [member[0] for member in members]))

    def extract_zip(self, name, members):
        z = zipfile.ZipFile(name)
        for entry in z.namelist():
            matches = [member for member in members if re.search(member[0] + '$', entry)]
            if matches:
                if len(matches) > 1:
                    raise Error("Multiple patterns match '%s': %s" % (entry, matches))
                members.remove(matches[0])
                local = matches[0][1]
                if local is None:
                    local = entry.split('/')[-1]
                local = os.path.join(self.prefix, local.strip())
                if os.path.exists(local) and not self.force:
                    self.output.write("  ... %s already exists.\n" % local)
                else:
                    self.make_prefix()
                    # Only works for text files.
                    out = open(local, 'w')
                    self.output.write("  ... %s from %s.\n" % (local, entry))
                    src = z.open(entry)
                    while True:
                        s = src.read(4096)
                        if not s:
                            break
                        out.write(s)
                    out.close()
                    src.close()
        if members:
            raise Error("Couldn't find these members in '%s': %s" % (name, [member[0] for member in members]))

    def extract_gzip(self, name, members):
        import gzip
        import shutil

        if len(members) > 1:
            raise Error("Simple compressed file, %r, is only allowed exactly one member", name)

        if not members:
            basename = name.split('/')[-1]
            # remove final extension (should be '.gz')
            basename = '.'.join(basename.split('.')[:-1])
            members.append((basename, None))

        (local, _) = members[0]
        local = os.path.join(self.prefix, local)
        if os.path.exists(local) and os.path.getsize(local) == 0:
            self.output.write("%s is empty; removing it.\n" % local)
            os.remove(local)
        if os.path.exists(local) and not self.force:
            self.output.write("  ... %s already exists.\n" % local)
            return

        self.make_prefix()
        out = open(local, 'wb')
        inp = gzip.open(name)
        shutil.copyfileobj(inp, out)
        inp.close()
        out.close()


def copy_progress(source, destination, progress):
    """
    Copy the contents of open readable stream `source` to open
    writable stream `destination`. Progress is written to the
    open writable stream `progress`.

    Typically `source` will be a remote file fetched with
    urllib2.urlopen, and `destination` will be a local disk
    file. `progress` will be stderr.
    """

    try:
        content_length = int(source.info()['Content-Length'])
    except (AttributeError, KeyError, TypeError, ValueError):
        content_length = None

    got = 0
    while True:
        chunk = source.read(8000)
        got += len(chunk)
        if content_length is not None:
            outof = '/%d [%d%%]' % (
                content_length, 100 * got // content_length)
        else:
            outof = ''
        progress.write("\r  %d%s" % (got, outof))
        if not chunk:
            break
        destination.write(chunk)

    progress.write('\n')
    progress.flush()
    return 0


class Error(Exception):
    """Some sort of problem with fetch."""


# Guido's main, http://www.artima.com/weblogs/viewpost.jsp?thread=4829
class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


def main(argv=None):
    if argv is None:
        argv = sys.argv
    write_list = False
    kwargs = dict()
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "", ["help", "list", "force", "store=", "config="])
            for o, a in opts:
                if o in ('--help',):
                    print(__doc__)
                    return 0
                if o == '--list':
                    write_list = True
                if o == '--force':
                    kwargs.update(force=True)
                if o == '--config':
                    kwargs.update(config_file=a)
                if o == '--store':
                    kwargs.update(prefix=a)
        except getopt.error as msg:
            raise Usage(msg)
        kwargs.update(requests=args)
        fetcher = Fetcher(**kwargs)
        if write_list:
            fetcher.list_things()
        else:
            fetcher.fetch()
    except Usage as err:
        print(err.msg, file=sys.stderr)
        print("for help use --help", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
