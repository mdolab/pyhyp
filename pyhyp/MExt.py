import tempfile
import imp
import os
import shutil
import sys


def _tmp_pkg(dirName):
    """
    Create a temporary package.

    Returns (name, path)
    """
    while True:
        path = tempfile.mkdtemp(dir=dirName)
        name = os.path.basename(path)
        try:
            imp.find_module(name)
            # if name is found, delete and try again
            os.rmdir(path)
        except Exception:  # noqa: E722
            break
    init = open(os.path.join(path, "__init__.py"), "w")
    init.close()
    return name, path


class MExt(object):
    """
    Load a unique copy of a module that can be treated as a "class instance".
    """

    def __init__(self, name, path=None, debug=False):
        tmpdir = tempfile.gettempdir()
        self.name = name
        self.debug = debug
        # first find the "real" module on the "real" syspath
        srcfile, srcpath, srcdesc = imp.find_module(name, path)
        # now create a temp directory for the bogus package
        self._pkgname, self._pkgdir = _tmp_pkg(tmpdir)
        # copy the original module to the new package
        shutil.copy(srcpath, self._pkgdir)
        # add the directory containing the new package to the search path
        sys.path.append(tmpdir)
        # import the module
        # __import__ returns the package, not the sub-module
        self._pkg = __import__(self._pkgname, globals(), locals(), [self.name])
        # remove the bogus directory from sys.path
        sys.path.remove(tmpdir)
        # return the module object
        self._module = getattr(self._pkg, self.name)
        # now add the module's stuff to this class
        self.__dict__.update(self._module.__dict__)

    def __del__(self):
        # remove module if not in debug mode
        if not self.debug:
            del sys.modules[self._module.__name__]
            del sys.modules[self._pkg.__name__]

            # now try to delete the files and directory
            shutil.rmtree(self._pkgdir)
            # make sure the original module is loaded -
            # otherwise python crashes on exit
            # if MExt objects have not been explicitly 'del'd,
            # and __del__ is occurring on python shutdown, the import will fail
            # and the exception is caught here
            try:
                __import__(self.name)
            except ImportError:
                pass
