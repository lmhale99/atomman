# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
import sys

# Boolean checks of Python version
ispython2 = sys.version_info[0] == 2
ispython3 = sys.version_info[0] == 3

# Python 2 settings
if ispython2:
    stringtype = basestring
    inttype = (int, long)
    unicode = unicode
    long = long
    int = long
    range = xrange
    
    def iteritems(d):
        for key, value in d.iteritems():
            yield key, value
    
# Python 3 settings
elif ispython3:
    stringtype = str
    unicode = str
    inttype = int
    int = int
    long = int
    range = range
    
    def iteritems(d):
        for key, value in d.items():
            yield key, value
    
else:
    raise ValueError("Unsupported Python version")