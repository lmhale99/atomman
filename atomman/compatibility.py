from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
                        
import sys

# Python 2 settings
if sys.version_info[0] == 2:
    stringtype = basestring
    unicode = unicode
    long = long
    range = xrange
    
    def iteritems(d):
        for key, value in d.iteritems():
            yield key, value
    
# Python 3 settings
elif sys.version_info[0] == 3:
    stringtype = str
    unicode = str
    long = int
    range = range
    
    def iteritems(d):
        for key, value in d.items():
            yield key, value
    
else:
    raise ValueError("Unsupported Python version")