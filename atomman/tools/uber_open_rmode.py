from cStringIO import StringIO
import os

class uber_open_rmode():
    """
    context manager class for reading data from file-like objects, file names, 
    and data strings in the same manner.
    """
    
    def __init__(self, data):
        self.data = data
        
    def __enter__(self):

        if hasattr(self.data, 'read'):
            self.open_file = self.data
            self.to_close = False
        
        elif os.path.isfile(self.data):
            self.open_file = open(self.data)
            self.to_close = True
        
        else:
            self.open_file = StringIO(self.data)
            self.to_close = True
        
        return self.open_file
    
    def __exit__(self, *args):
        if self.to_close:
            self.open_file.close()