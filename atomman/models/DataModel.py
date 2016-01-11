#Standard library imports
import json
from collections import OrderedDict

#External library imports
import xmltodict

class DataModel():
    def __init__(self, model=None):
        
        if model is None:
            self.data = OrderedDict()
        else:
            try:
                self.loads(model)
            except:
                self.load(model)
    
    def __str__(self):
        return self.dumps()
    
    def load(self, file_name):
        if file_name[-5:] == '.json':
            with open(file_name,'r') as f:
                self.data = json.load(f, object_pairs_hook = OrderedDict, parse_int = long) 
                
        elif file_name[-4:] == '.xml':
            with open(file_name, 'r') as f:
                self.data = self.__to_numbers(xmltodict.parse(f.read()))
        
        else:
            raise ValueError('Unsupported file type for DataModel')
    
    def loads(self, value):
        try:
            self.data = json.loads(value, object_pairs_hook = OrderedDict, parse_int = long) 
        except:
            self.data = self.__to_numbers(xmltodict.parse(value))
    
    def find(self, key):
        return [val for val in self.__gen_dict_extract(key, self.data)]
        
    def dump(self, file_name):
        if file_name[-5:] == '.json':
            with open(file_name,'w') as f:
                f.write(json.dumps( self.data, indent = 4, separators = (',',': ') ))
                
        elif file_name[-4:] == '.xml':
            with open(file_name, 'w') as f:
                f.write(xmltodict.unparse(self.data, pretty=True))
        
        else:
            raise ValueError('Unsupported file type for DataModel')

    def dumps(self, format='json'):
        if format == 'json':
            return json.dumps( self.data, indent = 4, separators = (',',': ') )
        elif format == 'xml':
            return xmltodict.unparse(self.data, pretty=True)
        else:
            raise ValueError('Unsupported format type for DataModel')
            
    def ismodel(self, key):
        if key in self.data:
            return True
        else:
            return False
    
    def __to_numbers(self, var):
        if hasattr(var, 'iteritems'):
            for k, v in var.iteritems():
                var[k] = self.__to_numbers(v)
            return var
        elif hasattr(var, '__iter__'):
            for i in xrange(len(var)):
                var[i] = self.__to_numbers(var[i])
            return var
        elif isinstance(var, (str, unicode)):   
            if var == 'None' or var == 'null':
                return None
            elif var == 'True' or var == 'true':
                return True
            elif var == 'False' or var == 'false':
                return False
            else:
                try:
                    return long(var)
                except:
                    try:
                        return float(var)
                    except:
                        return var
        else:
            return var
            

    def __gen_dict_extract(self, key, var):
        if hasattr(var,'iteritems'):
            for k, v in var.iteritems():
                if k == key:
                    yield v
                if isinstance(v, dict):
                    for result in self.__gen_dict_extract(key, v):
                        yield result
                elif isinstance(v, list):
                    for d in v:
                        for result in self.__gen_dict_extract(key, d):
                            yield result