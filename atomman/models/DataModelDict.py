from collections import OrderedDict
import json
import xmltodict
from copy import deepcopy

class DataModelDict(OrderedDict):
    
    def __init__(self, *args, **kwargs):
        OrderedDict.__init__(self)
        if len(args) == 1:
            try:
                if 'parse_float' not in kwargs:
                    parse_float = None
                else:
                    parse_float = kwargs['parse_float']
                if 'parse_int' not in kwargs:
                    parse_int = None
                else:
                    parse_int = kwargs['parse_int']
                self.load(args[0], parse_float=parse_float, parse_int=parse_int)
            except:
                self.update(*args, **kwargs)
        else:
            self.update(*args, **kwargs)
            
    def load(self, in_model, parse_float=None, parse_int=None):
        #load from string
        if isinstance(in_model, (str, unicode)):
            
            #try json interpreter
            try:
                self.update(json.loads(in_model, object_pairs_hook=DataModelDict, parse_int=parse_int, parse_float=parse_float))
            except:
                
                #try xml interpreter
                if True:
                    if parse_float is None:
                        parse_float = float
                    if parse_int is None:
                        parse_int = int
                    hold = xmltodict.parse(in_model)
                    self.update(self.__xml_val_parse(hold, parse_int, parse_float))
                else:
                    raise ValueError('String is not acceptable json or xml')
        
        #load from file
        elif hasattr(in_model, 'read'):
            #try json interpreter
            try:
                self.update(json.load(in_model, object_pairs_hook=DataModelDict, parse_int=parse_int, parse_float=parse_float))
            except:
                in_model.seek(0)
                #try xml interpreter
                if True:
                    if parse_float is None:
                        parse_float = float
                    if parse_int is None:
                        parse_int = int
                    hold = xmltodict.parse(in_model)
                    self.update(self.__xml_val_parse(hold, parse_int, parse_float))
                else:
                    raise ValueError('File is not acceptable json or xml')
        else:
            raise TypeError('Can only load json/xml strings or file-like objects')
            
    def find(self, key):
        return [val for val in self.__gen_dict_extract(key, self)]       
    
    def json(self, fp=None, indent=None, separators=None):
        if fp is None:
            return json.dumps( self, indent=indent, separators=separators)
        else:
            json.dump(self, fp, indent=indent, separators=separators)
    
    def xml(self, fp=None, indent=None, full_document=True):
        if indent is None:
            indent = ''
            newl = ''
        else:
            indent = ''.join([' ' for i in xrange(indent)])
            newl = '\n'

        return xmltodict.unparse(self, output=fp, pretty=True, indent=indent, newl=newl, full_document=full_document)

    def ismodel(self, key):
        if key in self:
            return True
        else:
            return False
    
    def __xml_val_parse(self, var, parse_int, parse_float):
        if hasattr(var, 'iteritems'):
            for k, v in var.iteritems():
                var[k] = self.__xml_val_parse(v, parse_int, parse_float)
            return DataModelDict(var)
        elif hasattr(var, '__iter__'):
            for i in xrange(len(var)):
                var[i] = self.__xml_val_parse(var[i], parse_int, parse_float)
            return var
        elif isinstance(var, (str, unicode)):   
            if var == '':
                return None
            elif var == 'True' or var == 'true':
                return True
            elif var == 'False' or var == 'false':
                return False
            elif var == '-Infinity' or var == '-Inf' or var == '-inf':
                return float('-Inf')            
            elif var == 'Infinity' or var == 'Inf' or var == 'inf':
                return float('Inf')            
            elif var == 'NaN' or var == 'nan':
                return float('NaN')
            else:
                try:
                    return parse_int(var)
                except:
                    try:
                        return parse_float(var)
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
    
    def key_to_html(self, key, var):
        var = deepcopy(var)
        if hasattr(var,'iteritems'):
            for k, v in var.iteritems():
                if k == key:
                    if isinstance(v, dict):
                        var[k] = v.xml(full_document=False)
                    elif isinstance(v, list):
                        for d in xrange(len(v)):
                            var[k][d] = v[d].xml(full_document=False)
                        var[k] = ''.join(var[k])
                else:
                    if isinstance(v, dict):
                        var[k] = self.key_to_html(key, v)
                    elif isinstance(v, list):
                        for d in xrange(len(v)):
                            var[k][d] = self.key_to_html(key, v[d])
        return var    
        
    
    
    