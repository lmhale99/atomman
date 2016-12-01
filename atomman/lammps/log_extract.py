import pandas as pd

#External library imports
from StringIO import StringIO

import numpy as np

#atomman imports
from DataModelDict import DataModelDict as DM

def log_extract(log_info):
    """Parses a LAMMPS screen output/log file and returns a data model containing the information."""
    
    if isinstance(log_info, (str, unicode)):
        log_info = StringIO(log_info)
    
    headers = []
    footers = []
    i = 0
    #for all lines in file/output
    for line in log_info:
        terms = line.split()
        
        if 'Memory usage per processor =' in line:
            headers.append(i+1)
        elif 'Loop time of' in line:
            footers.append(i)
        i += 1
    log_info.seek(0)
    
    #Create DataModelDict root
    log_dict = DM()
    log_dict['LAMMPS-log-thermo-data'] = DM()
    
    #for all lines in file/output
    for header, footer in zip(headers, footers):
        #print header, footer
        df = pd.read_csv(log_info, header=header,skipfooter=i-footer, sep='\s+', engine='python', skip_blank_lines=False)
        log_info.seek(0)            
        #print df
        thermo = DM()
        for j in df:
            thermo[str(j)] = df[j].values.tolist()
        
        simulation = DM([('thermo', thermo)])
        log_dict['LAMMPS-log-thermo-data'].append('simulation', simulation)
                
    return log_dict     