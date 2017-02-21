import pandas as pd
import warnings
import numpy as np

#atomman imports
from atomman.tools import uber_open_rmode
from DataModelDict import DataModelDict as DM

from .Log import Log


def log_extract(log_info):
    """Parses a LAMMPS screen output/log file and returns a data model containing the information."""
    
    warnings.simplefilter('always')
    warnings.warn('log_extract function is replaced with the Log class', DeprecationWarning)
    
    #Generate a Log object
    log = Log(log_info)
    return log.model
    
    #Handle file names, strings and open file-like objects equivalently
    #with uber_open_rmode(log_info) as log_info:
    
    #    headers = []
    #    footers = []
    #    i = 0
        
        #for all lines in file/output
    #    for line in log_info:
            
            #skip blank lines
    #        if len(line.split()) == 0:
    #            continue
                
            #This is listed before both run and minimize simulations    
    #        if 'Memory usage per processor =' in line:
    #            headers.append(i+1)
            
            #This follows both run and minimize simulations
    #        elif 'Loop time of' in line:
    #            footers.append(i-1)
            
    #        i += 1
        
        #Add last line to footers for incomplete logs
    #    footers.append(i)
        
    #    log_info.seek(0)
        
        #Create DataModelDict root
    #    log_dict = DM()
    #    log_dict['LAMMPS-log-thermo-data'] = DM()
        
        #for all lines in file/output
    #    for header, footer in zip(headers, footers):

            #Read thermo data
    #        df = pd.read_csv(log_info, header=header, nrows=footer-header, sep='\s+', engine='python', skip_blank_lines=True)
    #        log_info.seek(0)            

            #Convert to DataModelDict
    #        thermo = DM()
    #        for j in df:
    #            thermo[str(j)] = df[j].values.tolist()
            
            #Append simulation results to DataModelDict root
    #        simulation = DM([('thermo', thermo)])
    #        log_dict['LAMMPS-log-thermo-data'].append('simulation', simulation)
                
    #return log_dict     