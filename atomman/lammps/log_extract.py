
#External library imports
import numpy as np

#atomman imports
from atomman.models import DataModelDict

def log_extract(log_info):
    """Parses a LAMMPS screen output/log file and returns a data model containing the information."""
    
    #Create DataModelDict root
    log_dict = DataModelDict()
    log_dict['LAMMPS-log-thermo-data'] = DataModelDict()
    log_dict['LAMMPS-log-thermo-data']['simulation'] = []
    
    #Initialize necessary parameters
    thermoread = False
    thermolist = None
    
    #Convert string to list if needed
    if isinstance(log_info, (str, unicode)):
        log_info = log_info.split('\n')    
    
    #for all lines in file/output
    for line in log_info:
        terms = line.split()
     
        #If the line has terms
        if len(terms)>0:
 
            #thermoread indicates time to read data
            if thermoread:
                
                #save values in list if they are numbers
                if terms[0].isdigit():
                    if thermolist is None:
                        thermolist = np.array([terms], dtype='float64')
                    else:
                        thermolist = np.append(thermolist, np.array([terms], dtype='float64'), axis=0)                    
                
                #if not values, transform to dictionary
                else:
                    thermo = DataModelDict()
                    for i in xrange(len(headers)):
                        thermoval = thermolist[:,i]
                        #if np.allclose(thermoval, np.asarray(thermoval, dtype='int32')):
                        #    thermoval = np.asarray(thermoval, dtype='int32')
                        thermo[headers[i]] = list(thermoval)
                        
                    #add dictionary to simulation list
                    log_dict['LAMMPS-log-thermo-data']['simulation'].append(DataModelDict([('thermo',thermo)]))                    
                    
                    #reset thermoread and thermolist
                    thermolist = None
                    thermoread = False    
                    
                    
            #If the line starts with Step, then save headers and set thermoread
            elif terms[0] == 'Step':
                headers = terms
                thermoread = True
                   
    #downsize simulation if it only has one term
    if len(log_dict['LAMMPS-log-thermo-data']['simulation']) == 1:
        log_dict['LAMMPS-log-thermo-data']['simulation'] = log_dict['LAMMPS-log-thermo-data']['simulation'][0]
    
    return log_dict     