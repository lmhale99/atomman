import pandas as pd

import numpy as np

#atomman imports
from atomman.tools import uber_open_rmode
from DataModelDict import DataModelDict as DM

class Log(object):
    """Object for representing a LAMMPS log output"""
    
    def __init__(self, log_info=None):
        """Intitializes a Log object."""
        
        #Initialize simulations array
        self.__simulations = []
        
        #Read log data if supplied
        if log_info is not None:
            self.read(log_info)
        
    def read(self, log_info, append=True):
        """Parses a LAMMPS screen output/log file for thermodynamic data."""
        
        #Remove existing data if append is False
        if append is False:
            self.__simulations = []
        
        #Handle file names, strings and open file-like objects equivalently
        with uber_open_rmode(log_info) as log_info:
        
            headers = []
            footers = []
            i = 0
            
            #for all lines in file/output
            for line in log_info:
                
                #skip blank lines
                if len(line.split()) == 0:
                    continue
                    
                #This is listed before both run and minimize simulations    
                if 'Memory usage per processor =' in line:
                    headers.append(i+1)
                
                #This follows both run and minimize simulations
                elif 'Loop time of' in line:
                    footers.append(i-1)
                
                i += 1
            
            #Add last line to footers for incomplete logs
            footers.append(i)
            
            #Reset file pointer
            log_info.seek(0)
            
            #for all lines in file/output
            for header, footer in zip(headers, footers):
                
                #Initialize simulation data dictionary
                sim = {}
                    
                #Read thermo data and reset file pointer
                sim['thermo'] = pd.read_csv(log_info, header=header, nrows=footer-header, sep='\s+', engine='python', skip_blank_lines=True)
                log_info.seek(0)

                #Append simulation results
                self.__simulations.append(sim)        
                

    @property
    def simulations(self):
        return self.__simulations
    
    def flatten(self, style='last'):
        """Combines all simulations into one 
        
        Keyword Arguments:
        style -- specifies which values to use for duplicate time steps
            - first uses the values from the earlier simulation
            - last uses the values from the later simulation
        """
        
        for sim in self.simulations:
            if 'thermo' in sim and len(sim['thermo']) > 0:
                assert 'Step' in sim['thermo'], 'All simulation thermos must have Step key in order to flatten'
        merged_df = self.simulations[0]['thermo']
        for i in xrange(1, len(self.simulations)):
            if 'thermo' in self.simulations[i]:
                thermo = self.simulations[i]['thermo']
                if   style == 'first':
                    merged_df = pd.concat([merged_df, thermo[thermo.Step > merged_df.Step.max()]], ignore_index=True)
                elif style == 'last':
                    merged_df = pd.concat([merged_df[merged_df.Step < thermo.Step.min()], thermo], ignore_index=True)
                else:
                    raise ValueError('Unsupported style')
        self.__simulations = [{'thermo':merged_df}]
        
    def model(self, flatten=None):
        """Returns a DataModelDict containing the information."""
        if flatten is not None:
            self.flatten(flatten)
        
        #Create DataModelDict root
        log_model = DM()
        log_model['LAMMPS-log-thermo-data'] = DM()
        
        #Loop over all simulations
        for sim in self.simulations:
            sim_model = DM()
            
            #Convert to DataModelDict
            sim_model['thermo'] = DM()
            thermo = sim['thermo']
            for j in thermo:
                sim_model['thermo'][str(j)] = thermo[j].values.tolist()
                
            #Append simulation results to DataModelDict root
            log_model['LAMMPS-log-thermo-data'].append('simulation', sim_model)
            
        return log_model
    