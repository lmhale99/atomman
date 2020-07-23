# coding: utf-8

# atomman imports
from .. import load_system_model
from ... import Library
from ...tools import screen_input

import pandas as pd

def load(potential, family=None, method='dynamic', standing='good'):
    library = Library()
    
    mquery = {}

    # Add potential to query
    mquery['relaxed-crystal.potential-LAMMPS.key'] = potential.key

    # Add family to query
    if family is not None:
        mquery['relaxed-crystal.system-info.family'] = family

    # Add relaxation method to query
    if method is not None:
        mquery['relaxed-crystal.method'] = method

    # Add standing to query
    if standing is not None:
        mquery['relaxed-crystal.standing'] = standing
    
    records = library.potdb.get_records(template='relaxed_crystal', mongoquery=mquery)
    
    if len(records) == 1:
        record = records[0]
        
    elif len(records) > 1 and len(records) < 100:
        print(f'{len(records)} matching crystals found.')
        
        # Build DataFrame
        df = []
        for record in records:
            series = {} 
            series['family'] = record['relaxed-crystal']['system-info']['family']
            series['a'] = float(record['relaxed-crystal']['atomic-system']['box']['avect']['value'][0])
            series['ecoh'] = float(record['relaxed-crystal']['cohesive-energy']['value'])
            series['symbols'] = ''.join(record['relaxed-crystal']['atomic-system'].aslist('atom-type-symbol'))
            series['method'] = record['relaxed-crystal']['method']
            series['standing'] = record['relaxed-crystal']['standing']
            df.append(series)
        df = pd.DataFrame(df)
        
        # Sort
        df = df.sort_values(['ecoh'])
        
        # Create list header
        header = '#  family               symbols  alat    Ecoh    '
        if method is None:
            header += 'method  '
        if standing is None:
            header += 'standing'
        print(header)
        
        # Generate options list
        for i in range(len(df)):
            series = df.iloc[i]     
            
            row = f'{i+1:2} {series.family:20.20} {series.symbols:8.8} {series.a:7.4f} {series.ecoh:7.4f} '
            
            if method is None:
                row += f'{series.method:7.7} '
            if standing is None:
                row += f'{series.standing:4.4}'
            print(row)

        choice = int(screen_input('Select which #:'))
        record = records[df.iloc[choice-1].name]
        
    elif len(records) >= 100:
        raise ValueError('More than 100 matching crystals found')
    else:
        raise ValueError('No matching crystals found')
    
    ucell = load_system_model(record)
        
    return ucell