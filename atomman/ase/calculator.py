import importlib

# Manually specifies module paths for known external Calculators or unconventionally named ones
known_ase_calculators = {
    'AseQmmmManyqm': 'ase.calculators.ase_qmmm_manyqm',
}

def calculator(classname:str, *args, **kwargs):
    """
    Wrapper that allows for the initialization of an ase calculator in a
    standard way.

    Parameters
    ----------
    classname : str
        The name of the calculator class to load. For calculators defined in
        ase, giving only the class name (e.g. 'EAM', 'Vasp') should work. Otherwise,
        giving the full package.module.class path will always work.
    *args, **kwargs : any
        Any of the positional or key-word arguments associated with initializing
        an object of the specified Calculator class.
    """
    
    if '.' not in classname:
        
        if classname in known_ase_calculators:
            # Get modulename if known
            modulename = known_ase_calculators[classname]
        else:
            # Guess modulename if not known
            modulename = f'ase.calculators.{classname.lower()}'
    else:
        # Split full class name
        terms = classname.split('.')
        modulename = '.'.join(terms[:-1])
        classname = terms[-1]

    # Load module and get class
    module = importlib.import_module(modulename)
    cls = getattr(module, classname)

    # Init a new object of the class using args, kwargs
    return cls(*args, **kwargs)