from . import has_ase

if has_ase:
    import ase
    import ase.optimize
else:
    from . import dummy as ase
    
def optimizer(name: str,
              atoms: ase.Atoms,
              **kwargs):
    """
    Wrapper method for the ase.optimize module allowing for optimizers to be specified by name.
    """
    if not has_ase:
        raise ModuleNotFoundError('this module requires ase to be installed')

    if name.lower().startswith('scipy'):
        # Set root and find optimizer names for scipy methods
        rootmodule = ase.optimize.sciopt
        optimizer_names = []
        for attr in dir(ase.optimize.sciopt):
            if attr.startswith('SciPy') and 'Optimizer' not in attr:
                optimizer_names.append(attr)
    else:
        # Set root and find optimizer names for default ase methods
        rootmodule = ase.optimize
        optimizer_names = ase.optimize.__all__

    # Match string name to optimizer class
    match = False
    for optimizer_name in optimizer_names:
        if name.lower() == optimizer_name.lower():
            match = True
            break
    if not match:
        raise ValueError(f'unrecognized optimizer style {name}')

    # Get the optimizer class
    optclass = getattr(rootmodule, optimizer_name)

    # Build the optimizer with the given kwargs
    return optclass(atoms, **kwargs)