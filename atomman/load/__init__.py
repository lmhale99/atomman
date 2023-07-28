# coding: utf-8
from importlib import import_module, resources

__all__ = ['load', 'FileFormatError']

# Set global load_styles dict
load_styles = {}
failed_load_styles = {}
class FileFormatError(Exception):
    pass

def load(style, *args, **kwargs):
    """
    Load a System from another format.
    
    Parameters
    ----------
    style : str
        Indicates the format of the content to load as an atomman.System
    args 
        Any positional-dependent arguments to pass to the underlying load methods.
    kwargs
        Any keyword arguments to pass to the underlying load methods.
        
    Returns
    -------
    system : atomman.System
        The system object associated with the data model.
    """
    
    if style in load_styles:
        return load_styles[style](*args, **kwargs)
    elif style in failed_load_styles:
        raise failed_load_styles[style]
    else:
        raise ValueError(f'Unsupported load style {style}')

def set_load_styles():
    """
    Imports and sets the load styles.  Should be called after importing the
    iprPy.load submodule.
    """
    # Define subfolders to ignore
    ignorelist = ['__pycache__']

    # Dynamically import calculation styles
    if hasattr(resources, 'files'):
        styles = (resource.name for resource in resources.files(__name__).iterdir())
    else:
        styles = resources.contents(__name__)
    for style in styles:

        # Only import subfolders not in ignorelist
        if '.' not in style and style not in ignorelist:
            
            # Import module and set to load_styles
            try:
                module = import_module(f'.{style}', __name__)
                load_styles[style] = module.load
            except ModuleNotFoundError as e:
                failed_load_styles[style] = e

set_load_styles()