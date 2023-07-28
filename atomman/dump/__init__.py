# coding: utf-8
from importlib import import_module, resources

__all__ = ['dump', 'set_dump_styles']

# Set global dump_styles dict
dump_styles = {}
failed_dump_styles = {}

def dump(style, system, **kwargs):
    """
    Convert a System to another format.

    Parameters
    ----------
    style : str
        Indicates the format of the content to dump the atomman.System as.
    system : atomman.System
        The system to convert.
    kwargs : any, optional
        Any extra keyword arguments to pass to the underlying dump methods.

    Returns
    -------
    str, object or tuple
        Any content returned by the underlying dump methods.
    """

    if style in dump_styles:
        return dump_styles[style](system, **kwargs)
    elif style in failed_dump_styles:
        raise failed_dump_styles[style]
    else:
        raise ValueError(f'Unsupported dump style {style}')

def set_dump_styles():
    """
    Imports and sets the dump styles.  Should be called after importing the
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

            # Import module and set to dump_styles
            try:
                module = import_module(f'.{style}', __name__)
                dump_styles[style] = module.dump
            except ModuleNotFoundError as e:
                failed_dump_styles[style] = e
