# coding: utf-8
def screen_input(prompt=''):
    """
    Replacement input function that is compatible with the mingw terminal.
    
    Parameters
    ----------
    prompt : str, optional
        The screen prompt to use for asking for the input.
        
    Returns
    -------
    str
        The user input.
    """
    # Flush prompt to work interactively with mingw
    print(prompt, end=' ', flush=True)
    return input()