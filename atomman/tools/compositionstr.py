# coding: utf-8
# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

import numpy as np

from ..compatibility import unicode

def compositionstr(symbols, counts):
    """
    Generates a reduced composition string based on symbols and counts.
    
    Parameters
    ----------
    symbols : list
        The model symbol for a site.
    count : list
        How many sites are occupied by each symbol.
    
    Returns
    -------
    str
        The reduced composition string.
    """
    assert len(symbols) == len(counts), 'Symbols and counts must be the same length'
    
    sym_dict = {}
    for symbol, count in zip(symbols, counts):
        if symbol in sym_dict:
            sym_dict[symbol] += count
        else:
            sym_dict[symbol] = count
    
    gcd = np.gcd.reduce(list(sym_dict.values()))
    
    composition =''
    for symbol in sorted(sym_dict):
        count = sym_dict[symbol] // gcd
        if sym_dict[symbol] > 0:
            composition += symbol
            if count != 1:
                composition += unicode(count)
    
    return composition