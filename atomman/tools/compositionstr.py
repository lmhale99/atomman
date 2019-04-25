# coding: utf-8
# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

from ..compatibility import range, int, unicode

def compositionstr(symbols, counts):
    """
    Generates a composition string based on symbols and their counts.
    
    Parameters
    ----------
    symbols : list
        All unique element model symbols.
    count : list
        How many unique sites are occupied by each symbol.
    
    Returns
    -------
    str
        The composition string.
    """
    primes = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47]
    
    sym_dict = {}
    for i in range(len(symbols)):
        sym_dict[symbols[i]] = counts[i]
    
    for prime in primes:
        if max(sym_dict.values()) < prime:
            break
        
        while True:
            breaktime = False
            for value in sym_dict.values():
                if value % prime != 0:
                    breaktime = True
                    break
            if breaktime:
                break
            for key in sym_dict:
                sym_dict[key] /= prime
    
    composition =''
    for key in sorted(sym_dict):
        if sym_dict[key] > 0:
            composition += key
            if sym_dict[key] != 1:
                composition += unicode(int(sym_dict[key]))
    
    return composition