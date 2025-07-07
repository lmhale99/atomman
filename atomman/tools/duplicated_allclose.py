from typing import Union

import numpy as np
import pandas as pd

__all__ = ['duplicated_allclose', 'newrows_allclose']

def duplicated_allclose(dataframe: pd.DataFrame,
                        dcols: list,
                        fcols: dict,
                        keep: Union[str, bool] = 'first'):
    """
    Determine duplicate rows in a dataframe using both exact comparisons of
    values in the dcols columns, and close comparisons of values in the fcols
    columns within given absolute tolerances.  This allows for quick
    identification of duplicate entries when some of the values being checked
    are inexact floats.
    
    Implementation steps:
    1. Uses pandas.DataFrame.duplicated to ignore any rows with unique dcols
       combinations.
    2. Sorts the remaining rows on all dcols and fcols.
    3. Compares the values in each row of the sorted array to the values in
       the previous row to determine if they match.  For the fcols, this is
       done by computing the diff between each row and checking if the absolute
       value of the diff is less than the target tolerance.
    4. Identifies any rows where all values match the previous row.  Also, the
       list of any rows where all values match the next row is easily obtained
       by shifting the "previous match" list.  These lists are then used to 
       identify duplicates based on the keep argument value.
    
    Notes and cautions: 
    1. Since the algorithm only compares neighboring rows, false duplicates may
       be possible in the case where tolerance ranges overlap.  For instance,
       consider a value in rows 1, 2, 3.  The value in row 2 is considered a
       match for row 1 if it is within the tolerance, and similarly for row 3
       based on row 2.  However, it is then possible that the difference in
       values between rows 1 and 3 exceed the tolerance. Thus, row 3 is
       considered a duplicate even though it is distinct from row 1.
    2. Because of the sorting, it is not guaranteed that when you use
       keep='first' or keep='last' that the row being kept is always uniquely
       identified.  In the case of keep='first', the row not considered a
       duplicate will be the one with the slightly smaller values for the
       fcols.
       
    
    Parameters
    ----------
    dataframe : pandas.DataFrame
        The dataframe to search for duplicates
    dcols : list
        The column names that are to be tested for exact duplicates.
    fcols : dict
        The column names (keys) that are to be tested for close duplicates
        using absolute tolerances (values).
    keep : str or bool, optional
        Determines which duplicates (if any) to mark. 
        'first' : Mark duplicates as True except for the first occurrence.
        'last' : Mark duplicates as True except for the last occurrence.
        False : Mark all duplicates as True.

    Returns
    -------
    pandas.Series
        Boolean series indicating which rows are considered duplicates.
    """
    # Initialize isduplicated for all rows of the original dataframe as False
    isduplicated = pd.Series(np.full(len(dataframe), False, dtype=bool),
                             index=dataframe.index)
    
    # Filter out any rows that have unique dcols combinations
    if dcols:
        dcols_dupes = dataframe[dataframe.duplicated(subset=dcols, keep=False)]
    else:
        dcols_dupes = dataframe
    
    # Quick escape if all dcols are unique
    if len(dcols_dupes) == 0:
        return isduplicated

    # Sort the dupes by the dcols and fcols terms
    dcols_dupes = dcols_dupes.sort_values(by=dcols + list(fcols.keys()))
    
    # Compute the diff between fcols values for each neighboring row
    # fillna with inf so first row is inf instead of nan
    fdiff = dcols_dupes[fcols.keys()].diff().fillna(np.inf)
    
    # Identify which element diffs are within the tolerance for the fcols column name
    fcompare = fdiff.apply(func=lambda x: np.absolute(x) <= fcols[x.name], axis=0)
    
    def issame(column: pd.Series):
        """
        DataFrame apply method for identifying which elements in a column have
        the same value as the previous element in the column.  Note the first
        row is always set to False as there is no previous column to compare to.
        """
        column = column.values
        same = column[1:] == column[:-1]
        return np.append([False], same)
    
    # Identify which dcols elements match with the previous rows.
    dcompare = dcols_dupes[dcols].apply(issame)
    
    # Combine dcompare and fcompare and find which rows have all values matching
    # the previous row.
    allcompare = pd.concat([dcompare, fcompare], axis=1).all(axis=1)
    
    # dupes_index is index of all elements that may be dupes
    dupes_index = allcompare.index
    
    # Values of allcompare indicate which rows of dcols_dupes are same as previous
    sameasprevious = allcompare.values
    
    # Shift sameasprevious by one to get sameasnext
    sameasnext = np.append(sameasprevious[1:], [False])
    
    # Update isduplicated with the identified duplicate rows
    if keep == 'first':
        isduplicated.loc[dupes_index] = sameasprevious
        
    elif keep == 'last':
        isduplicated.loc[dupes_index] = sameasnext
        
    elif keep is False:
        isduplicated.loc[dupes_index] = np.any(np.vstack([sameasprevious, sameasnext]), axis=0)
        
    else:
        raise ValueError('keep must be "first", "last", or False')
    
    return isduplicated




def newrows_allclose(dataframe1: pd.DataFrame,
                     dataframe2: pd.DataFrame,
                     dcols: list,
                     fcols: dict):
    """
    Identifies which rows of dataframe2 are new and unique in comparison to
    other rows of dataframe2 and a reference dataframe1.
    """
    # Get count of dataframe1 rows
    len_dataframe1 = len(dataframe1)

    # Check rows in dataframe2 against itself, keeping first unique instances
    isduplicated = duplicated_allclose(dataframe2, dcols, fcols, keep='first')

    # Quick exit if dataframe1 is empty
    if len_dataframe1 == 0:
        return ~isduplicated

    # Pull out all unique rows of dataframe2
    tocheck = dataframe2[~isduplicated]

    # Combine dataframe1 with tocheck and check for duplicates
    dataframe12 = pd.concat([dataframe1, tocheck], ignore_index=True)
    isduplicated12 = duplicated_allclose(dataframe12, dcols, fcols, keep=False)

    # Extract the new isduplicated values associated with the tocheck rows
    existsin1 = isduplicated12.values[len_dataframe1:]

    # Update isduplicated with the tocheck rows
    isduplicated.loc[tocheck.index] = existsin1

    return ~isduplicated
