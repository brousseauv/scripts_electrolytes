from itertools import groupby
import numpy as np
''' Miscellanous useful functions '''

def sort_consecutive_groups(arr):
    ''' This functions sorts array (integer) elements into consecutive 
        groups and returns a list containing the first number from 
        each of these groups
    '''

    my_order = lambda x: x  # returns the value of a given element
    grouped = groupby(enumerate(list(arr)), key=lambda x: x[0] - my_order(x[1]))
    start = []
    end = []

    for key, group in grouped:
        g = np.asarray(list(group))
        # Prevent groups of single elements to be selected
        if np.shape(g)[0] > 1:
            start.append(g[0,1])
            end.append(g[-1,1])
    return start, end
