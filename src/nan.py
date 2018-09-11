def array_isnan(array):
    import numpy as np
    for i in range(len(array)):
        if (np.isnan(array[i])):
            return True
    return False