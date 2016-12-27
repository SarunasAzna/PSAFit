from math import log
import numpy as np


def logDerivative(x, data):
    """ Function for calculating logarithmic(10) derivative dy/d(log(x)).

    :param x Array : x values
    :param data 2D Array : y values
    :return 2D Numpy Array : values of logarithmic(10) derivative
    """
    ndata, nx = data.shape
    derivative = []
    for k in data:
        der = []
        for pos, y in enumerate(k):
            if pos < (len(k) - 1):
                d1 = (k[pos+1] - k[pos])/(log(x[pos+1], 10) - log(x[pos], 10))
            d2 = (k[pos] - k[pos-1])/(log(x[pos], 10) - log(x[pos-1], 10))
            if pos == 0:
                der.append(d1)
            elif pos == (len(k) - 1):
                der.append(d2)
            else:
                der.append(0.5*(d1 + d2))
        derivative.append(der)
    derivative = np.array(derivative)
    return derivative
