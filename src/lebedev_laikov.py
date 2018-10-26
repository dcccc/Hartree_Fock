# To use this script, first compile the code using:
#
#     cc -fPIC -shared -o liblebedevlaikov.so src/Lebedev-Laikov.c
#
# http://ccl.net/cca/software/SOURCES/FORTRAN/Lebedev-Laikov-Grids/Lebedev-Laikov.F
#
import ctypes
import numpy as np

LIB = ctypes.CDLL("./liblebedevlaikov.so")

LDNS = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434,
        590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890,
        4334, 4802, 5294, 5810]

def ld(n):
    '''Returns (x, y, z) coordinates along with weight w.'''
    if n not in LDNS:
        raise ValueError("n = {} not supported".format(n))
    xyzw = np.zeros((4, n), dtype=np.float64)
    c_double_p = ctypes.POINTER(ctypes.c_double)
    n_out = ctypes.c_int(0)
    getattr(LIB, "ld{:04}_".format(n))(
        xyzw[0].ctypes.data_as(c_double_p),
        xyzw[1].ctypes.data_as(c_double_p),
        xyzw[2].ctypes.data_as(c_double_p),
        xyzw[3].ctypes.data_as(c_double_p),
        ctypes.byref(n_out))
    assert n == n_out.value
    return(xyzw)

# def test():
#     for n in LDNS:
#         xyzw = ld(n)
#         # print(abs(w.sum() - 1.0) < 1e-13)
#         # print(all(abs(np.linalg.norm([x, y, z], axis=0) - 1.0) < 1e-15))
#         # print(x,y,z,w)
#         np.savetxt("111.csv",xyzw.T,delimiter=',')

# test()
