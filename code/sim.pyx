import numpy as np
cimport numpy as np
cimport cython

ctypedef np.float32_t DTYPE_t

@cython.boundscheck(False) 
@cython.wraparound(False) 
cpdef ssim(np.ndarray[DTYPE_t, ndim=2] im1, \
            np.ndarray[DTYPE_t, ndim=2] im2, \
            float data_range=1, \
            bint full=False, \
            float K1=0.01):

    cdef float C1
    cdef np.ndarray[DTYPE_t, ndim=2] S = np.empty_like(im1)
    cdef float mssim 

    C1 = (K1 * data_range) ** 2

    cdef Py_ssize_t dim1 = S.shape[0]
    cdef Py_ssize_t dim2 = S.shape[1]
    cdef Py_ssize_t i,j
    for i in range(dim1):
        for j in range(dim2):
            S[i,j] = ((2 * (im1[i,j] * im2[i,j])) + C1) / (im1[i,j] ** 2 + im2[i,j] ** 2 + C1)

    mssim = S.mean()

    if full:
        return mssim, S
    else:
        return mssim

################################################################