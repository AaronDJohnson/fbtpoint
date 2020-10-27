import numpy as np
from mpmath import mpc, re, im, mp, mpf


def complex_cont_frac(a, b, tol=10**-(mp.dps)):
    """
    Computes a complex continued fraction using a sum of convergents.

    Parameters:
        a (complex array): input 1xn array (first element is 0)
        b (complex array): input 1xn array (first element is 0)

    Keyword Args:
        tol (float): defaults to mpmath precision

    Returns:
        res (complex)
    """
    n = len(a)
    A = np.zeros(n + 1) * mpc("1")
    B = np.zeros(n + 1) * mpc("1")
    x = np.zeros(n + 1) * mpc("1")
    A[0] = mpc("1")
    A[1] = b[0]
    B[0] = mpc("0")
    B[1] = mpc("1")
    x[1] = A[1] / B[1]
    # print(A[1])
    A[2] = b[1] * A[1] + a[1] * A[0]
    B[2] = b[1] * B[1] + a[1] * B[0]
    x[2] = A[2] / B[2]
    # print(A[2])
    A[3] = b[2] * A[2] + a[2] * A[1]
    B[3] = b[2] * B[2] + a[2] * B[1]
    x[3] = A[3] / B[3]
    # print(A[3])

    for i in range(3, n - 1):
        A[i + 1] = b[i] * A[i] + a[i] * A[i - 1]
        # print("A = ", A[i + 1])
        B[i + 1] = b[i] * B[i] + a[i] * B[i - 1]
        # print(B[i + 1])
        x[i + 1] = A[i + 1] / B[i + 1]
        # print(i)
        # print(x[i + 1])
        # print(x[i])
        re_err = abs(re(x[i + 1]) - re(x[i]))
        im_err = abs(im(x[i + 1]) - im(x[i]))
        # print(tol)
        if re_err < tol and im_err < tol:
            return x[i + 1]


def contfracK(a, b):
    """
    We know how long the vectors are in this case and we need the entire
    thing so this is simpler than the other one. Just calculate backwards!

    Parameters:
        a (float array): input 1xn array
        b (float array): input 1xn array

    Returns:
        res (float)
    """
    if len(a) != len(b):
        print("There is a mismatch between the len of vectors a and b")
        return 0.0
    result = np.zeros(len(a)) * mpf("1")
    result[0] = a[-1] / b[-1]
    end = len(b) - 1
    # print('a =', a[-1])
    # print('b =', b[-1])
    # print('res =', result[0])
    for i in range(1, len(result)):
        end = len(b) - 1
        # print('a =', a[end - i])
        # print('b =', b[end - i])
        # print('res =', result[i - 1])
        result[i] = (result[i - 1] + b[end - i])**(-1) * a[end - i]
        # print('res =', result[i])
    return result[-1]












