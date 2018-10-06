"""This module contains helper functions for all the other modules used in this
    project
"""

def product(lst):
    """ Return product of all numbers in a list
    """
    prod = 1
    for _ in lst:
        prod *= _
    return prod

def gcd(a, b):
    """ Returns the greatest common divisor of the two numbers
    """
    while b:
        a, b = b, a % b
    return a

def kth_iroot(n, k):
    """ Return the integer k-th root of a number by Newton's method
    """
    u = n
    s = n + 1
    while u < s:
        s = u
        t = (k - 1) * s + n // pow(s, k - 1)
        u = t // k
    return s

def isqrt(n):
    """ Return the square root of a number rounded down to the closest integer
    """
    if n < 0:
        raise ValueError("Square root of negative number!")
    x = int(n)
    if n == 0:
        return 0
    a, b = divmod(x.bit_length(), 2)
    n = 2 ** (a + b)
    while True:
        y = (n + x // n) // 2
        if y >= x:
            return x
        x = y