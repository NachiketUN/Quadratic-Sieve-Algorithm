from random import randrange

""" This module is used to test the primality of a number using the
    Rabin-Miller primality test
"""

def is_probable_prime(a):
    """ Perform Rabin-Miller primality test to determine whether given number
        is prime. Return True if number is very likely to be a prime, and False
        if it is definitely composite

    Arguments:
        a - Number to be tested
    """
    if a == 2:
        return True

    if a == 1 or a % 2 == 0:
        return False

    return rabin_miller_primality_test(a, 50)

def rabin_miller_primality_test(a, iterations):
    """ Rabin Miller primality test

    Arguments:
        a - Number to be tested
        iterations - Number of iterations allowed
    """
    r, s = 0, a - 1

    while s % 2 == 0:
        r += 1
        s //= 2

    for _ in range(iterations):
        n = randrange(2, a - 1)
        x = pow(n, s, a)
        if x == 1 or x == a - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, a)
            if x == a - 1:
                break
        else:
            return False
    return True