import helpers
from random import randint

def is_probable_prime(a):
    """ Perform Rabin-Miller primality test to determine whether given number
        is prime. Return True if number is very likely to be a prime, and False
        if it is definitely composite
    """
    if a == 2:
        return True

    if a == 1 or a % 2 == 0:
        return False

    return rabin_miller_primality_test(a, 50)

def rabin_miller_primality_test(a, iterations):
    """ Rabin Miller primality test
    """
    r, s = 0, a - 1

    while s % 2 == 0:
        r += 1
        s //= 2

    for _ in range(iterations):
        n = randint(2, a - 1)
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

def check_perfect_power(n):
    """ Check if the given integer is a perfect power. If yes, return (r, b)
        such that r^b = n, otherwise return None.
        Assume that global small_primes has already been initialised and that n
        does not have any prime factors from small_primes.
    """
    prime = small_primes[-1]
    for p in small_primes:
        pth_root = helpers.kth_iroot(n, p)
        if pth_root < prime:
            break
        if pth_root ** p == n:
            return (pth_root, p)
    return None

def check_factor(n, i, factors):
    """ Checks whether 'i' is a factor of 'n' and adds 'i' to 'factors' if true
        by trial division
    """
    while n % i == 0:
        n //= i
        factors.append(i)
        if is_probable_prime(n):
            factors.append(n)
            n = 1
    return n

def find_small_primes(n, upper_bound):
    """ Perform trial division on the given number using all the primes up
        to the upper bound. Initialize the global variable 'small_primes' with
        a list of all the primes <= upper_bound.
    """
    print("Trial division and initializing small primes...")
    global small_primes
    is_prime = [True] * (upper_bound + 1)
    is_prime[0:2] = [False] * 2
    factors = []
    small_primes = []
    max_i = helpers.isqrt(upper_bound)
    rem = n
    for i in range(2, max_i + 1):
        if is_prime[i]:
            small_primes.append(i)
            rem = check_factor(rem, i, factors)
            if rem == 1:
                return factors, 1

            for j in range(i ** 2, upper_bound + 1, i):
                is_prime[j] = False

    for i in range(max_i + 1, upper_bound + 1):
        if is_prime[i]:
            small_primes.append(i)
            rem = check_factor(rem, i, factors)
            if rem == 1:
                return factors, 1

    print("Primes initialised!")
    return factors, rem

def find_prime_factors(n):
    """ Return one or more prime factors of the given number n.
        Assume that n is not a prime and does not have very small factors, and
        that the global small_primes has already been initialised. Do not
        return duplicate factors.
    """
    print("Checking whether {} is a perfect power ...".format(n))
    perfect_power = check_perfect_power(n)
    if perfect_power:
        print("{} is {}^{}".format(n, perfect_power[0], perfect_power[1]))
        factors = perfect_power[0]
    else:
        print("Not a perfect power")
        digits = len(str(n))
        if digits <= 30:
            print("Using Brent's variant of Pollard's rho factorization " + \
                  "algorithm to factorise {} ({} digits)".format(n, digits))
            factors = [brent_factorise(n)]
        else:
            print("Using Self-Initialising Quadratic Sieve to factorise " + \
                  "{} ({} digits)".format(n, digits))
            #factors = siqs_factorise(n)

    prime_factors = []
    for f in set(factors):
        for pf in find_all_prime_factors(f):
            prime_factors.append(pf)

    return prime_factors

def find_all_prime_factors(n):
    """ Return all prime factors of the given number n.
        Assume that n does not have very small factors and that the global
        small_primes has already been initialised.
    """
    rem = n
    factors = []

    while rem > 1:
        if is_probable_prime(rem):
            factors.append(rem)
            break

        for f in find_prime_factors(rem):
            print("Prime factor found: {}".format(f))
            assert is_probable_prime(f)
            assert rem % f == 0
            while rem % f == 0:
                rem //= f
                factors.append(f)

    return factors

def _pollard_brent_func(c, n, x):
    """ Return f(x) = (x^2 + c) % n
        Assume c < n
    """
    y = (x ** 2) % n + c
    if y >= n:
        y -= n

    assert y >= 0 and y < n
    return y

def brent_factorise(n, iterations=None):
    """ Perform Brent's variant of Pollard's rho factorization algorithm to
        attempt to find a non-trivial factor of the given number number, n.
        If iterations > 0, return None if no factors are found within its range
    """
    y, c, m = (randint(1, n - 1) for _ in range(3))
    r, q, g = 1, 1, 1
    i = 0
    while g == 1:
        x = y
        for _ in range(r):
            y = _pollard_brent_func(c, n, y)
        k = 0
        while k < r and g == 1:
            ys = y
            for _ in range(min(m, r - k)):
                y = _pollard_brent_func(c, n, y)
                q = (q * abs(x - y)) %  n
            g = helpers.gcd(q, n)
            k += m
        r *= 2
        if iterations:
            i += 1
            if i == iterations:
                return None

    if g == n:
        while True:
            ys = _pollard_brent_func(c, n, ys)
            g = helpers.gcd(abs(x - ys), n)
            if g > 1:
                break
    return g

def pollard_brent_iterator(n, factors):
    """ Iterator function for Brent's variant of Pollard's rho factorization
        algorithm to find all small prime factors. Restart every time a factor
        is found.
        Return 1 if all prime factors are found, or otherwise the remaining
        factor
    """
    rem = n
    while True:
        if is_probable_prime(n):
            factors.append(n)
            rem = 1
            break

        digits = len(str(n))
        if digits < 45:
            iterations = 20
        else:
            iterations = 25

            f = brent_factorise(rem, iterations)
            if f and f < rem:
                if is_probable_prime(f):
                    print("Brent's (Pollard's rho): Prime factor found: " + \
                        "{}".format(f))
                    factors.append(f)
                    rem //= f
                else:
                    print("Brent's (Pollard's rho): Composite factor " + \
                          "found: {}".format(f))
                    rem_f = pollard_brent_iterator(f, factors)
                    rem = (rem // f) * rem_f
            else:
                print("No more small factors found")
                break
    return rem

def factorise(n):
    if type(n) != int or n < 1:
        raise ValueError("Number must be a POSITIVE INTEGER")

    print("Factorizing {} ({} digits)...".format(n, len(str(n))))

    if n == 1:
        return []

    if is_probable_prime(n):
        return [n]

    factors, rem = find_small_primes(n, 1000000)

    if factors:
        print("Prime factors found so far:")
        factors_temp = []
        for _ in factors:
            if _ not in factors_temp:
                factors_temp.append(_)
        print(*factors_temp, sep=', ')
    else:
        print("No small factors found!")

    if rem != 1:
        digits = len(str(rem))
        if digits > 30:
            print("Attempting Quick Pollard's rho (Brent's variation) to " + \
                  "find slightly larger factors...")
            rem = pollard_brent_iterator(rem, factors)
        if rem > 1:
            for f in find_all_prime_factors(rem):
                factors.append(f)

    factors.sort()
    assert helpers.product(factors) == n
    for p in factors:
        assert is_probable_prime(p)
    return factors

def main():
    n = int(input("Enter the number to be factorized: "))
    result = factorise(n)
    new_result = []

    for _ in result:
        if _ not in new_result:
            new_result.append(_)

    print("\nPrime factors: {}".format(new_result))

if __name__ == '__main__':
    main()