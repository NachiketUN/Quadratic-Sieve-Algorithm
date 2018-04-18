
# utility function to find pow(base, exponent) % modulus
# def pow(base,exponent,mod):
# 	result=1
# 	base=base%mod
# 	while exponent > 0:

# 		if exponent%2==1:
# 			result=(result * base)%mod

# 		exponent=exponent>>1
# 		base=(base*base)%mod

# 	return result

def gcd(a,b):
	if b==0:
		return a
	else:
		return gcd(b, a%b)

# Returns k such that b^k = 1 (mod p)
def order(p,q):
	if gcd(p,q) !=1:
		print(p,' and ',q,'are not co-prime')
		return -1
	k=3
	while True:
		if pow(q,k,p)==1:
			return k
		k+=1

# function return p - 1 (= x argument) as x * 2^e,
# where x will be odd sending e as reference because
# updation is needed in actual e
def convertx2e(x):
	e=0
	while x%2==0:
		x/=2
		e+=1
	return int(x),e

# Main function for finding the modular square root
def STonelli(n,p):
	# a and p should be coprime for finding the modular
	# square root
	if gcd(n, p) != 1:
	
		print("a and p are not coprime\n")
		return -1


    # If below expression return (p - 1) then modular
	# square root is not possible
	#printf("%d\n",pow(n, (p - 1) / 2, p) );
	if pow(n, (p - 1) // 2, p) == (p - 1):

		print("no sqrt possible\n")
		return -1


	# expressing p - 1, in terms of s * 2^e, where s
	# is odd number
	s,e = convertx2e(p - 1)

	# finding smallest q such that q ^ ((p - 1) / 2)
	# (mod p) = p - 1
	q=2
	while True:

		# q - 1 is in place of (-1 % p)
		if pow(q, (p - 1) // 2, p) == (p - 1):
			break
		q+=1

	
	# Initializing variable x, b and g
	x = pow(n, (s + 1) // 2, p)
	b = pow(n, s, p)
	g = pow(q, s, p)

	r = e
 
	# keep looping until b become 1 or m becomes 0
	while True:

		for m in range(r):
		
			if order(p, b) == -1:
				return -1
			
			# finding m such that b^ (2^m) = 1
			if order(p, b) == 2**m:
				break
		if m == 0:
			return x

		# updating value of x, g and b according to
		# algorithm
		x = (x * pow(g, 2**(r - m - 1), p)) % p
		g = pow(g, 2**(r - m), p)
		b = (b * g) % p

		if b == 1:
			return x
		r = m

def main():
	n=2
	p=113
	print(STonelli(n,p))


if __name__ == '__main__':
	main()