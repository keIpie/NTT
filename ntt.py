#!/usr/bin/env sage

from sage.all import *
import itertools

def test_ntt():
    print("Implementation of Number Theoretic Transform")
    #test_ntt_forward()
    #test_ntt_backward()
    #test_ntt_radix_2()
    test_ntt_radix_3()
    #test_transform_radix2_vs_naive()
    #test_ntt_multiply()
    #test_ntt_multiply_wrapped()
    #test_ntt_multiply_3()
    #test_simple_convolution()

def ntt(input, q):
    print("Performing Number Theoretic Transform")
    w = 0
    if Zmod(q).multiplicative_group_is_cyclic():
        w = Zmod(q).multiplicative_generator()
    else:
        print("The group is not cyclic.")
        return []
    output = []
    for i in range(len(input)):
        tmp = 0
        for (j, val) in enumerate(input):
            tmp = (tmp + val*w**(i*j)) % q
        output.append(tmp)
    return output

def ntt_w(input, w, q):
    print("Performing Number Theoretic Transform")
    output = []
    for i in range(len(input)):
        tmp = 0
        for (j, val) in enumerate(input):
            tmp = (tmp + val*w**(i*j)) % q
        output.append(tmp)
    return output

def inverse_ntt_w(input, w, q):
    print("Performing Inverse Number Theoretic Transform")
    output = ntt_w(input, inverse_mod(int(w), q), q)
    scaler = inverse_mod(len(input), q)
    return [(val*scaler % q) for val in output]

def reverse(x, bits):
	y = 0
	for i in range(bits):
		y = (y << 1) | (x & 1)
		x >>= 1
	return y

def ntt_radix_2(vector, root, mod):
    n = len(vector)
    levels = n.bit_length() - 1
    if 1 << levels != n:
        raise ValueError("Length is not a power of 2")

    print(vector)
    for i in range(n):
        j = reverse(i, levels)
        if j > i:
            vector[i], vector[j] = vector[j], vector[i]
    print(vector)

    powtable = []
    temp = 1
    for i in range(n // 2):
        powtable.append(temp)
        temp = (temp * root) % mod

    size = 2

    while size <= n:
        halfsize = size // 2
        tablestep = n // size
        for i in range(0, n, size):
            k = 0
            for j in range(i, i + halfsize):
                l = j + halfsize
                left = vector[j]
                right = vector[l] * powtable[k]
                vector[j] = (left + right) % mod
                vector[l] = (left - right) % mod
                k += tablestep
        size *= 2

def ntt_radix_3(vector, root, mod):
    print("To be done... not working yet")

    n = 9
    l = n // 3

    powtable = []
    temp = 1
    for i in range(l):
        powtable.append(temp)
        temp = (temp * root) % mod

    p = 2
    m = 1
    for t in range(0, p+1, 1):
        for j in range(0,l,1):
            for k in range(0,m,1):
                c0 = vector[k+j*m]
                c1 = vector[k+j*m+l*m]
                c2 = vector[k+j*m+2*l*m]
                d0 = c1 + c2
                d1 = c0 - int(GF(mod)(2)**(-1)) * d0
                d2 = - int(GF(mod)(2)**(-1)) * (c1 - c2)
                vector[k+3*j*m] = (c0 + d0) % mod
                vector[k+3*j*m+m] = ((d1 + d2) * powtable[j]) % mod
                vector[k+3*j*m+2*m] = ((d1 - d2) * powtable[j]**2) % mod
        l = l//3
        m = m*3

def ntt_multiply(input1, input2, w, q):
    if not (0 < len(input1) == len(input2)):
        raise ValueError()
    if any((val < 0) for val in itertools.chain(input1, input2)):
        raise ValueError()
    tmp1 = ntt_w(input1, w, q)
    tmp2 = ntt_w(input2, w, q)
    tmp_out = [(x * y % q) for (x, y) in zip(tmp1, tmp2)]
    return inverse_ntt_w(tmp_out, w, q)

def ntt_multiply_wrapped(input1, input2, w, q):
    if not (0 < len(input1) == len(input2)):
        raise ValueError()
    if any((val < 0) for val in itertools.chain(input1, input2)):
        raise ValueError()
    input1_wrapped = []
    input2_wrapped = []
    for i in range(0,len(input1),1):
        input1_wrapped.append(((-w)**i * input1[i]) %q)
        input2_wrapped.append(((-w)**i * input2[i]) %q)
    tmp1 = ntt_w(input1_wrapped, w, q)
    tmp2 = ntt_w(input2_wrapped, w, q)
    tmp_out = [(x * y % q) for (x, y) in zip(tmp1, tmp2)]
    result = inverse_ntt_w(tmp_out, w, q)
    result_wrapped = []
    for i in range(0,len(result),1):
        result_wrapped.append(((-w)**(-i) * result[i]) %q)
    return result_wrapped

def find_modulus(length, minimum):
    if length < 1 or minimum < 1:
        raise ValueError()
    start = (minimum - 1 + length - 1) // length
    for i in itertools.count(max(start, 1)):
        n = i * length + 1
        assert n >= minimum
        if is_prime(n):
            return n

def find_primitive_root(degree, totient, mod):
    if not (1 <= degree <= totient < mod):
        raise ValueError()
    if totient % degree != 0:
        raise ValueError()
    gen = Zmod(mod).multiplicative_generator()
    root = pow(gen, totient // degree, mod)
    assert (0 <= int(root) < int(mod)), "root outside bounds"
    return root

def test_ntt_forward():
    actual = ntt_w([6, 0, 10, 7, 2], 3, 11)
    expect = [3, 7, 0, 5, 4]
    print(expect, actual)

def test_ntt_backward():
    actual = inverse_ntt_w([3, 7, 0, 5, 4], 3, 11)
    expect = [6, 0, 10, 7, 2]
    print(expect, actual)

def test_ntt_radix_2():
    input = [6, 0, 10, 7, 1, 2, 3, 1]
    maxval = max(val for val in input)
    mod = find_modulus(len(input), maxval+1)
    root = find_primitive_root(len(input), mod-1, mod)
    actual = ntt_w(input, root, mod)
    ntt_radix_2(input, root, mod)
    print(input, actual)

def test_ntt_radix_3():
    input = [6, 0, 10, 7, 1, 2, 3, 1, 4]
    maxval = max(val for val in input)
    mod = find_modulus(len(input), maxval+1)
    root = find_primitive_root(len(input), mod-1, mod)
    actual = ntt_w(input, root, mod)
    ntt_radix_3(input, root, mod)
    print(input, actual)

def test_transform_radix2_vs_naive():
    TRIALS = 300
    for _ in range(TRIALS):
        veclen = 2**randint(1,8)
        maxval = randint(1,100) + 1
        vec = [randint(1,maxval + 1) for _ in range(veclen)]
        mod = find_modulus(len(vec), maxval+1)
        root = find_primitive_root(len(vec), mod-1, mod)
        temp = ntt_w(vec, root, mod)
        ntt_radix_2(vec, root, mod)
        assert (temp == vec), "there are errors"

def test_ntt_multiply():
    input1 = [3, 7, 0, 5, 4]
    input2 = [1, 0, 2, 3, 7]
    maxval = max(val for val in itertools.chain(input1, input2))
    minmod = maxval**2 * len(input1) + 1
    q = find_modulus(len(input1), minmod)
    w = find_primitive_root(len(input1), q-1, q)
    actual = ntt_multiply(input1, input2, w % q, q)
    R = PolynomialRing(GF(q),'y')
    y = R.gens()[0]
    n = w.multiplicative_order()
    I = R.quotient(y**n-1)
    x = I.gens()[0]
    f1 = I(input1)
    f2 = I(input2)
    expect = (f1*f2).lift().coeffs()
    print(expect, actual)
    return True

def test_ntt_multiply_wrapped():
    input1 = [3, 7, 0, 5, 4]
    input2 = [1, 0, 2, 3, 7]
    maxval = max(val for val in itertools.chain(input1, input2))
    minmod = maxval**2 * len(input1) + 1
    q = find_modulus(len(input1), minmod)
    w = find_primitive_root(len(input1), q-1, q)
    actual = ntt_multiply_wrapped(input1, input2, w % q, q)
    R = PolynomialRing(GF(q),'y')
    y = R.gens()[0]
    n = w.multiplicative_order()
    I = R.quotient(y**n+1)
    x = I.gens()[0]
    f1 = I(input1)
    f2 = I(input2)
    expect = (f1*f2).lift().coeffs()
    print(expect, actual)
    return True

def test_ntt_multiply_3():
    input1 = [3, 7, 0, 5, 4, 1, 6, 2, 8]
    input2 = [1, 0, 2, 3, 7, 1, 1, 0, 5]
    n = len(input1)
    assert n % 3 == 0, "input length is not a multiple of 3"
    k = int(n/3)
    maxval = max(val for val in itertools.chain(input1, input2))
    minmod = maxval**2 * len(input1) + 1
    q = find_modulus(len(input1), minmod)
    w = find_primitive_root(len(input1), q-1, q)
    actual = ntt_multiply(input1, input2, w % q, q)
    R = PolynomialRing(GF(q),'y')
    y = R.gens()[0]
    I = R.quotient(y**(2*k)+y**k+1)
    x = I.gens()[0]
    f1 = I(input1)
    f2 = I(input2)
    expect = (f1*f2).lift().coeffs()
    fa = I(actual)
    actuala = fa.lift().coeffs()
    print(expect, actuala)
    return True

def test_simple_convolution():
    mod = 673
    root = 326
    vec0 = ntt_w([4, 1, 4, 2, 1, 3, 5, 6], root, mod)
    vec1 = ntt_w([6, 1, 8, 0, 3, 3, 9, 8], root, mod)
    vec2 = [(x * y % mod) for (x, y) in zip(vec0, vec1)]
    actual = inverse_ntt_w(vec2, root, mod)
    expect = [123, 120, 106, 92, 139, 144, 140, 124]
    R = PolynomialRing(GF(mod),'y')
    y = R.gens()[0]
    n = Zmod(mod)(root).multiplicative_order()
    I = R.quotient(y**n-1)
    x = I.gens()[0]
    f1 = I([4, 1, 4, 2, 1, 3, 5, 6])
    f2 = I([6, 1, 8, 0, 3, 3, 9, 8])
    expect2 = (f1*f2).lift().coeffs()
    print(expect, expect2, actual)
