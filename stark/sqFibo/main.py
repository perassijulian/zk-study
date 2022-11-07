from field import FieldElement
from polynomial import Polynomial, interpolate_poly, X
from hashlib import sha256
from channel import serialize
from merkle import MerkleTree
from channel import Channel

def FriCommit(cp, domain, cp_eval, cp_merkle, channel):
    fri_polys = [cp]
    fri_domains = [domain]
    fri_layers = [cp_eval]
    fri_merkles = [cp_merkle]
    while fri_polys[-1].degree() > 0:
        print(fri_merkles[-1].root)
        beta = channel.receive_random_field_element()
        next_poly, next_domain, next_layer = next_fri_layer(fri_polys[-1], fri_domains[-1], beta)
        fri_polys.append(next_poly)
        fri_domains.append(next_domain)
        fri_layers.append(next_layer)
        fri_merkles.append(MerkleTree(next_layer))
        channel.send(fri_merkles[-1].root)
    channel.send(str(fri_polys[-1].poly[0]))
    return fri_polys, fri_domains, fri_layers, fri_merkles

def next_fri_layer(poly, domain, beta):
    next_poly = next_fri_polynomial(poly, beta)
    next_domain = next_fri_domain(domain)
    next_layer = [next_poly(d) for d in next_domain]
    return next_poly, next_domain, next_layer

def next_fri_polynomial(poly,  beta):
    odd_coefficients = poly.poly[1::2]
    even_coefficients = poly.poly[::2]
    odd = beta * Polynomial(odd_coefficients)
    even = Polynomial(even_coefficients)
    return odd + even

def next_fri_domain(fri_domain):
    return [x**2 for x in fri_domain[:len(fri_domain) // 2]]

def CP_eval(cp):
    return[cp(e) for e in eval_domain]

def get_CP(channel):
    alpha0 = channel.receive_random_field_element()
    alpha1 = channel.receive_random_field_element()
    alpha2 = channel.receive_random_field_element()
    return alpha0*p0 + alpha1*p1 + alpha2*p2

# TRACE
a = [FieldElement(1), FieldElement(3141592)]
while len(a) < 1023:
    a.append(a[-2]*a[-2] + a[-1]*a[-1] )

# testing
assert len(a) == 1023, 'The trace must consist of exactly 1023 elements.'
assert a[0] == FieldElement(1), 'The first element in the trace must be the unit element.'
for i in range(2, 1023):
    assert a[i] == a[i - 1] * a[i - 1] + a[i - 2] * a[i - 2], f'The FibonacciSq recursion rule does not apply for index {i}'
assert a[1022] == FieldElement(2338775057), 'Wrong last element!'
print('Trace testing success!')

# MAKING GENERATOR
subgroup_size = 1023 + 1
multiplicative_group_size = FieldElement.k_modulus - 1
k = multiplicative_group_size / subgroup_size
g = FieldElement.generator() ** k
G = [ g ** i for i in range(subgroup_size)]

# testing that g and G are correct.
assert g.is_order(1024), 'The generator g is of wrong order.'
b = FieldElement(1)
for i in range(1023):
    assert b == G[i], 'The i-th place in G is not equal to the i-th power of g.'
    b = b * g
    assert b != FieldElement(1), f'g is of order {i + 1}'
    
if b * g == FieldElement(1):
    print('Generator testing success!')
else:
    print('g is of order > 1024')

# INTERPOLATING 
f = interpolate_poly(G[:-1], a) #except for G[-1], since a is one element shorter

# MAKING LARGER DOMAIN
eval_domain_size = subgroup_size * 8
k = multiplicative_group_size / eval_domain_size
h = FieldElement.generator() ** k
H = [ h ** i for i in range(eval_domain_size)]
eval_domain = [ FieldElement.generator() * h for h in H]

# testing
assert len(set(eval_domain)) == len(eval_domain)
w = FieldElement.generator()
w_inv = w.inverse()
assert '55fe9505f35b6d77660537f6541d441ec1bd919d03901210384c6aa1da2682ce' == sha256(str(H[1]).encode()).hexdigest(),\
    'H list is incorrect. H[1] should be h (i.e., the generator of H).'
for i in range(8192):
    assert ((w_inv * eval_domain[1]) ** i) * w == eval_domain[i]
print('Testing H list success!')

# EVALUATING IN LARGER DOMAIN
f_eval =  [ f(x) for x in eval_domain ]

# testing against a precomputed hash.
assert '1d357f674c27194715d1440f6a166e30855550cb8cb8efeb72827f6a1bf9b5bb' == sha256(serialize(f_eval).encode()).hexdigest()
print('Evaluating in extended domain success!')

# GETTING ROOT
f_merkle = MerkleTree(f_eval)
assert f_merkle.root == '6c266a104eeaceae93c14ad799ce595ec8c2764359d7ad1b4b7c57a4da52be04'
print('Success!')

# APPLYING FIAT-SHAMIR AND SENDING TO VERIFIER
channel = Channel()
channel.send(f_merkle.root)

# CONSTRUCTION CONSTRAIN POLYNOMIALS
numer0 = f - 1
denom0 = X - 1 
p0 = numer0 / denom0

numer1 = f - 2338775057
denom1 = X - g**1022
p1 = numer1 / denom1

numer2 = f(g**2 * X) - f(g * X)**2 - f**2
denom2 = (X**1024 - 1) / ((X-g**1023) * (X-g**1022) * (X-g**1021))
p2 = numer2 / denom2

# testing polynomials
assert p0(2718) == 2509888982
assert p1(5772) == 232961446
assert p2(31415) == 2090051528
assert p2.degree() == 1023, f'The degree of the third constraint is {p2.degree()} when it should be 1023.'
print('Contraints made correctly!')

# MAKING COMPOSITION POLYNOMIAL AND COMMITING
cp = get_CP(channel)
cp_eval = CP_eval(cp)
cp_merkle = MerkleTree(cp_eval)

channel.send(cp_merkle.root)

# assert cp_merkle.root == 'a8c87ef9764af3fa005a1a2cf3ec8db50e754ccb655be7597ead15ed4a9110f1', 'Merkle tree root is wrong.'
print('Success!')
print(channel.proof)

# FRI FOLDING ALGORITHM
print('STARTING FRI')

next_domain = next_fri_domain(eval_domain)
assert '5446c90d6ed23ea961513d4ae38fc6585f6614a3d392cb087e837754bfd32797' == sha256(','.join([str(i) for i in next_domain]).encode()).hexdigest()
print('Next domain calculated ok')

test_poly = Polynomial([FieldElement(2), FieldElement(3), FieldElement(0), FieldElement(1)])
test_domain = [FieldElement(3), FieldElement(5)]
beta = FieldElement(7)
next_p, next_d, next_l = next_fri_layer(test_poly, test_domain, beta)
assert next_p.poly == [FieldElement(23), FieldElement(7)]
assert next_d == [FieldElement(9)]
assert next_l == [FieldElement(86)]
print('Next layer test calculated ok')

fri_polys, fri_domains, fri_layers, fri_merkles = FriCommit(cp, eval_domain, cp_eval, cp_merkle, channel)
print(channel.proof)
assert len(fri_layers) == 11, f'Expected number of FRI layers is 11, whereas it is actually {len(fri_layers)}.'
assert len(fri_layers[-1]) == 8, f'Expected last layer to contain exactly 8 elements, it contains {len(fri_layers[-1])}.'
assert all([x == FieldElement(-1138734538) for x in fri_layers[-1]]), f'Expected last layer to be constant.'
assert fri_polys[-1].degree() == 0, 'Expacted last polynomial to be constant (degree 0).'
assert fri_merkles[-1].root == '1c033312a4df82248bda518b319479c22ea87bd6e15a150db400eeff653ee2ee', 'Last layer Merkle root is wrong.'
assert channel.state == '61452c72d8f4279b86fa49e9fb0fdef0246b396a4230a2bfb24e2d5d6bf79c2e', 'The channel state is not as expected.'
print('Success!')