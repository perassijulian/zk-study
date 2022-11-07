from channel import Channel
from field import FieldElement
from merkle import MerkleTree
from polynomial import interpolate_poly, Polynomial

def FriCommit(cp, domain, cp_eval, cp_merkle, channel):
    fri_polys = [cp]
    fri_doms = [domain]
    fri_layers = [cp_eval]
    merkles = [cp_merkle]
    while fri_polys[-1].degree() > 0:
        alpha = channel.receive_random_field_element()
        next_poly, next_dom, next_layer = next_fri_layer(fri_polys[-1], fri_doms[-1], alpha)
        fri_polys.append(next_poly)
        fri_doms.append(next_dom)
        fri_layers.append(next_layer)
        merkles.append(MerkleTree(next_layer))
        channel.send(merkles[-1].root)
    channel.send(str(fri_polys[-1].poly[0]))
    return fri_polys, fri_doms, fri_layers, merkles

def next_fri_domain(domain):
    return [x ** 2 for x in domain[:len(domain) // 2]]

def next_fri_polynomial(poly, alpha):
    odd_coefficients = poly.poly[1::2]
    even_coefficients = poly.poly[::2]
    odd = Polynomial(odd_coefficients).scalar_mul(alpha)
    even = Polynomial(even_coefficients)
    return odd + even


def next_fri_layer(poly, dom, alpha):
    next_poly = next_fri_polynomial(poly, alpha)
    next_dom = next_fri_domain(dom)
    next_layer = [next_poly.eval(x) for x in next_dom]
    return next_poly, next_dom, next_layer

print('Making Square Fibonacci list (y)...')
t = [FieldElement(1), FieldElement(3141592)]
while len(t) < 1023:
    t.append(t[-2] * t[-2] + t[-1] * t[-1])

print('Making generator (x)...')
subgroup_size = 1023 + 1
multiplicative_group_size = FieldElement.k_modulus - 1
k = multiplicative_group_size / subgroup_size
g = FieldElement.generator() ** k
points = [g ** i for i in range(1024)]

print('Making extended domain generator (x)...')
eval_domain_size = subgroup_size * 8
k = multiplicative_group_size / eval_domain_size
h_gen = FieldElement.generator() ** k
h = [h_gen ** i for i in range(eval_domain_size)]
domain = [FieldElement.generator() * x for x in h]

print('Interoplating values (y) in domain (x)...')
p = interpolate_poly(points[:-1], t)

print('Evaluating interpolated poly in extended domain...')
ev = [p.eval(d) for d in domain]

print('Getting merkle tree and commiting it...')
mt = MerkleTree(ev)
ch = Channel()
ch.send(mt.root)

print('Construction constrain polynomials...')
numer0 = p - Polynomial([FieldElement(1)])
denom0 = Polynomial.gen_linear_term(FieldElement(1))
q0, r0 = numer0.qdiv(denom0)
numer1 = p - Polynomial([FieldElement(2338775057)])
denom1 = Polynomial.gen_linear_term(points[1022])
q1, r1 = numer1.qdiv(denom1)
inner_poly0 = Polynomial([FieldElement(0), points[2]])
final0 = p.compose(inner_poly0)
inner_poly1 = Polynomial([FieldElement(0), points[1]])
composition = p.compose(inner_poly1)
final1 = composition * composition
final2 = p * p
numer2 = final0 - final1 - final2
coef = [FieldElement(1)] + [FieldElement(0)] * 1023 + [FieldElement(-1)]
numerator_of_denom2 = Polynomial(coef)
factor0 = Polynomial.gen_linear_term(points[1021])
factor1 = Polynomial.gen_linear_term(points[1022])
factor2 = Polynomial.gen_linear_term(points[1023])
denom_of_denom2 = factor0 * factor1 * factor2
denom2, r_denom2 = numerator_of_denom2.qdiv(denom_of_denom2)
q2, r2 = numer2.qdiv(denom2)

print('Making composition polynomial...')
cp0 = q0.scalar_mul(ch.receive_random_field_element())
cp1 = q1.scalar_mul(ch.receive_random_field_element())
cp2 = q2.scalar_mul(ch.receive_random_field_element())
cp = cp0 + cp1 + cp2

print('Evaluation composition polynomial...')
cp_ev = [cp.eval(d) for d in domain]

print('Getting MT and commiting CP...')
cp_mt = MerkleTree(cp_ev)
ch.send(cp_mt.root)

print('FRI algorithm...')
fri_polys, fri_domains, fri_layers, fri_merkles = FriCommit(cp, domain, cp_ev, cp_mt, ch)
print(ch.proof)
assert len(fri_layers) == 11, f'Expected number of FRI layers is 11, whereas it is actually {len(fri_layers)}.'
assert len(fri_layers[-1]) == 8, f'Expected last layer to contain exactly 8 elements, it contains {len(fri_layers[-1])}.'
assert all([x == fri_layers[-1][-1] for x in fri_layers[-1]]), f'Expected last layer to be constant.'
assert fri_polys[-1].degree() == 0, 'Expacted last polynomial to be constant (degree 0).'
# assert fri_merkles[-1].root == '1c033312a4df82248bda518b319479c22ea87bd6e15a150db400eeff653ee2ee', 'Last layer Merkle root is wrong.'
# assert ch.state == '61452c72d8f4279b86fa49e9fb0fdef0246b396a4230a2bfb24e2d5d6bf79c2e', 'The channel state is not as expected.'
print('Success!')