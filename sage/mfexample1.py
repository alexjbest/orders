x = QQ['x'].gens()[0]

E = EllipticCurve([0,-1,-1,0,0]) # DOES NOT WORK
print E.discriminant()

f = x**3 - x**2 - 24*x + 32
K = NumberField(f,'a')

disc = 3#E.discriminant()
index = sqrt((2**8 *disc)/K.disc())
print index
ords = orders_of_index(K.maximal_order(), index)

curves = []
for o in ords:
    if len(o.ring_generators()) == 1: # Monogenic
        coeffs = o.ring_generators()[0].minpoly().coeffs()
        curves.append(EllipticCurve([0,coeffs[2],0,coeffs[1],coeffs[0]]).minimal_model())
