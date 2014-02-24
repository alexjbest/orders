x = QQ['x'].gens()[0]

#E = EllipticCurve([0,0,1,-1,0])
#P = x**3 - x
E = EllipticCurve([0,-1,-1,0,0]) # DOES NOT WORK
P = x**3 + E.a2()*x**2 + E.a4()*x + E.a6()
Q = E.a1()*x + E.a3()

R = 4*P + Q**2
R0 = R/4
f = 16*R.subs({x:x/4})

print E.discriminant() == -2**(-4)*R.disc()

K = NumberField(f,'a')

disc = E.discriminant() # 37
index = sqrt((2**8 *disc)/K.disc())
print index
ords = orders_of_index(K.maximal_order(), index)

curves = []
for o in ords:
    if len(o.ring_generators()) == 1: # Monogenic
        coeffs = o.ring_generators()[0].minpoly().coeffs()
        curves.append(EllipticCurve([0,coeffs[2],0,coeffs[1],coeffs[0]]).minimal_model())
