x = QQ['x'].gens()[0]

E = EllipticCurve([0,0,1,-1,0])
#P = x**3 - x
#E = EllipticCurve([0,-1,-1,0,0]) # DOES NOT WORK
P = x**3 + E.a2()*x**2 + E.a4()*x + E.a6()
Q = E.a1()*x + E.a3()

R = 4*P + Q**2
R0 = R/4
f = 16*R.subs({x:x/4})
print f.disc().factor()
print R.disc().factor()

print 2**(8)*E.discriminant() == f.disc()

K = NumberField(f,'a')

disc = E.discriminant() # 37
index = sqrt(Integer((2**8 *disc)/K.disc()).abs())
print index
ords = orders_of_index(K.maximal_order(), index)

curves = []
for o in ords:
    print o.index_in(K.maximal_order())
    if len(o.ring_generators()) == 1: # Monogenic
        f = o.ring_generators()[0].minpoly().subs({x:4*x})/16
        print f
        coeffs = f.coeffs()
        #coeffs = o.ring_generators()[0].minpoly().subs({x:4*x})).coeffs()
        curves.append(EllipticCurve([0,coeffs[2]/4,0,coeffs[1]/4,coeffs[0]/4]).minimal_model())
