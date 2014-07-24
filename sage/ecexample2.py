def findec(K,discs):
    curves = []
    for disc in discs:
        print "Discriminant: " + str(disc)
        index = Integer(sqrt(2**8*disc/K.disc()))
        print "Index: " + str(index)
        ords = orders_of_index_via_hnf(K.maximal_order(), index)	
        print "Orders found: " + str(len(ords))
        c = 1
        for o in ords:
            print "[" + str(c) + "/" + str(len(ords)) + "]" + str(o.index_in(K.maximal_order()))
            c += 1
            print o.ring_generators()
            monogenic = False #(len(o.ring_generators()) == 1)
            o1 = o.gen(1)
            o2 = o.gen(2)
            al = o.ring_generators()[0]
            gens = []
            for i in range(-20,20):
                for j in range(-20,20):
                    al = i*o1 + j*o2
                    if al == 0:
                        continue
                    if K.order(al) == o:
                        monogenic = True
                        gens.append(al)
            print al
            if monogenic: # len(o.ring_generators()) == 1: # Monogenic
                for g in gens:
                    pol = g.minpoly().subs({x:4*x})/(4**3)
                    coeffs = pol.coeffs()
                    E = EllipticCurve([0,coeffs[2],0,coeffs[1],coeffs[0]])
                    if not E in curves:
                        curves.append(E)
        print curves
    return curves

x = QQ['x'].gens()[0]
f = x**3 - x**2 + x + 1
K = NumberField(f,'a')
discs = [-11**5]
