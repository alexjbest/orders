def conductor(order):
    """
    Return the conductor of the order.
    """
    K = order.fraction_field()
    R = Integers()
    omega = K.maximal_order().basis()
    n = len(omega)
    M = matrix(R.fraction_field(), nrows = n*n, ncols = n)

    d = 1
    for i in range(n):
        for j in range(n):
            coords = order.coordinates(omega[i]*omega[j])
            for k in range(n):
                M[j*n + k, i] = coords[k]
                d = lcm(d,coords[k].denominator())
    print M
    print d
    
    H = (M).hermite_form(include_zero_rows = False)
    print H
    print d*H.inverse()
    
    beta = vector(omega)*d*H.inverse()
    print beta
    return K.ideal(list(beta))

def order_of_index(O,I):
    primes = dict() # primes[p][i] will contain all prime ideals with norm p^i.
    K = O.fraction_field()
    ZK = K.maximal_order()
    R = Integers() # This may need to be different when relative orders are looked for.

    for (p,e) in (I*I).factor():
        #print str(p)+"^"+str(e)
        primes[p] = dict()
        for i in range(1,e+1):
            primes[p][i] = []

        for (P,E) in K.ideal(p).factor():
            v = valuation(P.norm(), p)
            if v in primes[p]:
                primes[p][v].append(P)

    #print primes
    # We use the factorisation of prime ideals to find all ideals of norm I^2, which are all possibly conductors of O
    possible_conductors = [K.ideal(1)]
    for (p,e) in (I*I).factor():
        possible_factorisations = []
        for partition in Partitions(e, parts_in = primes[p].keys()):
            #print partition
            current_factorisation = [K.ideal(1)]
            for i in partition:
                new_factorisation = []
                for P in primes[p][i]:
                    for f in current_factorisation:
                        new_factorisation.append(f*P)
                current_factorisation = new_factorisation
            possible_factorisations += new_factorisation

        current_conductors = []
        for P in possible_factorisations:
            for f in possible_conductors:
                current_conductors.append(f*P)
        possible_conductors = current_conductors
        #print possible_conductors

    orders = []
    for f in possible_conductors:
        #print f
        #print f.basis()
        Of = K.order(f.basis()) # Use ring_generators here?!
        if conductor(Of) == f:
            if Of.index_in(ZK) == I:
                print Of
                orders.append(Of)
            else:
                raise Exception(str(Of) + " has index " + str(Of.index_in(ZK)) + " but " + str(f))
    
    return orders
            


