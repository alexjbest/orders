def conductor(order):
    """
    Return the conductor of the order.
    """
    
    K = order.fraction_field()
    R = Integers()
    ZK = K.maximal_order()
    omega = ZK.basis()
    n = order.rank()
    M = matrix(R.fraction_field(), nrows = n*n, ncols = n)
    
    d = 1
    for i in range(n):
        for j in range(n):
            coords = order.coordinates(ZK.gen(i)*ZK.gen(j))
            for k in range(n):
                M[j*n + k, i] = coords[k]
                d = lcm(d, coords[k].denominator())
    
    H = (d * M).change_ring(R).hermite_form(include_zero_rows = False)
    
    return K.ideal(list(vector(omega) * d * H.inverse()))

def orders_of_index(O,I):
    """
    Returns a list of orders with the given index.
    """
    K = O.fraction_field()
    ZK = K.maximal_order()
    R = Integers() # This may need to be different when relative orders are looked for.

    # We find all ideals of norm dividing I^2, which are all possibly conductors of our order.
    possible_conductors = []
    for d in divisors(I):
        possible_conductors += ideals_of_norm(K, I*d)
    
    print possible_conductors
    
    orders = []
    for f in possible_conductors:
        #print f
        #print f.basis()
        cur_orders = orders_with_conductor_and_index(f, I)
        if cur_orders:
            orders += cur_orders
        print str(len(cur_orders)) + " order(s) found with right index."
    return orders

def ideals_of_norm(K,N):
    # We use the factorisation of prime ideals to find all ideals with given norm N
    primes = dict() # primes[p][i] will contain all prime ideals with norm p^i.
    for (p,e) in N.factor():
        #print str(p)+"^"+str(e)
        primes[p] = dict()
        for i in range(1,e+1):
            primes[p][i] = []

        for (P,E) in K.factor(p): # for P in K.prime_factors(p) ?
            v = valuation(P.norm(), p)
            # replacing the condition below with v <= e would probably be faster, if it turns out we cannot rule out any k <= e as exponents.
            if v in primes[p]: # Otherwise the ideal is too small to be of use to us.
                primes[p][v].append(P)

    #print primes
    ideals = [K.ideal(1)]
    for (p,e) in N.factor():
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

        current_ideals = []
        for P in possible_factorisations:
            for f in ideals:
                current_ideals.append(f*P)
        ideals = current_ideals
    return ideals

def orders_with_conductor_and_index(f,I):
    K = f.number_field()
    Zk = K.maximal_order()
    naive_O = K.order(f.basis()) # Use ring_generators here?!
    if conductor(naive_O) == f:
        if naive_O.index_in(Zk) == I:
            return [naive_O]
        else: # Still not clever enough!
            orders = []
            quo = Zk.free_module().quotient(naive_O.free_module())
            r = quo.cardinality() / I
            for a in quo:
                if a.additive_order() == r:
                    O = K.order(naive_O.gens() + [Zk(a.lift())])
                    if O.index_in(Zk) == I:
                        if not O in orders:
                            orders.append(O)
            return orders
    else:
        pass
        #raise Exception("AAAH")
    return []

def cocyclic_orders_of_index(O,I):
    """
    Returns a list of orders with the given index.
    """
    K = O.fraction_field()
    ZK = K.maximal_order()
    R = Integers() # This may need to be different when relative orders are looked for.

    possible_conductors = ideals_of_norm(K, I*I)
    
    print possible_conductors
    
    orders = []
    for f in possible_conductors:
        q = ZK.free_module().quotient(f.free_module())
        if q.ngens() == 2:
            if q.gens()[0].order() == I: # ????
                orders.append(cocyclic_order_with_conductor(f))
        #print f
        #print f.basis()
    return orders

def cocyclic_order_with_conductor(f):
    K = f.number_field()
    Zk = K.maximal_order()
    O = K.order(f.basis()) # Use ring_generators here?!
    if conductor(O) == f:
        return O
    raise Exception("Order was not cocyclic.")
    return []
