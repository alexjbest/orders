def conductor(order): # Seems correct now
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

    # We find all ideals of norm dividing I^2 and divisible by I, which are all possibly conductors of our order.
    possible_conductors = []
    for d in divisors(I):
        possible_conductors += ideals_of_norm(K, I*d)
    
    print possible_conductors
    
    orders = []
    c = 1
    for f in possible_conductors:
        #print f
        cur_orders = orders_with_conductor_and_index(f, I)
        for o in cur_orders:
            if not o in orders:
                orders.append(o)
        print "[" + str(c) + "/" + str(len(possible_conductors)) + "] " + str(len(cur_orders)) + " order(s) found with right index."
        c += 1
    return orders

def ideals_of_norm(K,N): # Correct?, could use less memory?
    # We use the factorisation of prime ideals to find all ideals with given norm N
    primes = dict() # primes[p][i] will contain all prime ideals with norm p^i.
    ideals = [K.ideal(1)]
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
    
    #print possible_conductors
    
    orders = []
    for f in possible_conductors:
        q = ZK.free_module().quotient(f.free_module())
        if q.ngens() == 2:
            if q.gens()[0].order() == I: # ????
                O = cocyclic_order_with_conductor(f)
                if not O in orders:
                    orders.append(cocyclic_order_with_conductor(f))
        #print f
        #print f.basis()
    return orders

def cocyclic_order_with_conductor(f):
    K = f.number_field()
    O = K.order(f.basis()) # Use ring_generators here?!
    if conductor(O) == f:
        return O
    raise Exception("Order was not cocyclic.")
    return []

def expressions_as_product(n,k):
    if k == 1:
        yield [n]
    else:
        for d in divisors(n):
            for se in expressions_as_product(n/d,k-1):
                yield se + [d]

def hnf_matrices_with_det(I,n):
    matrices = []
    for exp in expressions_as_product(I,n):
        plain = [[(i == j)*exp[i] for i in range(n)] for j in range(n)]
        mats_for_exp = [plain]
        for i in range(1,n):
            d = exp[i]
            cur_matrices = []
            for j in range(d**i):
                if d == 1:
                    digi = [0 for l in range(i)]
                else:
                    digi = Integer(j).digits(base=d,padto=i)
                for m in mats_for_exp:
                    cur_matrices.append(deepcopy(m))
                    for k in range(i):
                        cur_matrices[-1][k][i] = digi[k]
            mats_for_exp = cur_matrices
        matrices += mats_for_exp
    return matrices

def orders_of_index_via_hnf(O,I):
    n = O.rank()
    K = O.fraction_field()
    orders = []
    for h in hnf_matrices_with_det(I,n-1):
        #if not h.row(0).list() == [1] + [0 for i in range(n-1)]:
        #    continue
        assert O.gen(0) == O(1)
        try:
            o = sage.rings.number_field.order.absolute_order_from_module_generators([O.gen(0)] + [sum([r[i]*O.gen(i+1) for i in range(n-1)]) for r in h],check_integral=False, check_rank=False) # TODO check the checks
            #o = sage.rings.number_field.order.absolute_order_from_module_generators([O(r) for r in h.rows()])
            if o.index_in(O) == I:
                orders.append(o)
            else:
                print "hmm"
        except:
            pass
    return orders

def orders_of_index_iterative(O,I):
    if I == 1:
        return [O]
    r = I.radical()
    ords_r = cocyclic_orders_of_index(O,r)
    ords = []
    for R in ords_r:
        ords += orders_of_index_iterative(R,Integer(I/r))
    return ords
