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
