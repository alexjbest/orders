def conductor(order):
    """
    Return the conductor of the order.
    """
    K = order.fraction_field()
    R = Integers()
    omega = K.maximal_order().basis()
    tau = order.basis()
    n = len(omega)
    b = [[0]*n for i in range(n)]
    M = matrix(R, nrows = n*n, ncols = n)

    for i in range(n):
        for j in range(n):
            coords = order.coordinates(omega[i]*omega[j])
            for k in range(n):
                M[j*n + k, i] = coords[k]
    print M
        
    #d = gcd()
    #H = (d*M).hermite_form()
    #beta = omega*d*H.inverse()
    #return K.ideal(beta)
