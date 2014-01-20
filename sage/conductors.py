def conductor(order):
    """
    Return the conductor of the order.
    """
    print "hi"
    K = order.fraction_field()
    omega = K.maximal_order().basis()
    tau = order.basis()
    n = len(omega)
    for i in range(n):
        for j in range(n):
            print order.coordinates(omega[i]*omega[j])
    #d = gcd()
    #H = (d*M).hermite_form()
    #beta = omega*d*H.inverse()
    #return K.ideal(beta)
