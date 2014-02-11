def random_test(K, n):
    Zk = K.maximal_order()
    for i in range(n):
        O = K.order([Zk.random_element() for j in range(2)])
        ords = order_of_index(Zk, O.index_in(Zk))
        if not O in ords:
            print "EROOOOOOOOOOOOOOR"
            print O
            print O.basis()
