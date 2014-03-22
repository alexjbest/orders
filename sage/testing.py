def random_test(K, n):
    Zk = K.maximal_order()
    for i in range(n):
        O = K.order([Zk.random_element() for j in range(2)]) # Construct a random order
        ords = orders_of_index(Zk, O.index_in(Zk)) # Try and find it
        if not O in ords:
            print "EROOOOOOOOOOOOOOR"
            print O
            print O.basis()

def problem_ideals(K, N):
    Zk = K.maximal_order()
    ideals = K.ideals_of_bdd_norm(N)
    for In in ideals.values():
        for I in In:
            O = K.order(I.basis())
            if not conductor(O) == I:
                print I
