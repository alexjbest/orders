x = QQ['x'].gens()[0]
K = NumberField(x**3 - x**2 + 1,'a')
counts = []
for i in range(1,60):
    counts.append(len(orders_of_index_via_hnf(K.maximal_order(), i)))
