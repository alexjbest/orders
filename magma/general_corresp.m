load "orderofindex.m";

R<x> := PolynomialRing(RationalField());

K := NumberField(x^4+ 5*x^2 + 1);
Zk := Integers(K);
O := sub<Zk|[K.1]>;
Index(Zk,O);
