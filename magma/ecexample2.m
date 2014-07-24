load "orderofindex.m";

R<x> := PolynomialRing(RationalField());

K := NumberField(x^3 -x^2 + x + 1);

Zk := Integers(K);
Basis(K);

IndexFormEquation(Zk,968);
