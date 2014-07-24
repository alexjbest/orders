load "orderofindex.m";

R<x> := PolynomialRing(RationalField());

K := NumberField(x^4 + 5*x + 1);

Zk := Integers(K);
Basis(K);
O := sub<Zk|[Zk!3*K.1^2,Zk!3*K.1^3]>;

"Discriminant:", Factorisation(Discriminant(K));
"Order:", O;
"Index:", Index(Zk,O);
"Conductor:", idealOfO(Zk,Conductor(O));
"Index of conductor:",Norm(idealOfO(Zk,Conductor(O)));
"Factorisation:",Factorisation(idealOfO(Zk,Conductor(O)));

orderOfIndex(K,Index(Zk, O));
