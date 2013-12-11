load "orderofindex.m";

R<x> := PolynomialRing(RationalField());
K := NumberField(20432*x^5+ 18*x^2 +  5*x+ 108);
DefiningPolynomial(K);
Zk := Integers(K);
I := 1;
repeat
    randElt := [K!Random(Zk,6) : i in [1..2]];
    O := Order(randElt);
    I := Index(Zk,O);
until I ne 1;

"Discriminant:", Factorisation(Discriminant(K));
"Order:", O;
"Index:", I;
"Conductor:", TwoElement(idealOfO(Zk,Conductor(O))), idealOfO(Zk,Conductor(O));
"Index of conductor:",Norm(idealOfO(Zk,Conductor(O)));
"Factorisation:",Factorisation(idealOfO(Zk,Conductor(O)));

result := orderOfIndex(K,11);
if #result gt 0 then
    result;
else;
    "No orders with given index";
end if;
