load "orderofindexrelative.m";

R<y> := PolynomialRing(RationalField());
K := NumberField(y^5+ 18*y^2 +  5*y+ 178);
P<x> := PolynomialRing(K);
L := NumberField(x^4+67*x^2 + 17);

L;
DefiningPolynomial(L);
ZL := Integers(L);
repeat
    RandElts := [L!Random(ZL,6) : i in [1..2]];
    O := sub<ZL|RandElts>;
    I := Index(ZL,O);
until I ne ideal<ZL|1>;

"Order:", O;
"Discriminant:", Discriminant(L);
"Index:", I;
"Bas", Basis(I,FieldOfFractions(O));
"Conductor:", Conductor(O);
"Index of conductor:",Norm(idealOfO(ZL,Conductor(O)));
"Factorisation:",Factorisation(idealOfO(ZL,Conductor(O)));

result := orderOfIndex(L,4);
if #result gt 0 then
    result;

    //O := result[2];
    //Basis(O, K);
    //cond := idealOfO(Zk,Conductor(O));
else;
    "No orders with given index";
end if;

