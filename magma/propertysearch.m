load "orderofindex.m";

R<x> := PolynomialRing(RationalField());

//K := NumberField(x^5+ 18*x^2 +  5*x+ 178);
K1 := NumberField(x^2- 5);
K2 := NumberField(x^2- 7);
K := Compositum(K1,K2);

DefiningPolynomial(K);
Zk := Integers(K);

procedure test()
    repeat
        randElt := [K!Random(Zk,5) : i in [1..1]];
        O := Order(randElt);
        I := Index(Zk,O);
    until I ne 1;

    "Discriminant:", Factorisation(Discriminant(K));
    "Index:", I;

    result := orderOfIndex(K,I);
    for o in result do
        #(Factorisation(idealOfO(Zk,Conductor(o))));
    end for;
end procedure;
test();
