load "orderofindex.m";

R<x> := PolynomialRing(RationalField());

//K := NumberField(x^3+ 18*x^2 +  5*x+ 178);
K1 := NumberField(x^2- 5);
K2 := NumberField(x^2- 7);
K := Compositum(K1,K2);

DefiningPolynomial(K);
Zk := Integers(K);

procedure test()
    repeat
        randElt := [K!Random(Zk,3) : i in [1..3]];
        O := Order(randElt);
        I := Index(Zk,O);
    until I ne 1;
    
    result := orderOfIndex(K,I);

    if #(result) eq 0 then
        ">>>>>BUUUG";
        "Order:", O;
        "Conductor", Conductor(O);
        "Index:", I;
    end if;
end procedure;

for i := 0 to 50 do
    test();
end for;
