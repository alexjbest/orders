load "orderofindex.m";

R<x> := PolynomialRing(RationalField());

K := NumberField(x^3 + x - 1);
Zk := Integers(K);
"Discriminant:", Factorisation(Discriminant(K));

//O:=sub<Zk|[2*Zk.2,4*Zk.3]>;

n := 0;
SO := 1;
C := 0;
repeat
    randElt := [K!Random(Zk,8) : i in [1..2]];
    O := Order(randElt);
    I := Index(Zk,O);
    N := Norm(idealOfO(Zk,Conductor(O)));
    if O ne sub<Zk|Basis(Conductor(O))> then
        "Order:", O;
        "Conductor:", idealOfO(Zk,Conductor(O));
        "Norm of conductor:", N;
        "Index:", I;
        C := C + 1;
        SO := O;
    end if;
    n := n + 1;
until n eq 100;
C;
