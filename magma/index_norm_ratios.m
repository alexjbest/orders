load "orderofindex.m";
R<x> := PolynomialRing(RationalField());

//K := NumberField(x^5+ 26*x^3+  5*x + 1);

K5 := NumberField(x^2- 5);
K13 := NumberField(x^2- 13);
K := Compositum(K5,K13);
Zk := Integers(K);

DefiningPolynomial(K);
n := 0;
SO := 1;
C := 0;
repeat
    randElt := [K!Random(Zk,8) : i in [1..2]];
    O := Order(randElt);
    //"Conductor:", idealOfO(Zk,Conductor(O));
    N := Norm(idealOfO(Zk,Conductor(O)));
    //"Norm of conductor:", N;
    I := Index(Zk,O);
    //"Index:", I;
    if I*I ne N then
        C := C + 1;
        SO := O;
        //N;
        //I^2;
        "R:", I^2/N;
        if (I*I)/N le 1 then
            ">>>>>>>>>>>>>";
        end if;
    end if;
    n := n + 1;
until n eq 500;
C;
//O := Order([1, 4*(K.1^3) + K.1, 4*(K.1^3) + K.1^2, 7*(K.1^3)]);
