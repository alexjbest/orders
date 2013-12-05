function overOrder(I) // Returns the smallest order containing I
    return Order(Basis(I,FieldOfFractions(Order(I))));
    //return Order(TwoElement(I));
end function;

function idealOfO(O,I)
    return ideal<O|[O!b : b in Basis(I,FieldOfFractions(O))]>;
end function;

function orderOfIndex(K, d)
    "Finding orders:";
    primeIdeals := AssociativeArray();
    Zk := Integers(K);

    for pe in Factorisation(d^2) do
        p := pe[1];
        e := pe[2];
        "Prime:", p;
        primeIdeals[p] := AssociativeArray();

        for i in [1..e] do
            primeIdeals[p][i] := [];
        end for;

        for PE in Factorisation(ideal<Zk|p>) do
            P := PE[1];
            E := PE[2];
            v := Valuation(Norm(P),p);
            if v in Keys(primeIdeals[p]) then // Otherwise the exponent is too high anyway
                Append(~primeIdeals[p][Valuation(Norm(P),p)],P); // Add P to the set of prime ideals of norm p^e
            end if;
        end for;
        for x in Keys(primeIdeals[p]) do x, primeIdeals[p][x]; end for;
    end for;

    "Partitions:";

    possConds := [ ideal<Zk|1> ];
    for pe in Factorisation(d^2) do
        p := pe[1];
        e := pe[2];
        p,"^",e;
        possFacts := [* *];

        for part in RestrictedPartitions(e, Keys(primeIdeals[p])) do
            part;
            curFact := [ideal<Zk|1>];
            for i in part do
                newFact := [* *];
                for P in primeIdeals[p][i] do
                    for f in curFact do
                        Append(~newFact, f*P);
                    end for;
                end for;
                curFact := newFact;
            end for;
            possFacts cat:= curFact;
        end for;
        curConds := [];
        for P in possFacts do
            for f in possConds do
                Append(~curConds, f*P);
            end for;
        end for;
        possConds := curConds;
    end for;

    out := [];
    
    "Trying conductors:";
    for f in possConds do
        f,Norm(f);
        if Conductor(overOrder(f)) eq f then
            Of := overOrder(f);
            if Index(Zk,Of) eq d then
                Append(~out, Of);
            end if;
        end if;
    end for;

    return out;
end function;

R<x> := PolynomialRing(RationalField());
K := NumberField(x^6 + 20432*x^5+ 18*x^2 +  5*x+ 108);
//K := NumberField([x^2 - 3, x^2 -7]);
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

result := orderOfIndex(K,33);
//if #result gt 0 then
//    result[1];
//end if;
