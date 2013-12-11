function overOrder(I,K) // Returns the smallest order containing I
    return Order(Basis(I,K));
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
        possFacts := [ ];

        for part in RestrictedPartitions(e, Keys(primeIdeals[p])) do
            part;
            curFact := [ideal<Zk|1>];
            for i in part do
                newFact := [ ];
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
        Of := overOrder(f,K);
        if Conductor(Of) eq f then
            if Index(Zk,Of) eq d then // Sanity check
                Append(~out, Of);
            else;
                error Error();
            end if;
        end if;
    end for;

    return out;
end function;
