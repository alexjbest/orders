R<x> := PolynomialRing(RationalField());
K := NumberField(x^4 + 5*x + 108);
DefiningPolynomial(K);
Ok := MaximalOrder(K);
//M:= MatrixAlgebra(IntegerRing(),6);
//D:=DiagonalMatrix(M,[2,2,2,2,1,1]);
//D;
//O := Order(Ok,D,1);
//O := sub<Ok|[Random(Ok,2)]>;
//O;

//Conductor(O);

function overOrder(I)
return Order(Basis(I,K));
end function;

function idealOfO(O,I)
return ideal<O|[O!b : b in Basis(I,K)]>;
end function;

print "Random conductor";
randElt := [K!Random(Ok,4) : i in [1..3]];
cond := idealOfO(Ok,Conductor(Order(randElt)));
cond;
Basis(cond, K);

Ord := overOrder(cond);
print "Computed order";
Type(Ord);
Ord;
Ord eq Order(randElt);
print "Computed conductor";
Conductor(Ord);
ncond:=idealOfO(Ok,Conductor(Ord));
ncond;
Basis(ncond, K);
Norm(ncond);
Norm(cond);
ncond eq cond;

//Norm(Conductor(O));
//Index(Ok, Ok*Conductor(O));
