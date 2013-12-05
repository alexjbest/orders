n:=10;
a := [0 : i in [1..n]];
print a;
k := 1;
a[1] := 0;
a[2] := n;
while k gt 0 do
x := a[k] + 1;
y := a[k+1] - 1;
k -:= 1;
while x le y and k lt n - 1 do
a[k+1] := x;
//y -:= x;
k +:= 1;
end while;
a[k+1] := x + y;
print a; //do something with first k + 1 of a here
end while;
