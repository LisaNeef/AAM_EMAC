

Y=y0:1:yf
M=1:1:12
D=day0:1:dayf
H=0:1:24
MN=0:1:60
S=0:1:60

V=[Y M D H MN S];
N=datenum(V)
datetick('x','YYYY')
