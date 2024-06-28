/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PARI macros and scripts used for the article 

"Effective multiplicative independence of~$3$ singular moduli"

Yuri Bilu, Sanoli Gun and Emanuele Tron 


Table of contents


A. General scripts

primdivs(m)       prime divisors of m
emax(n)           max m such that ph(m) <= n
vtom(v,m,n)       vector to matrix 
mtov(M)           matrix to vector
prodv(u,v)        vector of products
lexaddone(v,lens) adding one in lexicographic order
chebtheta(x)      Chebyshev's theta
sfome(a,X,m)      integers <=X composed of m distinct primes >=a
intome(a,X,m)     integers n<=X with omega(n)=m and with smallest prime divisor >=a
onlymod(v,m,a)    entries of a vector which are congruent to a mod m


B. Scripts on discriminants 

maxcl(X)            max class number of discriminats <= X in absolute value
                    (used in the proof of Proposition 2.1)
rhotwo(Delta)       2-rank 
fundamx(X)          fundamental discriminants not exceeding X
capitalpsi(m,Delta) h(m^2\Delta)/h(\Delta), see equation (2.4)
twoelfs(D)          finding (almost) 2-elementary Df^2 for a given almost 2-elementary D


C. Scripts for Section 2

twoeldisv(v)             counting 2-elementary  Df^2  with D in a  list
almdisv(v)               counting almost 2-elementary  Df^2  with D in a  list
twoeldisx(X), almdisx(X) counting (almost) 2-elementary discriminants below a given bound
                         (used in the proof of Proposition 2.9)
almdismx(m,X)            counting almost 2-elementary  discriminants  with fundamental discriminant having the given number of prime divisors
                         (used in the proof of Proposition 2.9)


D. Scripts for Section 5

kronmod(Delta,p)           kronecker, but with Delta defined as a residue class
relcln(e,Delta)            Psi as in (2.1), but with Delta defined as a residue class
calq(n)                    Q(n), see Section 2.4
denoms(n,ax,ex,ey)         possible denominators
noone(v)                   removing 1 on the left
exeyez()                   composing Table 2
systems(ex,ey,axone,axtwo) composing Table 3 and solving the systems 
dumbsystems(az)            solving systems from Table 4 
dumbsystemsbis(ex)         solving systems from Table 5



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/ 


/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
A. General scripts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

/*----------------------------------------
prime divisors of m
-------------------------------------------*/

primdivs(m)=

{

my(f,ell,p);

f=factor(m);

f[,1]

}





/*------------------------------------------
max m such that ph(m) <= n
--------------------------------------------*/


emax(n)=
{
my(m,N);

N=max(6,n^2);
m=1;

for (k=1,N,
if(eulerphi(k)<=n,
m=k;);
);
return(m);
}

/*
modtobmod(M)=
lift(M)"\\bmod"M.mod
*/


maxrepeat(v)=

{
my(m,n,mi);
m=1;
n=length(v);

for (i=1,n-1,

mi=1;

for(j=i+1,n,

if(v[i]==v[j],
mi++;
);

);

m=max(m,mi);

);

m
}

/*------------------------------------------
vector to matrix 

input: 
m, n positive integers; 
vector v of length mn

output: matrix mxn made of the entries of v row by row
-----------------------------------------------*/

vtom(v,m,n)=

{


if (length(v) !=m*n,

print ("bad dimensions");

return;

);

matrix(m,n,i,j,v[n*(i-1)+j])


}

/*------------------------------------
matrix to vector

input: a matrix M
output: a vector made of the entries of M row by row
----------------------------------------*/

mtov(M) =

{

my(m,n);

n=length(M);

m=length(M[,1]); 


vector(m*n,k, M[ceil(k/n),k-n*(ceil(k/n)-1)])

}

/*--------------------------------------------------
Vector of products

input: vector u,v of lengths m,n respectively
output: vector of length mn consisting of products u[i]*v[j]

-----------------------------------------------------*/

prodv(u,v)=

{

my(m,n,M);

m=length(u);
n=length(v);

M=matrix(m,n,i,j,u[i]*v[j]);

mtov(M)

}

/*---------------------------------------------------
adding one in lexicographic order

input:
v & lens vectors of positive integers of the same length satisfying 
v[i] \in {1,...,lens[i]}

output: vector next to v in the lexicographic order
if v=lens then the output is v 
-----------------------------------------------------*/

lexaddone(v,lens) = 


{

my(ii,n);

n=length(lens);

if (length(v) !=n,

print("bad data");

return;

);

for (i=1,n,

if (v[i]>lens[i],

print("bad data");

return;

);

);

if (v==lens,

return(v);

);

for (i=1,n,

if (v[i]<lens[i],

ii=i;

);

);

v[ii]++;

for (i=ii+1,n,

v[i]=1;

);

v


}


/*-------------------------------------
Chebyshev's theta
---------------------------------------*/

chebtheta(x)=

{

my(S);

S=0;

forprime(p=2,floor(x), 
S=S+log(p);
);

S

}

/*------------------------------------------------
integers <=X composed of m distinct primes >=a
------------------------------------------------*/

sfome(a,X,m)=

{

my(v);

if (m==0, return([1]););

/* if (a>b&&m>0, return([]);); */

v=[];

forprime(p=a, X^(1/m),
v=concat(v,p*sfome(p+1,X/p,m-1));
);

v

}

/*------------------------------------------------
integers n<=X with omega(n)=m and with smallest prime divisor >=a
------------------------------------------------*/

intome(a,X,m)=

{

my(v,e);

if (m==0, return([1]););


v=[];

forprime(p=a, X^(1/m),

e=floor((log(X/p))/log(2));



if(p^(e+1)<=X, e=e+1;); /* avoiding rounding error */

for(k=1,e,

v=concat(v,p^k*intome(p+1,X/p^k,m-1));
);

);

v

}

/*------------------------------------------------------------------
entries of a vector which are congruent to a mod m
-------------------------------------------------------------------*/

onlymod(v,m,a)=

{

my(k);

k=0;

for(i=1,length(v),

if (v[i]%m==a%m,
k++;
v[k]=v[i];
);

);

v[1..k]

}



/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
B. Scripts on discriminants
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

/*--------------------------------------------------
max class number of discriminats <= X(<10^10) in absolute value
used in the proof of Proposition 2.1
--------------------------------------------------*/

maxcl(X) =

{
my(Y,h);

Y=floor(X/4);

h=1;

for (n=1,Y,

if (qfbclassno(1-4*n)>h,

h= qfbclassno(1-4*n);

);

if (qfbclassno(-4*n)>h,

h= qfbclassno(-4*n);

);

);

return(h);

}


/*--------------------------------------------------
2-rank of Delta
--------------------------------------------------*/

rhotwo(Delta)=

{
my(r,mu,n);

if (Delta%4!=0&&Delta%4!=1, print("bad Delta"); return;); 

if (Delta%2==0, r=omega(-Delta)-1;);

if (Delta%2==1, r=omega(-Delta););

if (Delta%4==1, mu=r;);

if (Delta%4==0, n=-Delta/4;  

if (n%4==3, mu=r;);

if (n%4==1, mu=r+1;);

if (n%4==2, mu=r+1;);

if (n%8==4, mu=r+1;);

if (n%8==0, mu=r+2;);

);

return(mu-1); 


}


/*---------------------------------------------------
Fundamental discriminants not exceeding X
------------------------------------------------------*/

fundamx(X)=

{
my(v, count);

count=0;

v=vector(floor(X));

for (n=1,X,

if (n%4==3&& issquarefree(n),
count++; 
v[count]=-n;
);

if (n%16==4&&issquarefree(n/4),
count++; 
v[count]=-n;
);

if (n%16==8&&issquarefree(n/8),
count++; 
v[count]=-n;
);

);

v[1..count]

}





/*--------------------------------------------
h(m^2\Delta)/h(\Delta), see equation (2.4)
--------------------------------------------*/

capitalpsi(m,Delta)=

{
my(Psi,u,ell,p);

Psi=m;



p= primdivs(m); 

ell=length(p);



for(i=1,ell, 

 
Psi=Psi*(1-kronecker(Delta,p[i])/p[i]);
);


if (Delta==-3&&m>1,
Psi=Psi/3;
);

if (Delta==-4&&m>1,
Psi=Psi/2;
);

Psi

}

/*------------------------------------------------------------------
Finding (almost) 2-elementary Df^2 for a given almost 2-elementary D

input: 
- an almost 2-elementary discriminant D

output: 
- the lists of f such that Df^2 is 2-elementary and such that Df^2 is almost 2-elementary
- the maximal class number for each of the two lists 
------------------------------------------------------------------*/

twoelfs(D)=

{

my(twoels,twoelscount, alms,almscount,hfd,rho);

twoels=vector(8);
alms=vector(20);
twoelscount=0;
almscount=0;
maxclntwoel=0;
maxclnoalm=0;

h=qfbclassno(D); 

if (2^(rhotwo(D)+1)%h !=0, 
print ("bad D");
return();
);

if (D<-3,

fordiv(2^4*3*5*7*17,f,

hfd=h*capitalpsi(f,D);
rho=rhotwo(D*f^2);

if(hfd==2^rho,
twoelscount++;
almscount++;
twoels[twoelscount]=f;
alms[almscount]=f;
maxclntwoel=max(maxclntwoel,hfd);
maxclnoalm=max(maxclnoalm,hfd);
);

if(hfd==2^(rho+1),
almscount++;
alms[almscount]=f;
maxclnoalm=max(maxclnoalm,hfd);
);

);

);

if (D==-3,

for(f=1,23,

hfd=h*capitalpsi(f,D);
rho=rhotwo(D*f^2);

if(hfd==2^rho,
twoelscount++;
almscount++;
twoels[twoelscount]=f;
alms[almscount]=f;
maxclntwoel=max(maxclntwoel,hfd);
maxclnoalm=max(maxclnoalm,hfd);
);

if(hfd==2^(rho+1),
almscount++;
alms[almscount]=f;
maxclnoalm=max(maxclnoalm,hfd);
);

);

);

twoels=twoels[1..twoelscount];
alms=alms[1..almscount];

[twoels,maxclntwoel,alms,maxclnoalm]

}




/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C. Scripts for Section 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

/*-----------------------------------------------------
counting 2-elementary  Df^2  with D in a  list

input: v a list of fundamental discriminants

output:
- the number of 2-elementary  discriminants Df^2 with D  \in v
- the biggest absolute value of these discriminants
- the biggest class number of these discriminants
-------------------------------------------------------*/

twoeldisv(v)=

{
my(D,numfutwoel,maxfutwoel,clnofutwoel,numtwoel,maxtwoel,clnotwoel,u,ltwo);

/*
numfutwoel=0;
maxfutwoel=0;
clnofutwoel=0;
*/

numtwoel=0;
maxtwoel=0;
clnotwoel=0;

for(i=1, length(v),

D=v[i];
rho=rhotwo(D); 

if (D>-10^10, h=qfbclassno(D));
if (D<=-10^10, h=quadclassunit(D).no);

if (h==2^rho,

print(D);

u=twoelfs(D);
ltwo=length(u[1]); 

/*
numfutwoel=numfutwoel++;
maxfutwoel=max(maxfutwoel, abs(D));
clnofutwoel=max(clnotwoel,h);
*/


numtwoel=numtwoel+ltwo;
maxtwoel=max(maxtwoel, abs(D)*u[1][ltwo]^2);
clnotwoel=max(clnotwoel,u[2]); 


);

);

/* [numfutwoel,maxfutwoel,clnofutwoel,numtwoel,maxtwoel,clnotwoel] */

[numtwoel,maxtwoel,clnotwoel]

}






/*-----------------------------------------------------
counting almost 2-elementary  Df^2  with D in a  list

input: v a list of fundamental discriminants

output:
- the number of almost 2-elementary  discriminants Df^2 with D  \in v
- the biggest absolute value of these discriminants
- the biggest class number of these discriminants
-------------------------------------------------------*/

almdisv(v)=

{
my(D,numalm,maxalm,clnoalm,u,lalm);

/* numfualm=0; */

numalm=0;
maxalm=0;
clnoalm=0;

for(i=1, length(v),

D=v[i];
rho=rhotwo(D); 

if (D>-10^10, h=qfbclassno(D));
if (D<=-10^10, h=quadclassunit(D).no);

if (2^(rho+1)%h==0, 
 
print(D);

u=twoelfs(D);
lalm=length(u[3]); 


/*numfualm++;*/

numalm=numalm+lalm;
maxalm=max(maxalm, abs(D)*u[3][lalm]^2);
clnoalm=max(clnoalm,u[4]); 


);

);

/* [numfutwoel,maxfutwoel,clnofutwoel,numtwoel,maxtwoel,clnotwoel] */

[numalm,maxalm,clnoalm]

}

/*-------------------------------------------------------------------
Counting (almost) 2-elementary discriminants below a given bound

input: X>=3

output: 
- the number of (almost) 2-elementary Df^2 with |D| <= X
- the biggest absolute value of these discriminants
- the biggest class number of these discriminants
----------------------------------------------------------------------*/

twoeldisx(X)=twoeldisv(fundamx(X))

almdisx(X)=almdisv(fundamx(X))



/*-----------------------------------------------------
counting almost 2-elementary  discriminants  with fundamental discriminant having the given number of prime divisors

input: m a positive integer 
X>=3

output:
- the number of almost 2-elementary  discriminants Df^2 with omega(D)=m and |D|<=X
- the biggest absolute value of these discriminants
- the biggest class number of these discriminants
- same for the almost 2-elementary discriminants
-------------------------------------------------------------------------*/

almdismx(m,X) =

{
my(v,w);

v=sfome(3,X,m);
v=-onlymod(v,4,3);
w=almdisv(v);

v=sfome(3,X/4,m-1);
v=-4*onlymod(v,4,1);
w[1]=w[1]+almdisv(v)[1];
w[2]=max(w[2],almdisv(v)[2]);
w[3]=max(w[3],almdisv(v)[3]);

v=sfome(3,X/8,m-1);
v=-8*v;
w[1]=w[1]+almdisv(v)[1];
w[2]=max(w[2],almdisv(v)[2]);
w[3]=max(w[3],almdisv(v)[3]);

w


}

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C. Scripts for Section 5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

/*--------------------------------------------------------------------------------
kronecker, but with Delta defined as a residue class
input: Delta=Mod(r,m)
p a prime number, p|m

output kronecker(a/p)
--------------------------------------------------------------------------------*/


kronmod(Delta,p)=
{
my(m,r);

if (isprime(p)==0,
print(p" not prime");
return();
);

m=Delta.mod;
r=lift(Delta);

if (m%p!=0,
print ("bad prime");
return();
);

kronecker(r,p)
}

/*-----------------------------------------------
Psi as in (2.1), but with Delta defined as a residue class
input: Delta=Mod(r,m)
e an integer such that p|e => p|m, 2|e => 4|m, r=0,1 mod4

output: Psi(e,r)
---------------------------------------------------*/


relcln(e,Delta)=
{
my(h,m,r);

m=Delta.mod;
r=lift(Delta);
h=e;

fordiv(e,p,

if(isprime(p)==0,
next; 
);

if (p==2,

if(m%4!=0,
print(Delta "bad Delta");
return();
);

if (r%4==2 || r%4==3,
print(Delta "bad Delta");
return();
);

);

if(m%p!=0,
print(Delta "bad Delta");
return();
);

/*
print("h="h", p="p);
print("kronecker="kronecker(r,p));
*/

h=h*(1-kronecker(r,p)/p);

);

h

}


/*----------------------------------------------------------
Q(n), see Section 2.4
--------------------------------------------------------------*/

calq(n)=

{

my(v);

v=divisors(n);

for(i=1,length(v),
v[i]= v[i]^2/n;
);

v

}


/*-------------------------------------------------------
possible denominators
input: n,ax,ex,ey positive integers
output: possible denominators ay assuming that x,y are n-isogenous, see Corollary 2.18
---------------------------------------------------*/

denoms(n,ax,ex,ey)= 

{

my(den,dens,v,counter);

counter=0;

v=calq(n);

dens=vector(length(v),i,0);

num=ax*ey/ex;

for(i=1,length(v),

den=num*v[i];

if (den==floor(den),

counter++;

dens[counter]=den;

);

);

dens=dens[1..counter];

return(dens);

}


/*---------------------------------------------------------
removing 1 on the left
input: v a vector
output: v, if v[1]!=1
left translation of v, if v[1]=1
------------------------------------------------------------*/

noone(v)=

{

if(v[1]==1, 

v=v[2..length(v)];

);

v

}



/*------------------------------------------------------------
composing Table 2
--------------------------------------------------------*/


exeyez()=

{
my(ex,ey,ez,Exy,Exz,Eyz,E,Delta,Deltas,S,d,dx,dy,dz,list,thecounter,rem);

S=[1,2,3,4,6,12];
list=vector(100,i,[1,1,1,1,1,1,Mod(0,1),""]);
thecounter=0;

for(k=1,6,
for(j=1,k,
for(i=1,j,

ex=S[i];
ey=S[j];
ez=S[k];

if(gcd(gcd(ex,ey),ez)>1,
next;
);

Exy=lcm(ex,ey);
Exz=lcm(ex,ez);
Eyz=lcm(ey,ez);

E=lcm(Exy,ez);

if(E==1,
Deltas=[Mod(0,1)];
);

if(E==2,
Deltas=[Mod(1,8), Mod(0,4)];
);

if(E==3,
Deltas=[Mod(1,3)];
);

if(E==4,
Deltas=[Mod(1,8)];
);

if(E==6,
Deltas=[Mod(1,24),Mod(4,12)];
);

if(E==12,
Deltas=[Mod(1,24)];
);

for (n=1,length(Deltas),
Delta=Deltas[n];
rem="";

d=relcln(E,Delta);
dx=relcln(ex,Delta);
dy=relcln(ey,Delta);
dz=relcln(ez,Delta);
dxy=relcln(Exy,Delta);
dxz=relcln(Exz,Delta);
dyz=relcln(Eyz,Delta);

/*

print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");

print("ex="ex);
print("ey="ey);
print("ez="ez);
print("Exy="Exy);
print("Exz="Exz);
print("Eyz="Eyz);
print("E="E);
print("Deltas="Deltas);
print("Delta="Delta);
print("d="d);
print("dxy="dxy);
print("dxz="dxz);
print("dyz="dyz);
print("dx="dx);
print("dy="dy);
print("dz="dz);

*/

if(dxy<d,
/* print("K(x,y) not L"); */
next;
);

if(dxz<d,
/* print("K(x,z) not L"); */
next;
);

if(dyz<d,
/* print("K(y,z) not L"); */
next;
);

if(2*dx<d,
/* print("K(x) too small"); */
next;
);

if(2*dy<d,
/* print("K(y) too small"); */
next;
);

if(2*dz<d,
/* print("K(z) too small"); */
next;
);



if(dx<d,
if(ey!=ez,
print("ex="ex", ey="ey", ez="ez" is impossible"); 
next;
);
rem="n=r"; 
);

if(dy<d,
if(ex!=ez,
print("ex="ex", ey="ey", ez="ez" is impossible");
next;
);
rem="m=r";
);

if(dz<d,
if(ex!=ey,
print("ex="ex", ey="ey", ez="ez" is impossible");
next;
);
rem="m=n";
);





thecounter++;
/* print("thecounter="thecounter);*/

list[thecounter]=[ex,d/dx,ey,d/dy,ez,d/dz,Delta,rem];
/* print(list[thecounter]); */

);
);
);
);

list=list[1..thecounter];
return(list);
}



/*---------------------------------------------------------------
composing Table 3 and solving the systems
------------------------------------------------------------*/

systems(ex,ey,axone,axtwo) = 

{

my(Aone,Bone,Atwo,Btwo,D,az,azone,aztwo,n,numsys);

n=ex*ey;

az=noone(denoms(n,1,ex,ey));

azone=denoms(n,axone,ex,ey);

aztwo=denoms(n,axtwo,ex,ey);

numsys=0;

Aone= 1-1/axone;

Atwo= 1-1/axtwo;

for(i=1,length(az),

for(jy=1,length(azone),
for(jz=jy,length(azone),

for(ky=1,length(aztwo),
for(kz=ky,length(aztwo),

numsys++;

Bone= 1/az[i]-1/azone[jy]-1/azone[jz]; 

Btwo= 1/az[i]-1/aztwo[ky]-1/aztwo[kz]; 

D=Aone*Btwo-Bone*Atwo;


/* debug
print(Aone);
print(Bone);
print(Atwo);
print(Btwo);
*/ 

if (D==0, 

print("Bad Data");

print ("a(z)="az[i]);

print ("a(y^\sigma1)="azone[jy]);
print ("a(z^\sigma1)="azone[jz]);

print ("a(y^\sigma2)="aztwo[ky]);
print ("a(y^\sigma2)="aztwo[ky]);

return;

);

);););););

print(numsys" systems solved");

print("no non-trivial solution detected");

return;

}


/*--------------------------------------------------------
solving systems from Table 4 az=3 or 4
---------------------------------------------------------*/

dumbsystems(az) =

{
my(as,ass,matsys,lens,lengths,ind,numsys,v);

as=[1,1,az];

if (az==4,

ass=
[[2],[8],[2,8];
[4],[16],[1];
[2,8],[8,32],[2]];
);

if (az==3,

ass=
[[2],[8],[24];
[3],[3],[1];
[6,24],[24],[8]];
);

lengths=matrix(3,3,i,j,length(ass[i,j])); 

lens=mtov(lengths);

v=vector(9,i,1);
numsys=0;

for (i=1,9,

ind=vtom(v,3,3);

/* debugging

print(ind);

*/

matsys  = matrix (3,3,i,j, 
1/as[i] - 1/ass[i,j][ind[i,j]]); 

/* debugging

print(matsys);

*/

if (matdet(matsys)==0, 
print("anomaly");
print(ind); 
return;
);

/* debugging
print(i);

print(v);

*/

numsys++;

if (v==lexaddone(v,lens), 

print(numsys" systems solved");

print("no non-trivial solution detected");

return;
);

v=lexaddone(v,lens);

);

}



/*--------------------------------------------------------
solving systems from Table 5 ex=1 or 2
---------------------------------------------------------*/

dumbsystemsbis(ex) =

{
my(as,ass,matsys,lens,lengths,ind,numsys,v);

as=[1,2,2];

if (ex==1,

ass=
[[2],[1],[4];
[2],[4],[1];
[4,16],[8],[2,8,32]];
);

if (ex==2,

ass=
[[8],[1],[4];
[8],[4],[1];
[16,64],[8],[2,8,32]];
);

lengths=matrix(3,3,i,j,length(ass[i,j])); 

lens=mtov(lengths);

v=vector(9,i,1);
numsys=0;

for (i=1,9,

ind=vtom(v,3,3);

/* debugging

print(ind);

*/

matsys  = matrix (3,3,i,j, 
1/as[i] - 1/ass[i,j][ind[i,j]]); 

/* debugging

print(matsys);

*/

if (matdet(matsys)==0, 
print("anomaly");
print(ind); 
return;
);

/* debugging
print(i);

print(v);

*/

numsys++;

if (v==lexaddone(v,lens), 

print(numsys" systems solved");

print("no non-trivial solution detected");

return;
);

v=lexaddone(v,lens);

);

}





/*---------------------------------------------------------------
systems for case (5.49)
------------------------------------------------------------*/

systemsbis(ex) = 

{

my(Aone,Bone,Atwo,Btwo,D,aztwo,n,numsys);

n=ex*6;

print(n);



aztwo=denoms(n,9,ex,3);



numsys=0;

Aone= 1-1/3;

Atwo= 1-1/9;





for(ky=1,length(aztwo),
for(kz=ky,length(aztwo),

numsys++;

Bone= 1/2-1/54-1/54; 

Btwo= 1/2-1/aztwo[ky]-1/aztwo[kz]; 

D=Aone*Btwo-Bone*Atwo;




if (D==0, 

print("Bad Data");

return;

);

););

print(numsys" systems solved");

print("no non-trivial solution detected");

return;

}







