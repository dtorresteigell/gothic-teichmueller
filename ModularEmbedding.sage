load("./kobayashi.sage")
R.<t> = PowerSeriesRing(QQ)

Len=500;


print("\n");
print(" ########################################################################### ");
print(" #                                                                         # ");
print(" #                COMPUTING MODULAR EMBEDDING OF Δ(3,6,∞)                  # ");
print(" #                                                                         # ");
print(" ########################################################################### ");



###########################################
######## TRACE FIELD CALCULATIONS #########
###########################################
K2.<r2>=NumberField(x^2-2);
K.<r3>=K2.extension(x^2-3);
place_precision=50;
pl = [ p for p in K.places(prec=place_precision) if (p(r2)>0 and p(r3)>0)][0]
polK.<x> = PolynomialRing(K);

F=coerce(polK,(x^2-2)^2-3);
aux=[f[0] for f in F.roots() if abs(pl(f[0])-sqrt(sqrt(3)+2))<0.001][0]

V1sym=r3-1;
V1=sqrt(3)-1;

H2D=matrix([[exp(-pi*I/12)*exp(pi*I/4),0],[0,exp(pi*I/12)*exp(-pi*I/4)]])*matrix([[1,-I],[1,I]])/(1+I); ## sends \infty to exp(pi*I/3)
V1H=Moebius_action(H2D^-1,V1);

# Rotation of pi/3 around 0
M3=matrix([[exp(pi*I/3),0],[0,exp(-pi*I/3)]]);
M3t=matrix([[exp(5*pi*I/3),0],[0,exp(-5*pi*I/3)]]);   ## sigma(zeta_6) = zeta_6^5
M3=(M3/sqrt(det(M3))).simplify_full()
M3H= (H2D^-1)*M3*H2D;
M3Hsymb = matrix( [[1/2,r3/2],[-r3/2,1/2]]);

# Rotation of pi/6 around 0
R6=matrix([[exp(pi*I/6),0],[0,exp(-pi*I/6)]]);

# Rotation of pi/6 around V1
A=matrix([[1,-V1],[-V1,1]]);
A1=matrix([[1,V1],[V1,1]]);
M6=(A1*R6*A).simplify_full();
M6=(M6/sqrt(det(M6))).simplify_full()
## (Moebius_action(M6,V1)-V1).simplify_full()
M6=matrix([[sqrt(3)/2+I*(1/2+2*sqrt(3)/3), -I*(sqrt(3)/3+1)],[ I*(sqrt(3)/3+1) , sqrt(3)/2-I*(1/2+2*sqrt(3)/3)]]);
M6t=matrix([[-sqrt(3)/2+I*(1/2-2*sqrt(3)/3), -I*(-sqrt(3)/3+1)],[ I*(-sqrt(3)/3+1) , -sqrt(3)/2-I*(1/2-2*sqrt(3)/3)]]);


M6H= (H2D^-1)*M6*H2D;
## M6Hsymb= matrix(2,2,[  [r for r in coerce(polK,algdep(f,4)).roots() if abs(pl(r[0])-f)<0.01][0][0] for f in M6H.list()]);
M6Hsymb = matrix([[-1/2 , 5/6*r3 + 1], [ -1/2*r3 , r3 + 1/2]]);

## Galois conjugate ???
M6Htsymb = matrix( [[ -1/6*r3 - 1 , -2/3*r3 + 1/2],[2/3*r3 - 1/2 , -5/6*r3 + 1] ]);
M6Ht = matrix(2,2, [pl(f) for f in M6Htsymb.list()]);

V1rep = Moebius_action(H2D,fixed_point(M6H)[0]);
V1trep = Moebius_action(H2D,fixed_point(M6Ht)[0]);
g=coerce(polK, algdep(V1trep,2));

# V1tsymb=[f[0] for f in g.roots() if abs(pl(f[0])-V1trep)<0.001][0];
V1t = 1/2-sqrt(3)/2;
V1tH=Moebius_action(H2D^-1,V1t);

# parabolic element
Minf= (M3*M6).simplify_full()
    ## print(fixed_point(Minf)[0].n());


### Symbolic generators of Delta(\infty,3,6) in \bH
M0=matrix([[1,2+2*r3/3],[0,1]])
M3=matrix([[1/2,r3/2],[-r3/2,1/2]])
M6=matrix([[r3+1/2,5*r3/6+1],[-r3/2,-1/2]])

## Conjugating Delta(\infty,3,6) to the Hilbert modular group
L=matrix(2,2,[r3, 42*r3 + 42, -1/12*r3 + 1/12 ,-2*r3]);
Ls=matrix(2,2,[-r3, -42*r3 + 42, 1/12*r3 + 1/12 ,2*r3]);
RP= -9*(3 + r3); IP= 3*(1 + r3)
RPs= -9*(3 - r3); IPs= -3*(1 - r3)
A=matrix([[1/IP,-RP/IP],[0,1]]);
As=matrix([[1/IPs,-RPs/IPs],[0,1]]);
nA= matrix(2,2,[pl(t) for t in A.list()]);
nAs= matrix(2,2,[pl(t) for t in As.list()]);

RV6=-(1+1/r3); IV6=r3/3;
RV6s=-(1-1/r3); IV6s=r3/3;
A6=matrix([[1/IV6,-RV6/IV6],[0,1]]);
A6s=matrix([[1/IV6s,-RV6s/IV6s],[0,1]]);
nA6= matrix(2,2,[pl(t) for t in A6.list()]);
nA6s= matrix(2,2,[pl(t) for t in A6s.list()]);

H0=(A^-1*A6*M0*A6^-1*A)
H3=(A^-1*A6*M3*A6^-1*A)
H6=(A^-1*A6*M6*A6^-1*A)


x=var('x');

#################################################
######## HGDE CORRESPONDING TO Δ(3,6,∞) #########
#################################################

c = 2/3;  # 1/p = |1-c|
a = 1/4;  # 1/q = |c-a-b| , 1/r = |a-b|
b = 1/4;  #

print("\n Hypergeometric differential equation associated to Δ(3,6,∞) given by\n\n    L := L(",a,",",b,",",c,") = t*(1-t)*y'' + (",c," - (",a,"+",b,"+1)*t)*y' -",a,"*",b,"*y = 0");

########## 1st SOLUTION AROUND t=0 ##########

# Coefficients of HG series F(a,b,c|t)
U1= HGcoefficients(a,b,c,Len);
y1=R(U1,Len+1);
dy1=y1.derivative(1);
ddy1=y1.derivative(2);
Ly1 = t*(1-t)*ddy1+ (c - (a+b+1)*t)*dy1 -a*b*y1;
print("\n The first solution around 0 is the hypergeometric series y1(t) = F(",a,",",b,",",c,"; t )");

#############################################


########## 2nd SOLUTION AROUND t=0 ##########

# Coefficients of HG series F(1+a-c,1+b-c,2-c|t)
aa=1+a-c;
bb=1+b-c;
cc=2-c;
U2= HGcoefficients(aa,bb,cc,Len);
y2=R(U2,Len+1);
  # This series must be multiplied by t^(1/3) to get the actual solution

print(" The second solution around 0 is the hypergeometric series y2(t) = t^(1/3)*F(",aa,",",bb,",",cc,"; t )");

# First derivative of y2
dy2=y2/(3*t)+y2.derivative(1);

# Second derivative of y2
ddy2=-2*y2/(9*t^2) +y2.derivative(1)*2/(3*t)+y2.derivative(2);

  # These series must still be multiplied by t^(1/3) to get the actual solution


Ly2 = t*(1-t)*ddy2+ (c - (a+b+1)*t)*dy2 -a*b*y2;

######## DOUBLE CHECK ########
# Py2=x^(1/3)*sum(U2[i]*x^i for i in range(Len+1));
# dPy2=Py2.derivative(1);
# ddPy2=Py2.derivative(2);
# LPy2 = (x*(1-x)*ddPy2+ (c - (a+b+1)*x)*dPy2 -a*b*Py2).simplify_full();
##############################

#############################################

y1_inv=y1.inverse()
y2_inv=y2.inverse()

zt=y2*y1_inv;

print("\n The coordinate z around 0 is the quotient of the period maps f2(z)/f1(z), so it should correspond to y2(t(z))/y1(t(z))");

Psi=t*zt^3; ## Q in pdf
h=Psi.reverse();

print("\n The map Psi ......................");




print("\n");

###################################################
######## HGDE CORRESPONDING TO Δ(3,6,∞)^τ #########
###################################################

# IS THIS CORRECT????
c_t = 5/3;  #
a_t = 5/12;  #
b_t = 5/12;  #

print("\n Hypergeometric differential equation associated to Δ(3,6,∞)^τ given by \n\n    L^τ := L(",a_t,",",b_t,",",c_t,") = t*(1-t)*y'' + (",c_t," - (",a_t,"+",b_t,"+1)*t)*y' -",a_t,"*",b_t,"*y = 0");

########## 1st SOLUTION AROUND t=0 ##########

# Coefficients of HG series F(a_t,b_t,c_t|t)
U1t= HGcoefficients(a_t,b_t,c_t,Len);
y1_t=R(U1t,Len+1);
dy1_t=y1_t.derivative(1);
ddy1_t=y1_t.derivative(2);
Ly1_t = t*(1-t)*ddy1_t+ (c_t - (a_t+b_t+1)*t)*dy1_t -a_t*b_t*y1_t;
print("\n The first solution around 0 is the hypergeometric series y1_t(t) = F(",a_t,",",b_t,",",c_t,"; t )");

#############################################


########## 2nd SOLUTION AROUND t=0 ##########

# Coefficients of HG series F(1+a_t-c_t,1+b_t-c_t,2-c_t|t)
aa_t=1+a_t-c_t;
bb_t=1+b_t-c_t;
cc_t=2-c_t;
U2t= HGcoefficients(aa_t,bb_t,cc_t,Len);
y2_t=R(U2t,Len+1);

print(" The second solution around 0 is the hypergeometric series y2_t(t) = t^(-2/3)*F(",aa_t,",",bb_t,",",cc_t,"; t )");

# First derivative of y2_t
dy2_t=-2*y2_t/(3*t)+y2_t.derivative(1);

# Second derivative of y2_t
ddy2_t=10*y2_t/(9*t^2) -4*y2_t.derivative(1)/(3*t)+y2_t.derivative(2);

Ly2_t = t*(1-t)*ddy2_t+ (c_t - (a_t+b_t+1)*t)*dy2_t -a_t*b_t*y2_t;

######## DOUBLE CHECK ########
# Py2_t=x^(1/4)*sum(U2t[i]*x^i for i in range(Len+1));
# dPy2_t=Py2_t.derivative(1);
# ddPy2_t=Py2_t.derivative(2);
# LPy2_t = (x*(1-x)*ddPy2_t+ (c_t - (a_t+b_t+1)*x)*dPy2_t -a_t*b_t*Py2_t).simplify_full();
##############################

#############################################

y1_t_inv=y1_t.inverse()
y2_t_inv=y2_t.inverse()

zt_t=y2_t*y1_t_inv;

print("\n The coordinate z around 0 is the quotient of the period maps f2_t(z)/f1_t(z), so it should correspond to y2_t(t(z))/y1_t(t(z))");

Psi_t=t*zt_t^4;
h_t=Psi_t.reverse();

print("\n The map Psi_t ......................");




###################################################
################ MODULAR EMBEDDING ################
###################################################

Phi_bar=Psi_t.subs(t=h)

# Calling x the z-coordinate in the disc
x=var('x')
PolPhi=Phi_bar.polynomial();
PolPhi=PolPhi.subs(t=x^4)
PhiNew = PolPhi.power_series(QQ);
Phi=PhiNew.sqrt().sqrt()





###################################################

print("\n\n\n Just to sum things up: \n - a, b, c: parameters of the HGDE L corresponding to Δ(4,6,12); \n - y1(t), y2(t): basis of solutions of L near t=0; \n - h: map sending z^4 to t(z); \n - Psi: map sending t to z(t)^4, inverse of h; \n - a_t, b_t, c_t:- parameters of the HGDE L^τ corresponding to Δ(4,6,12)^τ; \n - y1_t(t), y2_t(t): basis of solutions of L^τ near t=0; \n - h_t: map sending z^4 to t(z); \n - Psi_t: map sending t to z(t)^4, inverse of h_t; \n ");