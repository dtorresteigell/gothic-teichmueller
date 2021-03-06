def Tgroup(place_precision=150):
	"""
	Computes the generators of the triangle group (4,6,12)
	Input: 
	*
	Output:
	* [] : array with the generators
	"""
	K.<r3>=QuadraticField(3);
	L.<r2>=K.extension(x^2-2);
	P=PolynomialRing(L,'y');
	h=((P.0)^2 - 3 - r3);
	F.<a>=L.extension(h);
	pl = [ p for p in F.places(prec=place_precision) if p(r3) > 0 and p(a) > 0 and p(r2) > 0][0];
	A=(r2/2)*matrix([[1,1],[-1,1]]);
	B=matrix([[r3/2,r3+a+1/2],[-r3+a-1/2,r3/2]]);
	C=(r2/4)*matrix([[-r3-2*a-1,-3*r3-2*a-1],[3*r3-2*a+1,-r3+2*a-1]]);
	Vec=[[r2, sqrt(2)], [r3, sqrt(3)], [a, sqrt(sqrt(3) + 3)]];
	return [[K,L,F,r2,r3,a,pl],[A,B,C]]


def symb_to_complex(x, Vec=[]):
	"""
	Redefines x from a symbolic number to an actual complex number
	Input: 
	* x : symbolic number
	* Vec(opt) : array with the equivalence between defined numbers and elements in C
	Output:
	* X : corresponding complex number
	"""
	Sx=str(x);
	loc = '{ ';
	for v in Vec:
		Sx=Sx.replace(str(v[0]),'('+str(v[1])+')');
		loc=loc+'\''+str(v[0])+'\':'+str(v[0])+',';
	loc=loc[:-1];
	loc=loc+ '}';
	X = sage_eval(Sx,locals=eval(loc));
	return X


def symb_to_complex_matrix(M, Vec=[]):
	"""
	Redefines M from a symbolic matrix to an actual matrix in M(2,C)
	Input: 
	* M : symbolic matrix
	* Vec(opt) : array with the equivalence between defined numbers and elements in C
	Output:
	* A : corresponding matric in M(2,C)
	"""
	A=[]
	for i in range(M.dimensions()[0]):
		for j in range(M.dimensions()[1]):
			A+=[symb_to_complex(M[i,j], Vec)];
	A = matrix(M.dimensions()[0],M.dimensions()[1],A);
	return A
	

def Moebius_action(M,z,Vec=[]):
	"""
	Computes the image of the complex number z under the action of the Moebius transformation associated to M
	Input: 
	* M : matrix in SL(2,R)
	* z : complex number
	* Vec(opt) : array with the equivalence between defined numbers and elements in R
	Output:
	* Mx : image of z under M
	"""
	A=symb_to_complex_matrix(M, Vec);
	Mx=(A[0,0]*z+A[0,1])/(A[1,0]*z+A[1,1]);
	Mx=Mx.simplify_full();
	return Mx

def Moebius_action_n(M,z):
	"""
	Computes the image of the complex number z under the action of the Moebius transformation associated to M
	Input: 
	* M : matrix in SL(2,R)
	* z : complex number
	* Vec(opt) : array with the equivalence between defined numbers and elements in R
	Output:
	* Mx : image of z under M
	"""
	A=matrix(CC,M);
	Mx=(A[0,0]*z+A[0,1])/(A[1,0]*z+A[1,1]);
	return Mx

def fixed_point(M, Vec=[]):
	"""
	Computes the fixed points in H of the Moebius transformation associated to M
	Input: 
	* M : matrix in SL(2,R)
	* Vec(opt) : array with the equivalence between defined numbers and elements in R
	Output:
	* F : array with the fixed point(s) of M
	"""
	F=[];
	X=var('X');
	A=symb_to_complex_matrix(M, Vec);
	S=solve([(A[0,0]*X+A[0,1])-X*(A[1,0]*X+A[1,1])],[X]);
	for s in S:
		if  (imag(s.rhs())).n() >= 0:
			F+=[(s.rhs()).simplify_full()];
	return F


def draw_group(Vec,Title='the group',prec=10):
	G=plot([]);
	Pts=[]
	G+=plot(points([[real(p.n()),imag(p.n())] for p in Vec],size=50,color='blue',aspect_ratio=1));
	Pts+=[ [real(p.n()),imag(p.n())] for p in Vec];
	for i in range(len(Vec)):
		G+=plot(hyperbolic_arc( ( floor(10^prec*real(Vec[i].n()))*10^(-prec), floor(10^prec*imag(Vec[i].n()))*10^(-prec)), (floor(10^prec*real(Vec[(i+1) % len(Vec)].n()))*10^(-prec), floor(10^prec*imag(Vec[(i+1) % len(Vec)].n()))*10^(-prec)), color='red'));
#	G+=[ hyperbolic_arc( (real(Vec[i].n()),imag(Vec[i].n())), (real(Vec[(i+1) % len(Vec)].n()),imag(Vec[(i+1) % len(Vec)].n())), color='red') for i in range(len(Vec))];
	show(plot(G),aspect_ratio=1, title='Fundamental domain of '+Title);
	return [G,Pts]


def galois(x,Vec):
	"""
	Computes the Galois conjugate of the element x by the action described in Vec
	Input: 
	* x : element in the field
	* Vec : array with pairs [a,sigma(a)]
	Output:
	* sigma_x : Galois conjugate element of x
	"""
	Sx=str(x);
	loc = '{';
	for v in Vec:
		Sx=Sx.replace(str(v[0]),'('+str(v[1])+')');
		loc=loc+'\''+str(v[0])+'\':'+str(v[0])+',';
	loc=loc[:-1];
	loc=loc+ '}';
	sigma_x = sage_eval(Sx,locals=eval(loc));
	return sigma_x


def galois_matrix(M,Vec):
	"""
	Computes the Galois conjugate of the matrix M by the action described in Vec
	Input: 
	* M : matrix with entries in the field
	* Vec : array with pairs [a,sigma(a)]
	Output:
	* sigma_M : Galois conjugate element of M
	"""
	A=[];
	for i in range(M.dimensions()[0]):
		for j in range(M.dimensions()[1]):
			A+=[galois(M[i,j],Vec)];
	sigma_M = matrix(M.dimensions()[0],M.dimensions()[1],A);
	return sigma_M


def rotate(V,n):
	"""
	Rotates V by n positions
	Input: 
	* V : vector
	* n : integer
	Output:
	* aux : vector rotated from V by n positions
	"""
	aux=vector([V[(j-n) % len(V)] for j in range(len(V))]);
	return aux


def matrix_to_pari(Mat, str="M"):
	"""
	Prints the code to define the matrix Mat in Pari
	Input: 
	* Mat : matrix
	* str (opt): name of the matrix in Pari
	Output:
	"""
	print("\n", str, " = matrix(",Mat.dimensions()[0],",",Mat.dimensions()[1],");");
	for i in range(Mat.dimensions()[0]):
		for j in range(Mat.dimensions()[1]):
			if Mat[i][j] != 0:
				print(str,"[",i+1,",",j+1,"]=",Mat[i][j]);
	return
	

def simp_mat(M):
	"""
	Simplifies the matrix M
	Input: 
	* M : matrix
	Output:
	* M : simplified matrix
	"""
	M=matrix([[v.simplify_full() for v in w] for w in M.rows()])
	return M


def rat(V, r3, R):
	"""
	Gives the rational representation of the element in the quaternion algebra (-1,3+√3) over Q(√3) corresponding to the vector V \n
	Input:  \n
	* V : vector [ [Vn[0],Vn[1]] for n in (0..3)]; the entries of the vectors [a,b] in V give the coefficients a+b*sqrt(3) of 1, i, j and ij respectively \n
	* r3 : rational representation of √3 \n
	* R : vector with the rational representations of 1, i, j and ij respectively \n
	Output: \n
	* rat : rational representation of (V0[0]+V0[1]*√3) + (V1[0]+V1[1]*√3)i + (V2[0]+V2[1]*√3)j + (V3[0]+V3[1]*√3)ij \n
	"""
	Id=identity_matrix(8,8);
	rat= (V[0][0]*Id + V[0][1]*r3)*R[0] + (V[1][0]*Id + V[1][1]*r3)*R[1] + (V[2][0]*Id + V[2][1]*r3)*R[2] + (V[3][0]*Id + V[3][1]*r3)*R[3];
	return rat

def HGcoefficients(a,b,c,Len=100):
	"""
	Gives a list of length Len with the coefficients of the hypergeometric series F(a,b,c;z), defined by
		( (a)_n * (b)_n ) / ( (c)_n * n! ), where (x)_n is the Pocchamber symbol.
	Input: 
	* a, b, c : parameters of the series
	* Len (opt): length of the series
	Output:
	* Coef : list with the coefficients
	"""
	Coef = [1];
	for i in range(Len):
		Coef+=[Coef[i]*(a+i)*(b+i)/((c+i)*(i+1))];
	return Coef

def HtoD(z):
	"""
	Isomorphism between the upper half-plane H and the unit disc D
	Matrix: HD = 1/(1+I)*matrix([[1,-I],[1,I]]) CUIDAO
	Input: 
	* z : point in H
	Output:
	* fz : image in D
	"""
	fz=(z-I)/(z+I);
	fz=fz.simplify_full().canonicalize_radical();
	return fz

def DtoH(z):
	"""
	Isomorphism between the unit disc D and the upper half-plane H
	Matrix: DH = (1/sqrt(2))*(-I)*matrix([[1,1],[1,-1]]) CUIDAO
	Input: 
	* z : point in D
	Output:
	* gz : image in H
	"""
	gz=-I*(z + 1)/(z - 1)
	gz=gz.simplify_full().canonicalize_radical();
	return gz
