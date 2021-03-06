"""

def EisensteinCalc(D,filename="../pari/Test"):
	#filename="../pari/Test";
	E = EisensteinForms(Gamma1(D),2);
	E.set_precision(500);
	Ser=E.eisenstein_series();
	Text1="default(breakloop,0);\ndefault(format,\"g.10\");\nallocatemem(2048000000);\n\nmysigma (n) = {\\\\ return sigma_1 for >0, -1/24 for 0\n  if (n==0, return (-1/24));\n  return(sigma(n));\n}\n\nH(D) = { \\\\ calculate H(2,D)\n  local(HH,e);\n  if ( (D%4!=0) && (D%4!=1), print(\"Invalid Discriminant!\"); return());\n  HH=0;\n  forstep (e=D%2,sqrt(D),2,\n    if (e==0,\n      HH+=mysigma((D-e^2)/4),\n      HH+=2*mysigma((D-e^2)/4));\n  );\n  HH/=-5;\n  if (issquare(D),HH-=D/10);\n  return(HH);\n};\n\neulercharX(D) = { \\\\ calculate chi(X_D)\n  local (DD,f,z,chi);\n  z=coredisc(D,1);\n  DD=z[1]; f=z[2]; \\\\ split into fundamental discriminant and square\n  chi=sumdiv(f,r,kronecker(DD,r)*moebius(r)/r^2);\n  chi*=-H(DD)/12;\n  chi*=2*f^3;\n  return(chi);\n};\n\nfor(j=4,50,eval(Str(Str(Str(eval(Str(zeta,j)),\"=exp(2*Pi*I/\"),j),\")\")))\n";
	Text2="\nEM = matrix("+str(len(Ser))+",500,i,j,polcoeff(EV[i],j-1,q));\n\n\\r LevelsShimura\n\nL=LevelsShimuraTN(D);\nvv=L[1]; CV=L[2];\n\nCVP = CV + vector(length(CV),i,1); \\\\ shift by one because of indexing\nV=vecextract(EM,CVP);\nK=matker(V~);\n\n{squaredivisors(x)=local(k,aux);\nk=divisors(x); aux=[];\nfor(j=2,length(k),if(issquare(k[j]), aux=setunion(aux,[k[j]])));\nreturn(aux);\n};\n\n{conductor(D)=local(f,L,aux);\nL=squaredivisors(D);\naux=[];\nfor(j=1,length(L),if(isdiscriminant(D/L[j])==1, aux=concat(aux,L[j]);));\nif(aux==[],return(1), if(aux==[aux[1]],return(round(sqrt(aux[1]))), f=aux[1]; for(j=2,length(aux), f=max(f,aux[j]));return(round(sqrt(f)))));\n}\n\n{eulercharX16(D)=local(f);\nf=conductor(D);\nif(f%6==1,fact=1,if(f%6==2,fact=3/2,if(f%6==3,fact=4/3,if(f%6==6,fact=2))));\nreturn(eulercharX(D)*fact);\n};\n\nT=(round(K)~*EM); \nif(length(T[,1])>1,print("Not unique solution"), \n\\\\ Volumes of the T_N\nT= -T/T[1]*eulercharX16(D)/4;\n\n\nF=vector(499);\nfor(j=1,499, sq=squaredivisors(j); F[j]=T[j+1]-sum(l=1,length(sq),F[j/sq[l]]));   \\\\ Volumes of the F_N\nfor(j=1,499, if(F[j]!=0, print(\" · vol(F_\",j,\") = \",F[j])))\n\nM=500 -1;\nL=LevelsShimura(D,M);\nvv=L[1]; CV=L[2];\nfor(j=1,M,if( bitxor(T[j]==0 , setsearch(CV,j-1)!=0)!=0, print(j,\" WRONG!\"))))\n\n\n";
	out = file(filename,'w');
	out.write(Text1+'\n');
	out.write("\nD= "+str(D)+";\n\n");
	out.write("{EV = ");
	out.write(str(Ser));
	out.write("};");
	out.write(Text2+'\n');
	out.close();





def EisensteinCalcChar(D,filename="../pari/Test"):
	#filename="../pari/Test";
	E=EisensteinForms(DirichletGroup(D)[0],2)
	E.set_precision(500);
	Ser=E.eisenstein_series();
	Text1="default(breakloop,0);\ndefault(format,\"g.10\");\nallocatemem(2048000000);\n\nmysigma (n) = {\\\\ return sigma_1 for >0, -1/24 for 0\n  if (n==0, return (-1/24));\n  return(sigma(n));\n}\n\nH(D) = { \\\\ calculate H(2,D)\n  local(HH,e);\n  if ( (D%4!=0) && (D%4!=1), print(\"Invalid Discriminant!\"); return());\n  HH=0;\n  forstep (e=D%2,sqrt(D),2,\n    if (e==0,\n      HH+=mysigma((D-e^2)/4),\n      HH+=2*mysigma((D-e^2)/4));\n  );\n  HH/=-5;\n  if (issquare(D),HH-=D/10);\n  return(HH);\n};\n\neulercharX(D) = { \\\\ calculate chi(X_D)\n  local (DD,f,z,chi);\n  z=coredisc(D,1);\n  DD=z[1]; f=z[2]; \\\\ split into fundamental discriminant and square\n  chi=sumdiv(f,r,kronecker(DD,r)*moebius(r)/r^2);\n  chi*=-H(DD)/12;\n  chi*=2*f^3;\n  return(chi);\n};\n\nfor(j=4,50,eval(Str(Str(Str(eval(Str(zeta,j)),\"=exp(2*Pi*I/\"),j),\")\")))\n";
	Text2="\nEM = matrix("+str(len(Ser))+",500,i,j,polcoeff(EV[i],j-1,q));\n\n\\r LevelsShimura\n\nL=LevelsShimuraTN(D);\nvv=L[1]; CV=L[2];\n\nCVP = CV + vector(length(CV),i,1); \\\\ shift by one because of indexing\nV=vecextract(EM,CVP);\nK=matker(V~);\n\n{squaredivisors(x)=local(k,aux);\nk=divisors(x); aux=[];\nfor(j=2,length(k),if(issquare(k[j]), aux=setunion(aux,[k[j]])));\nreturn(aux);\n};\n\n{conductor(D)=local(f,L,aux);\nL=squaredivisors(D);\naux=[];\nfor(j=1,length(L),if(isdiscriminant(D/L[j])==1, aux=concat(aux,L[j]);));\nif(aux==[],return(1), if(aux==[aux[1]],return(round(sqrt(aux[1]))), f=aux[1]; for(j=2,length(aux), f=max(f,aux[j]));return(round(sqrt(f)))));\n}\n\n{eulercharX16(D)=local(f);\nf=conductor(D);\nif(f%6==1,fact=1,if(f%6==2,fact=3/2,if(f%6==3,fact=4/3,if(f%6==6,fact=2))));\nreturn(eulercharX(D)*fact);\n};\n\nT=(K~*EM)[1,]; \\\\T=round(T);\n\\\\ Volumes of the T_N\nT= -T/T[1]*eulercharX16(D)/4;\n\n\nF=vector(499);\nfor(j=1,499, sq=squaredivisors(j); F[j]=T[j+1]-sum(l=1,length(sq),F[j/sq[l]]));   \\\\ Volumes of the F_N\nfor(j=1,499, if(F[j]!=0, print(\" · vol(F_\",j,\") = \",F[j])))\n\nM=500 -1;\nL=LevelsShimura(D,M);\nvv=L[1]; CV=L[2];\nfor(j=1,M,if( bitxor(T[j]==0 , setsearch(CV,j-1)!=0)!=0, print(j,\" WRONG!\")))\n\n\n";
	out = file(filename,'w');
	out.write(Text1+'\n');
	out.write("\nD= "+str(D)+";\n\n");
	out.write("{EV = ");
	out.write(str(Ser));
	out.write("};");
	out.write(Text2+'\n');
	out.close();

"""


def Quad(D,M):
	return det(M)/(-6)/D;

def BB(D,X,Y):
	return Quad(D,X+Y) - Quad(D,X) - Quad(D,Y);

def norm6gen(D):
	K.<r> = QuadraticField(D);
	if mod(D,4)==1:
		gam=(1+r)/2;
		gamc=(1-r)/2;
	else:
		gam=r/2;
		gamc=-r/2;
	B=[];
	for a in [0..5]:
		n=(a+gam)*(a+gamc);
		if mod(n,6)==0 and n!=0:
			B.append(a);
	return B;

def getlattices(D):
	Lat=[];
	K.<r> = QuadraticField(D);
	if mod(D,4)==1:
		gam=(1+r)/2;
		gamc=(1-r)/2;
	else:
		gam=r/2;
		gamc=-r/2;
	B=norm6gen(D);
	for b in B:
		gen=b+gam;
		genc=b+gamc;
		gen = r*gen;
		genc = -r*genc;
		M1 = matrix(K,2,2,[1,0,0,0]);
		M2 = matrix(K,2,2,[0,0,0,6*33]);
		M3 = matrix(K,2,2,[0,r*6,-r*6,0]);
		M4 = matrix(K,2,2,[0,gen,genc,0]);
		V=list(vector(ZZ,[0,1,0,0, 1,0,0,0, 0,0,-BB(D,M3,M3),-BB(D,M3,M4), 0,0,-BB(D,M3,M4),-BB(D,M4,M4)]));
		Lat.append(V);
	return Lat;

def extract(D,filename= "../pari/Test", name="V", P=10):
	load("./modforms/integer_lattice.py");
	out = file(filename,'w');
	Lat=getlattices(D);
	out.write('{'+name+'=[');
	B=norm6gen(D);
	for k in range(len(B)):
		lat=Lat[k];
		L=Lattice(matrix(4,4,lat));
		E=L.eisenstein_series(2,prec=P);
		V=[ [ [v[0][2],v[0][3]], [ [j[0],j[1]] for j in v[1].items()]  ] for v in E.items() ];
		out.write('['+str(B[k]) + ' , [');
		for v in V:
			out.write('[ '+str(v[0])+', [');
			for l in v[1]:
				out.write(str(l[0])+', '+str(l[1]) );
				if l!=v[1][len(v[1])-1]:
					out.write(' ; ');
				else:
					out.write(' ] ]');
			if v != V[len(V)-1]:
				out.write(',');
			else:
				out.write(']');
		out.write('  ] ');
		if lat!=Lat[len(B)-1]:
			out.write(',\n');
	out.write("]};");
	out.close();


def generateEisenstein(dd,DD,P=10):
	aux=[0,1,4,9,12,16];
	for D in [dd..DD]:
		if (mod(D,24) in aux) & (is_square(D)==False):
			filename="../pari/EisList"+str(D);
			name="EisList"+str(D);
			extract(D,filename,name,P);
			print("\n * Finished calculating Eisenstein series for D = ",D);


def TNvsGamma0(D,N):
	load("./modforms/integer_lattice.py");
	G=Gamma0(N);
	Lat=getlattices(D);
	B=norm6gen(D);
	for k in range(len(B)):
		lat=Lat[k];
		L=Lattice(matrix(4,4,lat));
		E=L.eisenstein_series(2,prec=P);
		V=[ [ [v[0][2],v[0][3]], [ [j[0],j[1]] for j in v[1].items()]  ] for v in E.items() ];
		for w in V:
			L=[v[0] for v in w[1]];
			if N/D in L:
				pos=L.index(N/D);
				vol=w[1][pos][1]*eulercharX16(D)/4;
				if -G.index()/12 == vol:
					print("  The volumes of Gamma0(",N,") and T_",N,"(nu) agree.");
					return(1);
				else:
					print("  The volumes of Gamma0(",N,") and T_",N,"(nu) DO NOT agree!!!!!!!!!");
					print("       · Vol(T_",N,"(nu)) = ",-vol);
					print("       · Vol(Gamma0(",N,")) = ",G.index()/12);
					return(0);


def FNvsGamma0(D,N,r,P=10):
	"""
	Checks if F_N(\nu) and Gamma_O(N) have the same volume, for the ideal corresponding to k in norm6gen(D)
	"""
	load("./modforms/integer_lattice.py");
	G=Gamma0(N);
	Lat=getlattices(D);
	B=norm6gen(D);
	k=B.index(r);
	lat=Lat[k];
	L=Lattice(matrix(4,4,lat));
	E=L.eisenstein_series(2,prec=P);
	V=[ [ [v[0][2],v[0][3]], [ [j[0],j[1]] for j in v[1].items()]  ] for v in E.items() ];
	for w in V:
		L=[v[0] for v in w[1]];
		if N/D in L:
			pos=L.index(N/D);
			vol=w[1][pos][1]*eulercharX16(D)/4;
			[n for n in N.divisors() if n.is_square() and n!=1];
			for n in [n for n in N.divisors() if n.is_square() and n!=1]:
				check=0; count=0;
				while check==0 and count<=len(V)-1:
					w=V[count];
					L=[v[0] for v in w[1]];
					if N/n/D in L:
						check=1;
						pos=L.index(N/n/D);
						vol-=w[1][pos][1]*eulercharX16(D)/4;
					count+=1;
			if -G.index()/12 == vol:
				print("  The volumes of Gamma0(",N,") and F_",N,"(nu) agree.");
				return(1);
			else:
				print("  The volumes of Gamma0(",N,") and F_",N,"(nu) DO NOT agree!!!!!!!!!");
				print("       · Vol(F_",N,"(nu)) = ",-vol);
				print("       · Vol(Gamma0(",N,")) = ",G.index()/12);
				return(0);


def eulercharX16(D):
	aux=[[12, 1/3], [24, 1], [28, 4/3], [33, 2], [40, 7/3], [48, 4], [52, 5], [57, 14/3], [60, 4], [72, 20/3], [73, 22/3], [76, 19/3], [84, 10], [88, 23/3], [96, 12], [97, 34/3], [105, 12], [108, 12], [112, 16], [120, 34/3], [124, 40/3], [129, 50/3], [132, 18], [136, 46/3], [145, 64/3], [148, 25], [153, 80/3], [156, 52/3], [160, 28], [168, 18], [172, 21], [177, 26], [180, 40], [184, 74/3], [192, 32], [193, 98/3], [201, 98/3], [204, 26], [208, 40], [216, 36], [217, 116/3], [220, 92/3], [228, 42], [232, 33], [240, 48], [241, 142/3], [244, 55], [249, 46], [252, 128/3], [264, 112/3], [265, 160/3], [268, 41], [273, 148/3], [276, 60], [280, 134/3], [288, 80], [292, 66], [297, 72], [300, 130/3], [304, 76], [312, 46], [313, 200/3], [316, 56], [321, 66], [328, 54], [336, 80], [337, 76], [340, 90], [345, 220/3], [348, 52], [352, 92], [360, 224/3], [364, 206/3], [369, 320/3], [372, 90], [376, 212/3], [384, 96], [385, 284/3], [388, 102], [393, 86], [396, 280/3], [408, 206/3], [409, 316/3], [412, 76], [417, 94], [420, 108], [424, 87], [432, 144], [433, 326/3], [436, 135], [444, 244/3], [448, 128], [456, 88], [457, 120], [460, 278/3], [465, 112], [468, 160], [472, 277/3], [480, 136], [481, 404/3], [489, 374/3], [492, 90], [496, 160], [504, 400/3], [505, 144], [508, 320/3], [513, 168], [516, 150], [520, 346/3], [528, 144], [532, 170], [537, 134], [540, 144], [544, 184], [552, 308/3], [553, 476/3], [556, 127], [561, 472/3], [564, 180], [568, 126], [577, 500/3], [580, 192], [585, 640/3], [588, 350/3], [592, 200], [600, 120], [601, 566/3], [604, 148], [609, 176], [612, 240], [616, 452/3], [624, 208], [628, 215], [633, 170], [636, 140], [640, 224], [648, 180], [649, 638/3], [652, 467/3], [657, 704/3], [660, 220], [664, 503/3], [672, 216], [673, 646/3], [681, 210], [684, 608/3], [688, 252], [696, 478/3], [697, 228], [700, 520/3], [705, 628/3], [708, 234], [712, 512/3], [720, 320], [721, 752/3], [724, 285], [732, 164], [736, 296], [744, 182], [745, 252], [748, 550/3], [753, 230], [756, 360], [760, 610/3], [768, 256], [769, 818/3], [772, 294], [777, 704/3], [780, 180], [792, 736/3], [793, 832/3], [796, 220], [801, 1040/3], [804, 294], [808, 207], [816, 312], [817, 842/3], [820, 340], [825, 260], [828, 800/3], [832, 320], [840, 620/3], [844, 737/3], [849, 294], [852, 300], [856, 757/3], [864, 432], [865, 324], [868, 348], [873, 1088/3], [876, 692/3], [880, 368], [888, 210], [889, 1040/3], [892, 240], [897, 296], [904, 812/3], [912, 336], [913, 994/3], [916, 405], [921, 330], [924, 248], [928, 396], [936, 1000/3], [937, 1030/3], [940, 282], [945, 432], [948, 350], [952, 268], [960, 384], [964, 426], [969, 348], [972, 324], [976, 440], [984, 842/3], [985, 384], [988, 292], [993, 334], [996, 414], [1000, 875/3]];
	d=[v[0] for v in aux];
	if D in d:
		return aux[d.index(D)][1];
	else:
		print("ERROR Discriminant not found. Make sure that D=0,1,4,9,12,16 mod 24 and D in [12..1000]");
	return 0;


def checkTNGamma0(D,NN=50):
	for N in [1..NN]:
		if is_squarefree(N):
			TNvsGamma0(D,N);



Levels= [[12, []], [24, [1]], [28, [1]], [33, [1]], [40, [1]], [48, [2]], [52, [2]], [57, [2]], [60, [1]], [72, [3]], [73, [3, 2, 1]], [76, [3]], [84, [2]], [88, [3, 1]], [96, [4]], [97, [4, 1, 3, 2]], [105, [4, 1, 1]], [108, [3]], [112, [4, 2]], [120, [5]], [124, [5, 1]], [129, [5, 2]], [132, [4]], [136, [5, 3]], [145, [6, 5, 4, 1, 1]], [148, [6, 2]], [153, [6, 3]], [156, [5]], [160, [6, 4]], [168, [7, 1]], [172, [7, 3]], [177, [7, 4, 1]], [180, [6]], [184, [7, 5]], [192, [8, 2]], [193, [8, 2, 7, 6, 3, 1]], [201, [8, 2, 5]], [204, [7]], [208, [8, 6]], [216, [9, 3]], [217, [9, 1, 8, 2, 7, 4, 1, 2]], [220, [9, 1, 5, 1]], [228, [8]], [232, [9, 1, 7]], [240, [10, 4]], [241, [10, 9, 1, 8, 2, 5, 3]], [244, [10, 6, 2]], [249, [10, 7, 1]], [252, [9]], [264, [11, 5]], [265, [11, 10, 9, 1, 6, 4, 1]], [268, [11, 7, 3]], [273, [11, 8, 2, 2]], [276, [10]], [280, [11, 9, 1, 1]], [288, [12, 6]], [292, [12, 8, 4]], [297, [12, 3, 9, 3]], [300, [11]], [304, [12, 10, 2]], [312, [13, 7]], [313, [13, 12, 3, 11, 8, 2, 6, 1]], [316, [13, 9, 1, 5]], [321, [13, 10, 4, 1]], [328, [13, 11, 3]], [336, [14, 8]], [337, [14, 13, 12, 3, 9, 1, 7, 2]], [340, [14, 10, 6]], [345, [14, 11, 5]], [348, [13, 1]], [352, [14, 12, 4]], [360, [15, 9]], [364, [15, 11, 7]], [369, [15, 12, 3, 6]], [372, [14, 2]], [376, [15, 13, 5]], [384, [16, 10]], [385, [16, 4, 1, 15, 14, 11, 9, 1, 4, 1, 1]], [388, [16, 12, 8]], [393, [16, 4, 1, 13, 7]], [396, [15, 3]], [408, [17, 11]], [409, [17, 16, 4, 1, 15, 12, 3, 10, 5, 2]], [412, [17, 13, 9, 1]], [417, [17, 14, 8, 2]], [420, [16, 4]], [424, [17, 15, 7, 1]], [432, [18, 12]], [433, [18, 2, 17, 16, 4, 1, 13, 11, 6, 3]], [436, [18, 2, 14, 10]], [444, [17, 5]], [448, [18, 2, 16, 8, 2]], [456, [19, 13]], [457, [19, 18, 2, 17, 14, 12, 3, 7, 4, 1]], [460, [19, 15, 11]], [465, [19, 16, 4, 1, 10, 1]], [468, [18, 6]], [472, [19, 17, 9, 1, 3]], [480, [20, 14]], [481, [20, 5, 19, 18, 2, 15, 13, 8, 2, 5]], [489, [20, 5, 17, 11, 2]], [492, [19, 7]], [496, [20, 18, 2, 10, 4]], [504, [21, 15]], [505, [21, 20, 5, 19, 16, 4, 1, 14, 9, 1, 6]], [508, [21, 17, 13, 1]], [513, [21, 18, 12, 3, 3]], [516, [20, 8]], [520, [21, 19, 11, 5]], [528, [22, 16]], [532, [22, 18, 2, 14, 2]], [537, [22, 19, 13, 4, 1]], [540, [21, 9]], [544, [22, 20, 12, 6]], [552, [23, 17]], [553, [23, 22, 21, 18, 2, 16, 4, 1, 11, 8, 2, 1]], [556, [23, 19, 15, 3]], [561, [23, 20, 5, 14, 5]], [564, [22, 10]], [568, [23, 21, 13, 7]], [577, [24, 6, 23, 22, 19, 17, 12, 3, 9, 1, 2]], [580, [24, 20, 16, 4]], [585, [24, 6, 21, 15, 6]], [588, [23, 11]], [592, [24, 22, 14, 8]], [600, [25, 19, 1]], [601, [25, 1, 24, 6, 23, 20, 5, 18, 2, 13, 10, 3]], [604, [25, 1, 21, 17, 5]], [609, [25, 1, 22, 16, 4, 1, 7]], [612, [24, 12]], [616, [25, 1, 23, 15, 9, 1]], [624, [26, 20, 2]], [628, [26, 22, 18, 2, 6]], [633, [26, 23, 17, 8, 2]], [636, [25, 1, 13]], [640, [26, 24, 16, 10]], [648, [27, 21, 3]], [649, [27, 3, 26, 25, 1, 22, 20, 5, 15, 12, 3, 5, 1]], [652, [27, 3, 23, 19, 7]], [657, [27, 24, 6, 18, 9]], [660, [26, 14]], [664, [27, 3, 25, 1, 17, 11]], [672, [28, 22, 4]], [673, [28, 7, 27, 3, 26, 23, 21, 16, 4, 1, 13, 6, 2]], [681, [28, 7, 25, 1, 19, 10]], [684, [27, 15]], [688, [28, 26, 18, 2, 12]], [696, [29, 23, 5]], [697, [29, 28, 7, 27, 3, 24, 6, 22, 17, 14, 7, 3]], [700, [29, 25, 21, 9, 1, 1]], [705, [29, 26, 20, 5, 11]], [708, [28, 16]], [712, [29, 27, 3, 19, 13]], [720, [30, 24, 6]], [721, [30, 29, 28, 7, 25, 1, 23, 18, 2, 15, 8, 2, 4, 1]], [724, [30, 26, 22, 10, 2]], [732, [29, 17]], [736, [30, 28, 20, 14]], [744, [31, 25, 1, 7]], [745, [31, 30, 29, 26, 24, 6, 19, 16, 4, 1, 9, 1, 5]], [748, [31, 27, 3, 23, 11, 3]], [753, [31, 28, 7, 22, 13, 1]], [756, [30, 18]], [760, [31, 29, 21, 15]], [768, [32, 26, 8]], [769, [32, 8, 2, 31, 30, 27, 3, 25, 1, 20, 5, 17, 10, 6]], [772, [32, 28, 24, 12, 4]], [777, [32, 8, 2, 29, 23, 14, 2]], [780, [31, 19]], [792, [33, 27, 9]], [793, [33, 32, 8, 2, 31, 28, 7, 26, 21, 18, 2, 11, 7]], [796, [33, 29, 25, 1, 13, 5]], [801, [33, 30, 24, 6, 15, 3]], [804, [32, 20]], [808, [33, 31, 23, 17, 1]], [816, [34, 28, 10]], [817, [34, 33, 32, 8, 2, 29, 27, 3, 22, 19, 12, 3, 8, 2]], [820, [34, 30, 26, 14, 6]], [825, [34, 31, 25, 16, 4, 1, 4, 1]], [828, [33, 21]], [832, [34, 32, 24, 18, 2, 2]], [840, [35, 29, 11]], [844, [35, 31, 27, 3, 15, 7]], [849, [35, 32, 8, 2, 26, 17, 5]], [852, [34, 22]], [856, [35, 33, 25, 1, 19, 3]], [864, [36, 30, 12]], [865, [36, 9, 4, 1, 35, 34, 31, 29, 24, 6, 21, 14, 10, 1]], [868, [36, 4, 32, 28, 16, 8]], [873, [36, 9, 33, 27, 18, 6]], [876, [35, 23]], [880, [36, 4, 34, 26, 20, 4]], [888, [37, 31, 13]], [889, [37, 36, 9, 4, 1, 35, 32, 8, 2, 30, 25, 1, 22, 15, 11, 2]], [892, [37, 33, 29, 17, 9, 1]], [897, [37, 34, 28, 7, 19, 7]], [904, [37, 35, 27, 3, 21, 5]], [912, [38, 32, 14]], [913, [38, 37, 36, 9, 4, 1, 33, 31, 26, 23, 16, 4, 1, 12, 3, 3]], [916, [38, 34, 30, 18, 2, 10]], [921, [38, 35, 29, 20, 5, 8, 2]], [924, [37, 25, 1, 1]], [928, [38, 36, 4, 28, 22, 6]], [936, [39, 33, 15]], [937, [39, 38, 37, 34, 32, 8, 2, 27, 3, 24, 6, 17, 13, 4, 1]], [940, [39, 35, 31, 19, 11]], [945, [39, 36, 9, 30, 21, 9]], [948, [38, 26, 2]], [952, [39, 37, 29, 23, 7]], [960, [40, 34, 16]], [964, [40, 36, 4, 32, 20, 12]], [969, [40, 10, 37, 31, 22, 10]], [972, [39, 27, 3]], [976, [40, 38, 30, 24, 8]], [984, [41, 35, 17]], [985, [41, 40, 10, 39, 36, 9, 4, 1, 34, 29, 26, 19, 15, 6, 1]], [988, [41, 37, 33, 21, 13]], [993, [41, 38, 32, 8, 2, 23, 11]], [996, [40, 28, 4]], [1000, [41, 39, 31, 25, 9, 1]]];

Levels=[[12, []], [24, [[0, 1]]], [28, [[1, 1]]], [33, [[4, 1]]], [40, [[2, 1]]], [48, [[0, 2]]], [52, [[1, 2]]], [57, [[4, 2]]], [60, [[3, 1]]], [72, [[0, 3]]], [73, [[3, 3], [5, 2], [0, 1]]], [76, [[1, 3]]], [84, [[3, 2]]], [88, [[2, 3], [4, 1]]], [96, [[0, 4]]], [97, [[3, 4], [3, 1], [5, 3], [0, 2]]], [105, [[4, 4], [4, 1], [1, 1]]], [108, [[3, 3]]], [112, [[2, 4], [4, 2]]], [120, [[0, 5]]], [124, [[1, 5], [5, 1]]], [129, [[4, 5], [1, 2]]], [132, [[3, 4]]], [136, [[2, 5], [4, 3]]], [145, [[3, 6], [5, 5], [0, 4], [0, 1], [2, 1]]], [148, [[1, 6], [5, 2]]], [153, [[4, 6], [1, 3]]], [156, [[3, 5]]], [160, [[2, 6], [4, 4]]], [168, [[0, 7], [0, 1]]], [172, [[1, 7], [5, 3]]], [177, [[4, 7], [1, 4], [1, 1]]], [180, [[3, 6]]], [184, [[2, 7], [4, 5]]], [192, [[0, 8], [0, 2]]], [193, [[3, 8], [3, 2], [5, 7], [0, 6], [2, 3], [3, 1]]], [201, [[4, 8], [4, 2], [1, 5]]], [204, [[3, 7]]], [208, [[2, 8], [4, 6]]], [216, [[0, 9], [0, 3]]], [217, [[3, 9], [3, 1], [5, 8], [5, 2], [0, 7], [2, 4], [2, 1], [3, 2]]], [220, [[1, 9], [1, 1], [5, 5], [1, 1]]], [228, [[3, 8]]], [232, [[2, 9], [2, 1], [4, 7]]], [240, [[0, 10], [0, 4]]], [241, [[3, 10], [5, 9], [5, 1], [0, 8], [0, 2], [2, 5], [3, 3]]], [244, [[1, 10], [5, 6], [1, 2]]], [249, [[4, 10], [1, 7], [4, 1]]], [252, [[3, 9]]], [264, [[0, 11], [0, 5]]], [265, [[3, 11], [5, 10], [0, 9], [0, 1], [2, 6], [3, 4], [3, 1]]], [268, [[1, 11], [5, 7], [1, 3]]], [273, [[4, 11], [1, 8], [1, 2], [4, 2]]], [276, [[3, 10]]], [280, [[2, 11], [4, 9], [4, 1], [2, 1]]], [288, [[0, 12], [0, 6]]], [292, [[1, 12], [5, 8], [1, 4]]], [297, [[4, 12], [4, 3], [1, 9], [4, 3]]], [300, [[3, 11]]], [304, [[2, 12], [4, 10], [2, 2]]], [312, [[0, 13], [0, 7]]], [313, [[3, 13], [5, 12], [5, 3], [0, 11], [2, 8], [2, 2], [3, 6], [5, 1]]], [316, [[1, 13], [5, 9], [5, 1], [1, 5]]], [321, [[4, 13], [1, 10], [4, 4], [4, 1]]], [328, [[2, 13], [4, 11], [2, 3]]], [336, [[0, 14], [0, 8]]], [337, [[3, 14], [5, 13], [0, 12], [0, 3], [2, 9], [2, 1], [3, 7], [5, 2]]], [340, [[1, 14], [5, 10], [1, 6]]], [345, [[4, 14], [1, 11], [4, 5]]], [348, [[3, 13], [3, 1]]], [352, [[2, 14], [4, 12], [2, 4]]], [360, [[0, 15], [0, 9]]], [364, [[1, 15], [5, 11], [1, 7]]], [369, [[4, 15], [1, 12], [1, 3], [4, 6]]], [372, [[3, 14], [3, 2]]], [376, [[2, 15], [4, 13], [2, 5]]], [384, [[0, 16], [0, 10]]], [385, [[3, 16], [3, 4], [3, 1], [5, 15], [0, 14], [2, 11], [3, 9], [3, 1], [5, 4], [5, 1], [0, 1]]], [388, [[1, 16], [5, 12], [1, 8]]], [393, [[4, 16], [4, 4], [4, 1], [1, 13], [4, 7]]], [396, [[3, 15], [3, 3]]], [408, [[0, 17], [0, 11]]], [409, [[3, 17], [5, 16], [5, 4], [5, 1], [0, 15], [2, 12], [2, 3], [3, 10], [5, 5], [0, 2]]], [412, [[1, 17], [5, 13], [1, 9], [1, 1]]], [417, [[4, 17], [1, 14], [4, 8], [4, 2]]], [420, [[3, 16], [3, 4]]], [424, [[2, 17], [4, 15], [2, 7], [4, 1]]], [432, [[0, 18], [0, 12]]], [433, [[3, 18], [3, 2], [5, 17], [0, 16], [0, 4], [0, 1], [2, 13], [3, 11], [5, 6], [0, 3]]], [436, [[1, 18], [1, 2], [5, 14], [1, 10]]], [444, [[3, 17], [3, 5]]], [448, [[2, 18], [2, 2], [4, 16], [2, 8], [4, 2]]], [456, [[0, 19], [0, 13]]], [457, [[3, 19], [5, 18], [5, 2], [0, 17], [2, 14], [3, 12], [3, 3], [5, 7], [0, 4], [0, 1]]], [460, [[1, 19], [5, 15], [1, 11]]], [465, [[4, 19], [1, 16], [1, 4], [1, 1], [4, 10], [1, 1]]], [468, [[3, 18], [3, 6]]], [472, [[2, 19], [4, 17], [2, 9], [2, 1], [4, 3]]], [480, [[0, 20], [0, 14]]], [481, [[3, 20], [3, 5], [5, 19], [0, 18], [0, 2], [2, 15], [3, 13], [5, 8], [5, 2], [0, 5]]], [489, [[4, 20], [4, 5], [1, 17], [4, 11], [1, 2]]], [492, [[3, 19], [3, 7]]], [496, [[2, 20], [4, 18], [4, 2], [2, 10], [4, 4]]], [504, [[0, 21], [0, 15]]], [505, [[3, 21], [5, 20], [5, 5], [0, 19], [2, 16], [2, 4], [2, 1], [3, 14], [5, 9], [5, 1], [0, 6]]], [508, [[1, 21], [5, 17], [1, 13], [5, 1]]], [513, [[4, 21], [1, 18], [4, 12], [4, 3], [1, 3]]], [516, [[3, 20], [3, 8]]], [520, [[2, 21], [4, 19], [2, 11], [4, 5]]], [528, [[0, 22], [0, 16]]], [532, [[1, 22], [5, 18], [5, 2], [1, 14], [5, 2]]], [537, [[4, 22], [1, 19], [4, 13], [1, 4], [1, 1]]], [540, [[3, 21], [3, 9]]], [544, [[2, 22], [4, 20], [2, 12], [4, 6]]], [552, [[0, 23], [0, 17]]], [553, [[3, 23], [5, 22], [0, 21], [2, 18], [2, 2], [3, 16], [3, 4], [3, 1], [5, 11], [0, 8], [0, 2], [2, 1]]], [556, [[1, 23], [5, 19], [1, 15], [5, 3]]], [561, [[4, 23], [1, 20], [1, 5], [4, 14], [1, 5]]], [564, [[3, 22], [3, 10]]], [568, [[2, 23], [4, 21], [2, 13], [4, 7]]], [577, [[3, 24], [3, 6], [5, 23], [0, 22], [2, 19], [3, 17], [5, 12], [5, 3], [0, 9], [0, 1], [2, 2]]], [580, [[1, 24], [5, 20], [1, 16], [5, 4]]], [585, [[4, 24], [4, 6], [1, 21], [4, 15], [1, 6]]], [588, [[3, 23], [3, 11]]], [592, [[2, 24], [4, 22], [2, 14], [4, 8]]], [600, [[0, 25], [0, 19], [0, 1]]], [601, [[3, 25], [3, 1], [5, 24], [5, 6], [0, 23], [2, 20], [2, 5], [3, 18], [3, 2], [5, 13], [0, 10], [2, 3]]], [604, [[1, 25], [1, 1], [5, 21], [1, 17], [5, 5]]], [609, [[4, 25], [4, 1], [1, 22], [4, 16], [4, 4], [4, 1], [1, 7]]], [612, [[3, 24], [3, 12]]], [616, [[2, 25], [2, 1], [4, 23], [2, 15], [4, 9], [4, 1]]], [624, [[0, 26], [0, 20], [0, 2]]], [628, [[1, 26], [5, 22], [1, 18], [1, 2], [5, 6]]], [633, [[4, 26], [1, 23], [4, 17], [1, 8], [1, 2]]], [636, [[3, 25], [3, 1], [3, 13]]], [640, [[2, 26], [4, 24], [2, 16], [4, 10]]], [648, [[0, 27], [0, 21], [0, 3]]], [649, [[3, 27], [3, 3], [5, 26], [0, 25], [0, 1], [2, 22], [3, 20], [3, 5], [5, 15], [0, 12], [0, 3], [2, 5], [3, 1]]], [652, [[1, 27], [1, 3], [5, 23], [1, 19], [5, 7]]], [657, [[4, 27], [1, 24], [1, 6], [4, 18], [1, 9]]], [660, [[3, 26], [3, 14]]], [664, [[2, 27], [2, 3], [4, 25], [4, 1], [2, 17], [4, 11]]], [672, [[0, 28], [0, 22], [0, 4]]], [673, [[3, 28], [3, 7], [5, 27], [5, 3], [0, 26], [2, 23], [3, 21], [5, 16], [5, 4], [5, 1], [0, 13], [2, 6], [3, 2]]], [681, [[4, 28], [4, 7], [1, 25], [1, 1], [4, 19], [1, 10]]], [684, [[3, 27], [3, 15]]], [688, [[2, 28], [4, 26], [2, 18], [2, 2], [4, 12]]], [696, [[0, 29], [0, 23], [0, 5]]], [697, [[3, 29], [5, 28], [5, 7], [0, 27], [0, 3], [2, 24], [2, 6], [3, 22], [5, 17], [0, 14], [2, 7], [3, 3]]], [700, [[1, 29], [5, 25], [1, 21], [5, 9], [5, 1], [1, 1]]], [705, [[4, 29], [1, 26], [4, 20], [4, 5], [1, 11]]], [708, [[3, 28], [3, 16]]], [712, [[2, 29], [4, 27], [4, 3], [2, 19], [4, 13]]], [720, [[0, 30], [0, 24], [0, 6]]], [721, [[3, 30], [5, 29], [0, 28], [0, 7], [2, 25], [2, 1], [3, 23], [5, 18], [5, 2], [0, 15], [2, 8], [2, 2], [3, 4], [3, 1]]], [724, [[1, 30], [5, 26], [1, 22], [5, 10], [1, 2]]], [732, [[3, 29], [3, 17]]], [736, [[2, 30], [4, 28], [2, 20], [4, 14]]], [744, [[0, 31], [0, 25], [0, 1], [0, 7]]], [745, [[3, 31], [5, 30], [0, 29], [2, 26], [3, 24], [3, 6], [5, 19], [0, 16], [0, 4], [0, 1], [2, 9], [2, 1], [3, 5]]], [748, [[1, 31], [5, 27], [5, 3], [1, 23], [5, 11], [1, 3]]], [753, [[4, 31], [1, 28], [1, 7], [4, 22], [1, 13], [4, 1]]], [756, [[3, 30], [3, 18]]], [760, [[2, 31], [4, 29], [2, 21], [4, 15]]], [768, [[0, 32], [0, 26], [0, 8]]], [769, [[3, 32], [3, 8], [3, 2], [5, 31], [0, 30], [2, 27], [2, 3], [3, 25], [3, 1], [5, 20], [5, 5], [0, 17], [2, 10], [3, 6]]], [772, [[1, 32], [5, 28], [1, 24], [5, 12], [1, 4]]], [777, [[4, 32], [4, 8], [4, 2], [1, 29], [4, 23], [1, 14], [4, 2]]], [780, [[3, 31], [3, 19]]], [792, [[0, 33], [0, 27], [0, 9]]], [793, [[3, 33], [5, 32], [5, 8], [5, 2], [0, 31], [2, 28], [2, 7], [3, 26], [5, 21], [0, 18], [0, 2], [2, 11], [3, 7]]], [796, [[1, 33], [5, 29], [1, 25], [1, 1], [5, 13], [1, 5]]], [801, [[4, 33], [1, 30], [4, 24], [4, 6], [1, 15], [4, 3]]], [804, [[3, 32], [3, 20]]], [808, [[2, 33], [4, 31], [2, 23], [4, 17], [2, 1]]], [816, [[0, 34], [0, 28], [0, 10]]], [817, [[3, 34], [5, 33], [0, 32], [0, 8], [0, 2], [2, 29], [3, 27], [3, 3], [5, 22], [0, 19], [2, 12], [2, 3], [3, 8], [3, 2]]], [820, [[1, 34], [5, 30], [1, 26], [5, 14], [1, 6]]], [825, [[4, 34], [1, 31], [4, 25], [1, 16], [1, 4], [1, 1], [4, 4], [4, 1]]], [828, [[3, 33], [3, 21]]], [832, [[2, 34], [4, 32], [2, 24], [4, 18], [4, 2], [2, 2]]], [840, [[0, 35], [0, 29], [0, 11]]], [844, [[1, 35], [5, 31], [1, 27], [1, 3], [5, 15], [1, 7]]], [849, [[4, 35], [1, 32], [1, 8], [1, 2], [4, 26], [1, 17], [4, 5]]], [852, [[3, 34], [3, 22]]], [856, [[2, 35], [4, 33], [2, 25], [2, 1], [4, 19], [2, 3]]], [864, [[0, 36], [0, 30], [0, 12]]], [865, [[3, 36], [3, 9], [3, 4], [3, 1], [5, 35], [0, 34], [2, 31], [3, 29], [5, 24], [5, 6], [0, 21], [2, 14], [3, 10], [5, 1]]], [868, [[1, 36], [1, 4], [5, 32], [1, 28], [5, 16], [1, 8]]], [873, [[4, 36], [4, 9], [1, 33], [4, 27], [1, 18], [4, 6]]], [876, [[3, 35], [3, 23]]], [880, [[2, 36], [2, 4], [4, 34], [2, 26], [4, 20], [2, 4]]], [888, [[0, 37], [0, 31], [0, 13]]], [889, [[3, 37], [5, 36], [5, 9], [5, 4], [5, 1], [0, 35], [2, 32], [2, 8], [2, 2], [3, 30], [5, 25], [5, 1], [0, 22], [2, 15], [3, 11], [5, 2]]], [892, [[1, 37], [5, 33], [1, 29], [5, 17], [1, 9], [1, 1]]], [897, [[4, 37], [1, 34], [4, 28], [4, 7], [1, 19], [4, 7]]], [904, [[2, 37], [4, 35], [2, 27], [2, 3], [4, 21], [2, 5]]], [912, [[0, 38], [0, 32], [0, 14]]], [913, [[3, 38], [5, 37], [0, 36], [0, 9], [0, 4], [0, 1], [2, 33], [3, 31], [5, 26], [0, 23], [2, 16], [2, 4], [2, 1], [3, 12], [3, 3], [5, 3]]], [916, [[1, 38], [5, 34], [1, 30], [5, 18], [5, 2], [1, 10]]], [921, [[4, 38], [1, 35], [4, 29], [1, 20], [1, 5], [4, 8], [4, 2]]], [924, [[3, 37], [3, 25], [3, 1], [3, 1]]], [928, [[2, 38], [4, 36], [4, 4], [2, 28], [4, 22], [2, 6]]], [936, [[0, 39], [0, 33], [0, 15]]], [937, [[3, 39], [5, 38], [0, 37], [2, 34], [3, 32], [3, 8], [3, 2], [5, 27], [5, 3], [0, 24], [0, 6], [2, 17], [3, 13], [5, 4], [5, 1]]], [940, [[1, 39], [5, 35], [1, 31], [5, 19], [1, 11]]], [945, [[4, 39], [1, 36], [1, 9], [4, 30], [1, 21], [4, 9]]], [948, [[3, 38], [3, 26], [3, 2]]], [952, [[2, 39], [4, 37], [2, 29], [4, 23], [2, 7]]], [960, [[0, 40], [0, 34], [0, 16]]], [964, [[1, 40], [5, 36], [5, 4], [1, 32], [5, 20], [1, 12]]], [969, [[4, 40], [4, 10], [1, 37], [4, 31], [1, 22], [4, 10]]], [972, [[3, 39], [3, 27], [3, 3]]], [976, [[2, 40], [4, 38], [2, 30], [4, 24], [2, 8]]], [984, [[0, 41], [0, 35], [0, 17]]], [985, [[3, 41], [5, 40], [5, 10], [0, 39], [2, 36], [2, 9], [2, 4], [2, 1], [3, 34], [5, 29], [0, 26], [2, 19], [3, 15], [5, 6], [0, 1]]], [988, [[1, 41], [5, 37], [1, 33], [5, 21], [1, 13]]], [993, [[4, 41], [1, 38], [4, 32], [4, 8], [4, 2], [1, 23], [4, 11]]], [996, [[3, 40], [3, 28], [3, 4]]], [1000, [[2, 41], [4, 39], [2, 31], [4, 25], [2, 9], [2, 1]]]];

for Lev in Levels:
	D=Lev[0];
	print("\n D = ",D);
	for N in Lev[1]:
		null=FNvsGamma0(D,N[1],N[0]);
