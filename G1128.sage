def mult(v,w):
    return vector([v[0]*w[0]-v[1]*w[1], v[0]*w[1]+v[1]*w[0]]);

def G1128(l=[5,1],place_precision=150,spin=0, NoPrint=True, Labels=False):
	"""
	Creates the flat surface associated to a quadrilateral of type (1,1,2,8), with side [2,8] of length l
	Input: 
	* l : length of the side [2,8], given by (sqrt(l[0])+l[1])/2
	* place_precision (opt): integer
	Output:
	* TS : triangulated surface corresponding to the flat surface associated to a quadrilateral of type (1,1,2,8).
	* K : the field Q(sqrt(l[0]))
	* mD : square root gen(Q(l))
	* pl : the homomorphism from K to the real field of precision place_precision with pl(mD) > 0
	"""	
	K.<mD> = NumberField(x^2-l[0]);
	l=(mD+l[1])/2;
	pl = [ p for p in K.places(prec=place_precision) if (p(mD)> 0)][0];
	
	Vaux=[ vector([1,0]), vector([1/2,1/2]), vector([-1/2,1/2]), vector([-1,0]), vector([-1/2,-1/2]), vector([1/2,-1/2]) ];
	Vaux1=[ vector([-1/2, 1/2]), vector([-1, 0]), vector([-1/2, -1/2]), vector([1/2, -1/2]), vector([1, 0]), vector([1/2, 1/2]) ];
	Vaux2=[ vector([1/2*3, 1/2]), vector([0, 1]), vector([-1/2*3, 1/2]), vector([-1/2*3, -1/2]), vector([0, -1]), vector([1/2*3, -1/2]) ];
	Vaux3=[ vector([-1/2, -1/2]), vector([1/2, -1/2]), vector([1, 0]), vector([1/2, 1/2]), vector([-1/2, 1/2]), vector([-1, 0]) ];

	V1 = [ (2*l+1)*v for v in Vaux];
	V2 = [ (2*l+1)*Vaux[k] + l*Vaux1[k] for k in [0..5]];
	V3 = [ 1/2*(2*l+2)*v for v in Vaux2];
	V4 = [ (2*l+1)*Vaux[k] + l*Vaux3[k] for k in [0..5]];
	V5 = [ (3*l/2+1)*v for v in Vaux];

	a=Permutation('(1,2,3,4,5,6)'); # renumbering V4
	V4=[V4[a(i)-1] for i in [1..6]];
	W1=[v.list() for v in V1]; W2=[v.list() for v in V2]; W3=[v.list() for v in V3]; W4=[v.list() for v in V4]; W5=[v.list() for v in V5]; 
	O=[0,0];
	
	Laux1 = [ [ W2[j], W4[j], W4[(j-1)%6] ]  for j in [0..5] ];
	Laux2 = [ [ W4[2*j], W4[2*j+1], W4[(2*j+2)%6] ]  for j in [0..2] ];
	Laux3 = [ [ W4[0], W4[2], W4[4] ] ];
	Laux4 = [ [ W4[0], (V1[1]+V4[2]-V1[3]).list(), W2[1] ] ];
	Laux5 = [ [ W2[4], W4[3], (V1[4]+V4[5]-V1[0]).list() ] ];
	Laux6 = [ [ W2[2], (V2[2]+Vaux[2]).list(), W4[2] ], [(V2[2]+Vaux[2]).list(), (V2[2]+Vaux[2]+Vaux[3]).list(), W4[2] ], [(V2[2]+Vaux[2]+Vaux[3]).list(), (V2[2]+Vaux[2]+Vaux[3]+Vaux[4]).list(), W4[2] ], [(V2[2]+Vaux[2]+Vaux[3]+Vaux[4]).list(), (V2[2]+Vaux[2]+Vaux[3]+Vaux[4]+Vaux[5]).list(), W4[2] ] ] ;

	Ts= Laux1 + Laux2 + Laux3 + Laux4 + Laux5 + Laux6;
	
	G1=[ [[2*k,0],[6+((k-1)%3),0]] for k in [0..2] ] + [ [[2*k+1,0],[6+k,2]] for k in [0..2] ] ;
	G2=[ [[6,1],[9,2]], [[7,1],[9,0]], [[8,1],[9,1]] ];
	G3=[ [[2*k,1],[11,k]] for k in [0..2] ] + [ [[2*k+1,1],[10,(k+1)%3]] for k in [0..2] ] ;
	G4=[ [[0,2],[15,2]], [[1,2],[15,0]], [[2,2],[12,1]], [[3,2],[12,2]], [[4,2],[13,2]], [[5,2],[14,2]], [[12,0],[13,1]], [[13,0],[14,1]], [[14,0],[15,1]] ];
	Gs = G1+G2+G3+G4;
	
	Points = [];
	Edges = plot([]);
	eps=0.2;
	if Labels:
		txt=sum(plot(text(str(2*j),(pl(Ts[2*j][0][0])+eps,pl(Ts[2*j][0][1])+eps),fontsize=10, aspect_ratio=1, horizontal_alignment="center", vertical_alignment="center")) for j in [0..len(Ts)/2-1]);
		txt+=sum(plot(text(str(2*j+1),(pl(Ts[2*j+1][0][0])-eps,pl(Ts[2*j+1][0][1])-eps),fontsize=10, aspect_ratio=1, horizontal_alignment="center", vertical_alignment="center")) for j in [0..len(Ts)/2-1]);
	else:
		txt=plot([]);
	for t in Ts:
		Edges+=line([(pl(t[0][0]),pl(t[0][1])),(pl(t[1][0]),pl(t[1][1]))]);
		Edges+=line([(pl(t[1][0]),pl(t[1][1])),(pl(t[2][0]),pl(t[2][1]))]);
		Edges+=line([(pl(t[2][0]),pl(t[2][1])),(pl(t[0][0]),pl(t[0][1]))]);
		Points.append([pl(t[0][0]),pl(t[0][1])]);
		Points.append([pl(t[1][0]),pl(t[1][1])]);
		Points.append([pl(t[2][0]),pl(t[2][1])]);
	P=point2d(Points,size=20)+Edges+txt;
	if (not(NoPrint)):
		P.show(aspect_ratio=1,title='Translation surface TS');
	Ts = [ [vector(v) for v in t] for t in Ts]
	Ts = [ triangle(t,pl) for t in Ts]
	TS = triangulated_surface(Ts,Gs,pl)
	TS = TS.delaunay()
	return [TS, K, mD, pl]
