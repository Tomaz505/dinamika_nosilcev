#   K O O R D I N A T E   V O Z L I S C
vozlisca::Array{Float64} = [
	0. 0.;
	0. -5.;
	0. -8.5;
	0. -12.;
	4. 0.;
	4. -5.;
	4. -8.5;
	4. -12.;
    ]# * [1 0; -0.005 1] imperfektnosti v z


#   E L E M E N T I   M E D   V O Z L I S C I
elementi::Array{Int64} = [
        1 2;
        2 3;
        3 4;
        5 6;
        6 7;
        7 8;

        1 6;
        2 7;
        3 8;
        5 2;
        6 3;
        7 4;

        2 6;
        3 7;
        4 8;
    ]



#   P O D A T K I   R A Č U N A
const ti::Float64 = 0.0
const dt::Float64 = 0.01
const tf::Float64 = 5.
const g::Vector{Float64}  = [0.; 0.]


#	K O N T R O L N I   P A R A M E T R I
                                      #1          #2
metoda_t_integracije::String    = ["midpoint", "timeelement"][1]
#tnodes                         = QuadInt(3,mtd = "chebyshev")[2] #Pomembno le za "timeelement"
Integracija::String 	        = ["gauss", "lobatto"][1]
nt = 2

const dv_norm_tol_exp::Int64	   = -6
const nwt_iter_max_count::Int64	   = 170








n_elem,n_voz,ElementDataIn,VozDataIn = datainit(elementi,vozlisca)
#	L A S T N O S T I   S T R U K T U R

# ElementDataIn
# n = st. ke. na nosilcu
# m = st. tock geometrije
#
#	param	tip			oblika		default
#
#	v	- Vector{Int64}		(2)		...
#	C 	- Matrix{Float64} 	(3x3)
#	M 	    - Matrix{Float64} 	(nx2)
#	div1	- Vector{Float64} 	(n+1)		[-1.; 1,]
#	div2	- Vector{Int64} 	(n)		    [4]
#	nInt	- Vector{Int64}   	(n)		    [20]
#	Ci	    - Bool			    ()		    false
#	pz	    - Function		    t->(nx2)	t->nothing
#	px	    - Function		    t->(nx2)	t->nothing
#	my	    - Function		    t->(nx2)	t->nothing
#	Px	    - Function		    t->(n+1)	t->nothing
#	Pz	    - Function		    t->(n+1)	t->nothing
#	My	    - Function		    t->(n+1)	t->nothing
#	Ib_geom	- Matrix{FLoat64}	(m,m)		[0.5 0.5; -0.5 0.5] == re_gramshchim([-1.;1.])
#	Kb	    - Matrix{Float64}	(m-2x2)		[]

# VozDataIn
#
#	x	- Float64		()		...
#	z	- Float64		()		...
#	i	- Int64			()		...
#	Supp	- Vector{Bool}		(3)		[false,false,false]
#	dir	- Float64		()		0.0

















# E L E M E N T I
@assignto :(ElementDataIn) [1,2,3,4,5,6] :( 7.85*[0.0131; 0.0001927] ) :(M)

@assignto :(ElementDataIn) [7,10] :( 7.85*[0.004073; 8.619e-6] ) :(M)
@assignto :(ElementDataIn) [8,11] :( 7.85*[0.00227; 3.745e-6] ) :(M)
@assignto :(ElementDataIn) [9,12] :( 7.85*[0.001226; 1.453e-6] ) :(M)

@assignto :(ElementDataIn) [13,14,15] :( 7.85*[1.0; 0.0000194] ) :(M)


@assignto :(ElementDataIn) [1,2,3,4,5,6] :( 210*10^6*[0.0131 0. 0.;0. 0.0131/2.6 0.; 0. 0. 0.0001927] ) :(C)

@assignto :(ElementDataIn) [7,10] :( 210*10^6*[0.004073 0. 0.;0. 0.004073/2.6 0.; 0. 0. 8.619e-6] ) :(C)
@assignto :(ElementDataIn) [8,11] :( 210*10^6*[0.00227 0. 0.;0. 0.00227/2.6 0.; 0. 0. 3.745e-6] ) :(C)
@assignto :(ElementDataIn) [9,12] :( 210*10^6*[0.001226 0. 0.;0. 0.001226/2.6 0.; 0. 0. 1.463e-6] ) :(C)

@assignto :(ElementDataIn) [13,14,15] :( 210*10^6*[0.00285 0. 0.;0. 0.00285/2.6 0.; 0. 0. 0.0000194] ) :(C)


#@assignto :(ElementDataIn) [1] :(t->[0.1, 0.1]*t) :(px)
#@assignto :(ElementDataIn) [2] :(t->[0.;400.]*t ) :(Px)
#@assignto :(ElementDataIn) [1] :(t->[0., 0.]  ) :(pz)
#@assignto :(ElementDataIn) [1] :(t->[0.;0.;0.;0.;0.;100. *t^2] ) :(Pz)
#@assignto :(ElementDataIn) [1] :(t->[0., 0.]  ) :(my)
#@assignto :(ElementDataIn) [1] :(t->[0.;0.;0.;0.;0.;10.]*t^2) :(My)


#@assignto :(ElementDataIn) [1,2] :( [-1.0; 0.; 1.0] ) :(div1)
#@assignto :(ElementDataIn) [1,2] :( [4;4] ) :(div2)
#@assignto :(ElementDataIn) [1,2] :( [6;6] ) :(nInt)
#@assignto :(ElementDataIn) [1] :( true ) :(Ci)


#@assignto :(ElementDataIn) [1] :( re_gramschmid([[-1.,1.,0.]])) :(Ib_geom)
#@assignto :(ElementDataIn) [1] :( [2.5 -0.1] ) :(Kb)
@assignto :(ElementDataIn) [7,8,9,10,11,12,13,14,15] :( (true,true) ) :(relese)


# V O Z L I Š Č A
@assignto :(VozDataIn) [1,5] :( Bool[0, 0, 1] ) :(Supp)
@assignto :(VozDataIn) [2,3,4,6,7,8] :( Bool[0, 1, 1] ) :(Supp)

#@assignto :(VozDataIn) [2] :( Bool[1, 0, 0] ) :(Supp)
@assignto :(VozDataIn) [2,6] :( t->[2.0/12*5/5.0*t;0.;0.] ) :(mot)
@assignto :(VozDataIn) [3,7] :( t->[2.0/12*8.5/5.0*t;0.;0.] ) :(mot)
@assignto :(VozDataIn) [4,8] :( t->[2.0/5.0*t;0.;0.] ) :(mot)
#@assignto :(VozDataIn) [1,5] :( t->[0.2*sin(pi*15*t);0.;0.] ) :(mot)

#@assignto :(VozDataIn) [1] :( pi/3. ) :(dir)





#plot(M.ux[2,:],load;linez=M.vx[2,:],label = :none,xlabel = latexify("u_x"),ylabel = latexify("F_x"),colorbar_title = latexify("v_x"))
