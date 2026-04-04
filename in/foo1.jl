#   K O O R D I N A T E   V O Z L I S C
vozlisca::Array{Float64} = [
    0.0  0.0;
    1.0  0.0;
    ]# * [1 0; -0.005 1] imperfektnosti v x glede na z


#   E L E M E N T I   M E D   V O Z L I S C I
elementi::Array{Int64} = [
    1 2;
    ]



#   P O D A T K I   R A Č U N A
const ti::Float64        = 0.0
const dt::Float64        = 0.1
const tf::Float64        = 10.0
const g::Vector{Float64} = [0.; 0.]

metoda_t_integracije::String    = ["midpoint", "timeelementP","timeelementT"][1]
tnodes                          = [0.;0.5;1.]
Integracija::String 	        = ["gauss", "lobatto"][1]
nt = 2

const dv_norm_tol_exp::Int64	   = -7
const nwt_iter_max_count::Int64	   = 30








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
#   dist    - Symbol            ()          :unifom
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
@assignto :(ElementDataIn) [1] :( [1.; 10.] ) :(M)
@assignto :(ElementDataIn) [1] :( 10^4*[1. 0. 0.;0. 1. 0.; 0. 0. 0.01] ) :(C)


#@assignto :(ElementDataIn) [1] :(t->[0.1, 0.1]*t) :(px)
#@assignto :(ElementDataIn) [1] :(t->[repeat([0.],1);8.0]*Int(t<=2.5) ) :(Px)
@assignto :(ElementDataIn) [1] :(t->[-500.  500.]*t  ) :(pz)
#@assignto :(ElementDataIn) [1] :(t->[repeat([0.],1);8.0*8/10]*Int(t<=2.5)) :(Pz)
#@assignto :(ElementDataIn) [1] :(t->[0., 0.]  ) :(my)
#@assignto :(ElementDataIn) [1] :(t->[repeat([0.],1);-80.]*Int(t<=2.5)) :(My)


@assignto :(ElementDataIn) [1] :( range(-1,1,length=2) |> collect ) :(div1)
@assignto :(ElementDataIn) [1] :( repeat([5],1) ) :(div2)
@assignto :(ElementDataIn) [1] :( :chebyshev2 ) :(dist)
@assignto :(ElementDataIn) [1] :( repeat([7],1) ) :(nInt)
#@assignto :(ElementDataIn) [1] :( true ) :(Ci)


#@assignto :(ElementDataIn) [1] :( re_gramschmid([[-1.,1.,0.]])) :(Ib_geom)
#@assignto :(ElementDataIn) [1] :( [2.5 -0.5] ) :(Kb)


# V O Z L I Š Č A
@assignto :(VozDataIn) [1,2] :( Bool[0, 0, 1] ) :(Supp)
#@assignto :(VozDataIn) [1] :( t->[0, 0, 0] ) :(mot)

#@assignto :(VozDataIn) [2] :( Bool[1, 0, 1] ) :(Supp)

#@assignto :(VozDataIn) [1] :( pi/3. ) :(dir)





