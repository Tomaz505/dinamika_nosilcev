#   K O O R D I N A T E   V O Z L I S C
vozlisca::Array{Float64} = [
        0. 0.;
	0. 1.
    ]


#   E L E M E N T I   M E D   V O Z L I S C I
elementi::Array{Int64} = [
        1 2
    ]



#   P O D A T K I   R A Č U N A 
const ti = 0.
const tf = 2.
const dt = 0.3
const g = 0.0


n_elem,n_voz,ElementDataIn,VozDataIn = datainit(elementi,vozlisca)
#	L A S T N O S T I   S T R U K T U R
#	K T U R

# ElementDataIn
# n = st. ke. na nosilcui
# m = st. tock geometrije
#
#	param	tip			oblika		default
#
#	v	- Vector{Int64}		(2)		...
#	C 	- Matrix{Float64} 	(3x3)		
#	M 	- Matrix{Float64} 	(nx2)		
#	div1	- Vector{Float64} 	(n+1)		[-1.; 1,]
#	div2	- Vector{Int64} 	(n)		[4]
#	nInt	- Vector{Int64}   	(n)		[20]
#	Ci	- Bool			()		false
#	pz	- Function		t->(nx2)	t->[0. 0.]
#	px	- Function		t->(nx2)	t->[0. 0.]
#	my	- Function		t->(nx2)	t->[0. 0.]
#	Ib_geom	- Matrix{FLoat64}	(m,m)		[0.5 0.5; -0.5 0.5] == re_gramshchim([-1.;1.])
#	Kb	- Matrix{Float64}	(m-2x2)		[]

# VozDataIn
#
#	x	- Float64		()		...
#	z	- Float64		()		...
#	i	- Int64			()		...
#	Supp	- Vector{Bool}		(3)		[false,false,false]
#	dir	- Float64		()		0.0

















# E L E M E N T I
@assignto :(ElementDataIn) [1] :( [0.0 0.0] ) :(M)
@assignto :(ElementDataIn) [1] :( [21000. 0. 0.;0. 17500. 0.; 0. 0. 0.175] ) :(C)


#@assignto :(ElementDataIn) [1] :(t->[0.0 0.0]) :(pz)
@assignto :(ElementDataIn) [1] :(t->[0.0001 0.0001]) :(px)
#@assignto :(ElementDataIn) [1] :(t->[0. 0.]) :(my)

#@assignto :(ElementDataIn) [1] :( [-1.; 1.] ) :(div1)
@assignto :(ElementDataIn) [1] :( [5] ) :(div2)
@assignto :(ElementDataIn) [1] :([20] ) :(nInt) 
#@assignto :(ElementDataIn) [1] :( false ) :(Ci)


#@assignto :(ElementDataIn) [1] :( re_gramschmid([[-1.,1.,0.]])) :(Ib_geom)
#@assignto :(ElementDataIn) [1] :( [0.5 0.5] ) :(Kb) 


# V O Z L I Š Č A
@assignto :(VozDataIn) [1] :( Bool[0, 0, 0] ) :(Supp)
#@assignto :(VozDataIn) [1] :( pi/3. ) :(dir)





