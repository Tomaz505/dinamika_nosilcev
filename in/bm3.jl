#   PRIMER: LOK Z ALTERNIRAJOCO OBTEZBO

#   K O O R D I N A T E   V O Z L I S C
vozlisca::Array{Float64} = [
    -5. 0.;
    5. 0.;
    ]# * [1 0; -0.005 1] imperfektnosti v x glede na z


#   E L E M E N T I   M E D   V O Z L I S C I
elementi::Array{Int64} = [
    1 2;
    ]



#   P O D A T K I   R A Č U N A
const ti::Float64 = 0.0
const dt::Float64 = 0.0001
const tf::Float64 = 0.06
const g::Vector{Float64}  = [0.; 0.]

metoda_t_integracije::String    = ["midpoint", "timeelementP","timeelementT"][1]
tnodes                          = [0.;0.5;1.]
Integracija::String 	        = ["gauss", "lobatto"][1]
nt = 2

const dv_norm_tol_exp::Int64	   = -8
const nwt_iter_max_count::Int64	   = 150








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
@assignto :(ElementDataIn) [1] :( 7850*[0.087; 0.003562] ) :(M)
@assignto :(ElementDataIn) [1] :( 2.1e11*[0.087 0. 0.;0. 0.087/2.6 0.; 0. 0. 0.003562] ) :(C)


#@assignto :(ElementDataIn) [1] :(t->[0.1, 0.1]*t) :(px)
#@assignto :(ElementDataIn) [1] :(t->[repeat([0.],19);8.]*Int(t<=2.5) ) :(Px)
#@assignto :(ElementDataIn) [1] :(t->[5.  5.]*t  ) :(pz)
@assignto :(ElementDataIn) [1] :(t->[zeros(7);8.0e7*sin(1000*t);zeros(8) ]) :(Pz)
#@assignto :(ElementDataIn) [1] :(t->[0., 0.]  ) :(my)
#@assignto :(ElementDataIn) [1] :(t->[repeat([0.],19);-80.]*Int(t<=2.5)) :(My)


@assignto :(ElementDataIn) [1] :( range(-1,1,length=9) |> collect ) :(div1)
@assignto :(ElementDataIn) [1] :( repeat([4],8) ) :(div2)
#@assignto :(ElementDataIn) [1] :( repeat([5],8) ) :(div2)

#@assignto :(ElementDataIn) [1] :( :chebyshev2 ) :(dist)
@assignto :(ElementDataIn) [1] :( repeat([15],8) ) :(nInt)
#@assignto :(ElementDataIn) [1] :( true ) :(Ci)
#ElementDataIn[1].Ci =2

@assignto :(ElementDataIn) [1] :( re_gramschmid([[-1.,1.,-2/3,-1/3,1/3,2/3]])) :(Ib_geom)
@assignto :(ElementDataIn) [1] :( [-10*sin(pi/9) -10*cos(pi/9)+sqrt(75); -10*sin(pi/18) -10*cos(pi/18)+sqrt(75) ;10*sin(pi/18) -10*cos(pi/18)+sqrt(75);10*sin(pi/9) -10*cos(pi/9)+sqrt(75) ] ) :(Kb)


# V O Z L I Š Č A
@assignto :(VozDataIn) [1,2] :( Bool[0, 0, 0] ) :(Supp)
#@assignto :(VozDataIn) [1] :( t->[0.; 0.; (t <= 5.0 ? 1.5/5.0*t : 1.5)] ) :(mot)
#@assignto :(VozDataIn) [2] :( Bool[1, 0, 1] ) :(Supp)
#@assignto :(VozDataIn) [1] :( pi/3. ) :(dir)



#= OUTPUT

2x
plot(time_st,(M.gamma3[20,:]-M.gamma3[21,:]);linecolor = :black, legend = :none,xlabel = latexify("t"),tickfontsize=16,xguidefontsize=16,yguidefontsize=16,size=(800,400))

uz1=M.uz[14,:]
uz2=M.uz[10,:]
plot(time_st,(uz1-uz2);linecolor = :black, legend = :none,xlabel = latexify("t"),tickfontsize=16,xguidefontsize=16,yguidefontsize=16,size=(800,400))

plot(time_st,M.uz[14,:];linecolor = :black, legend = :none,xlabel = latexify("t"),tickfontsize=16,xguidefontsize=16,yguidefontsize=16,size=(800,400))
=#


# macro iteration_hook()
# nothing
# end;
#
# macro precompute_hook()
# nothing
# end;

# macro postcompute_hook()
# nothing
# end;


proc_check=true
# include("../NonLinBeamRUN.jl")
#=
plt1 = plot(time_st,M.uz[14,:],minorgrid=false,yflip=true,lc=:black,label=latexify("8C^0_3"))


plt2=plot(vcat(map(i-> E[1].xInt[1] .+ [0.;cumsum(E[1].L)][i],1:length(E[1].indx))...),M.gamma3[:,end],label = false,linecolor = :black, xticks = 0:2.5:10|>collect,minorgrid=false)=#

ElementDataIn[1].div2 = [3,2,2,2,2,2,2,3]
ElementDataIn[1].Ci=2


include("../NonLinBeamRUN.jl")

plt1 = plot(plt1,time_st,M.uz[13,:],minorgrid=false,yflip=true,lc=:lawngreen, line=:dash,label=latexify("8C^2_5"))

plt2=plot(plt2,vcat(map(i-> E[1].xInt[1] .+ [0.;cumsum(E[1].L)][i],1:length(E[1].indx))...),M.gamma3[:,end],label = false,linecolor = :lawngreen, line=:dash, xticks = 0:2.5:10|>collect,minorgrid=false)

#plt1 = plot(time_st,M.uz[14,:],minorgrid=false,yflip=true,lc=:black,label=latexify("8C^0"))

#savefig(plt1,"~/0_git/dinamika_nosilcev/out/bm3_uz.tex")


