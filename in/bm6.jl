#   PRIMER: KONZOLA Z MOMENTOM NA ROBU

#   K O O R D I N A T E   V O Z L I S C
vozlisca::Array{Float64} = [
    0. 0.;
    0. -10/pi;
    0. 0.;
    ]# * [1 0; -0.005 1] imperfektnosti v x glede na z


#   E L E M E N T I   M E D   V O Z L I S C I
elementi::Array{Int64} = [
    1 2;
    2 3
    ]



#   P O D A T K I   R A Č U N A
const ti::Float64 = 0.0
const dt::Float64 = 0.01
const tf::Float64 = 1.0
const g::Vector{Float64}  = [0.; 0.]

metoda_t_integracije::String    = ["midpoint", "timeelementP","timeelementT"][1]
tnodes                          = [0.;0.5;1.]
Integracija::String 	        = ["gauss", "lobatto"][2]
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
@assignto :(ElementDataIn) [1,2] :( [0. ; 0.] ) :(M)
@assignto :(ElementDataIn) [1,2] :( 10^4*[1. 0. 0.;0. 1. 0.; 0. 0. 0.05] ) :(C)


#@assignto :(ElementDataIn) [1] :(t->[0.1, 0.1]*t) :(px)
#@assignto :(ElementDataIn) [1] :(t->[repeat([0.],19);8.]*Int(t<=2.5) ) :(Px)
#@assignto :(ElementDataIn) [1] :(t->[5.  5.]*t  ) :(pz)
#@assignto :(ElementDataIn) [1] :(t->[0.;50.0*t*Int(t<0.5)]) :(Pz)
#@assignto :(ElementDataIn) [1] :(t->[0., 0.]  ) :(my)
@assignto :(ElementDataIn) [2] :(t->[repeat([0.],3);-100*pi*t]) :(My)


@assignto :(ElementDataIn) [1,2] :( range(-1,1,length=3) |> collect ) :(div1)
@assignto :(ElementDataIn) [1,2] :( repeat([6],2) ) :(div2)
#@assignto :(ElementDataIn) [1] :( :chebyshev2 ) :(dist)
@assignto :(ElementDataIn) [1,2] :( repeat([15],2) ) :(nInt)
#@assignto :(ElementDataIn) [1] :( true ) :(Ci)
#@assignto :(ElementDataIn) [1] :($(:(0.1))) :(beta)
#ElementDataIn[1].beta = 0.5
#ElementDataIn[1].Ci=2


@assignto :(ElementDataIn) [1,2] :( re_gramschmid([[-1.,1.,0.],[-1.,0.,1.]])) :(Ib_geom)
@assignto :(ElementDataIn) [1] :( 10/pi*[ 0.5 -0.5;1.0*pi/4  0.0; 0. -1.0*pi/4; -1.0*pi/4 0.0] ) :(Kb)
@assignto :(ElementDataIn) [2] :( 10/pi*[ -0.5 -0.5;-1.0*pi/4 0.0; 0. 1.0*pi/4; 1.0*pi/4 0.0] ) :(Kb)


# V O Z L I Š Č A
@assignto :(VozDataIn) [1] :( Bool[0, 0, 0] ) :(Supp)
#@assignto :(VozDataIn) [1] :( t->[0.; 0.; (t <= 5.0 ? 1.5/5.0*t : 1.5)] ) :(mot)


#@assignto :(VozDataIn) [2] :( Bool[1, 0, 1] ) :(Supp)

#@assignto :(VozDataIn) [1] :( pi/3. ) :(dir)


# macro iteration_hook()
# nothing
# end;
#
# macro precompute_hook()
# nothing
# end;

#=
macro postcompute_hook()
    nothing
end;
=#




proc_check=false


include("../NonLinBeamRUN.jl")
#=
plt1 = plotVar3(M,E,ElementDataIn,VozDataIn,[1;1:20:101|>collect];init_konf=false,linecolor=:black,label=:none,aspectratio=:equal,minorgrid=false,yflip=true,yrange=(-8,1),xrange=(-3,11))

plt2 = plot(vcat(map(i-> E[1].xInt[1] .+ [0.;cumsum(E[1].L)][i],1:length(E[1].indx))...),M.gamma3[:,end],label = false,linecolor = :black, xticks = 0:2.5:10|>collect,minorgrid=false)

A1 = [M.ux[[2,5,8,11],end]  M.uz[[2,5,8,11],end]]




@assignto :(ElementDataIn) [1] :( range(-1,1,length=5) |> collect ) :(div1)
@assignto :(ElementDataIn) [1] :( [4;repeat([2],2);4] ) :(div2)
@assignto :(ElementDataIn) [1] :( repeat([15],4) ) :(nInt)
ElementDataIn[1].Ci=2


include("../NonLinBeamRUN.jl")

plotVar3(M,E,ElementDataIn,VozDataIn,[1;1:20:101|>collect];init_konf=false, p0=plt1,linecolor=:lawngreen,line=:dash,label=:none,aspectratio=:equal,minorgrid=false,yflip=true,yrange=(-8,1),xrange=(-3,11))

plt2 = plot(plt2,vcat(map(i-> E[1].xInt[1] .+ [0.;cumsum(E[1].L)][i],1:length(E[1].indx))...),M.gamma3[:,end],label = false,linecolor = :lawngreen,line=:dash, xticks =0:2.5:10|>collect,minorgrid=false)






@assignto :(ElementDataIn) [1] :( range(-1,1,length=5) |> collect ) :(div1)
@assignto :(ElementDataIn) [1] :( [5,6,6,5] ) :(div2)
@assignto :(ElementDataIn) [1] :( repeat([15],4) ) :(nInt)
ElementDataIn[1].Ci=0

A2 = [M.ux[[2,4,7,10],end]  M.uz[[2,4,7,10],end]]

include("../NonLinBeamRUN.jl")

plotVar3(M,E,ElementDataIn,VozDataIn,[1;1:20:101|>collect];init_konf=false, p0=plt1,linecolor=:crimson,line=:dot,label=:none,aspectratio=:equal,minorgrid=false,yflip=true,yrange=(-8,1),xrange=(-3,11))

plt2 = plot(plt2,vcat(map(i-> E[1].xInt[1] .+ [0.;cumsum(E[1].L)][i],1:length(E[1].indx))...),M.gamma3[:,end],label = false,linecolor = :crimson,line=:dot, xticks =0:2.5:10|>collect,minorgrid=false)


plt1=plot(plt1,[0,0],[0,0],label = latexify("4^3C^0"), ylabel = latexify("z"),yguidefontrotation=-90,xlabel = latexify("x"),lc=:black)
plt1=plot(plt1,[0,0],[0,0],label = latexify("4^5C^2"),lc=:lawngreen,line=:dash)
plt1=plot(plt1,[0,0],[0,0],label = latexify("4^5C^0"),lc=:crimson,line=:dot)


plt2=plot(plt2,[0,0],[0,0],label = latexify("4^3C^0"), ylabel = latexify("K_2"),yguidefontrotation=-90,xlabel = latexify("s"),lc=:black)
plt2=plot(plt2,[0,0],[0,0],label = latexify("4^5C^2"),lc=:lawngreen,line=:dash)
plt2=plot(plt2,[0,0],[0,0],label = latexify("4^5C^0"),lc=:crimson,line=:dot)
plt3 = plot(plt2,yrange=(0.58,0.69),xrange=(2.5,5.0))


#savefig(plt1,"~/0_git/dinamika_nosilcev/out/bm5_konf.tex")
#savefig(plt2,"~/0_git/dinamika_nosilcev/out/bm5_K2.tex")
#savefig(plt3,"~/0_git/dinamika_nosilcev/out/bm5_K22.tex")
A3 = [M.ux[[2,7,12,17],end]  M.uz[[2,7,12,17],end]]

=#
