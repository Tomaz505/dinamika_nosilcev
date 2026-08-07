#   PRIMER: LETEC SPAGET

#   K O O R D I N A T E   V O Z L I S C
vozlisca::Array{Float64} = [
    0.0 -8.0;
    3.0  -4.0;
    6.0 0.0
    ]# * [1 0; -0.005 1] imperfektnosti v x glede na z


#   E L E M E N T I   M E D   V O Z L I S C I
elementi::Array{Int64} = [
    1 2;
    2 3;
    ]



#   P O D A T K I   R A Č U N A
const ti::Float64        = 0.0
const dt::Float64        = 0.005
const tf::Float64        = 100.0
const g::Vector{Float64} = [0.; 0.]

metoda_t_integracije::String    = ["midpoint", "timeelementP","timeelementT"][1]
tnodes                          = [0.;0.5;1.]
Integracija::String 	        = ["gauss", "lobatto"][2]
nt = 2

const dv_norm_tol_exp::Int64	   = -8
const nwt_iter_max_count::Int64	   = 30








n_elem,n_voz,ElementDataIn,VozDataIn = datainit(elementi,vozlisca)
#	L A S T N O S T I   S T R U K T U R

# ElementDataIn
# n = st. ke. na nosilcu
# m = st. tock geometrije
#
#	param	tip			oblika		default
#
#	v	    - Vector{Int64}		(2)		    ...
#	C 	    - Matrix{Float64} 	(3x3)       10^4*[1. 0. 0.;0. 1. 0.;0. 0. 0.1]
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
@assignto :(ElementDataIn) [1] :( [1.; 1.] ) :(M)
@assignto :(ElementDataIn) [2] :( [1.; 10.] ) :(M)
@assignto :(ElementDataIn) [1,2] :( 1.0e6*[1. 0. 0.;0. 1. 0.; 0. 0. 0.001] ) :(C)


#@assignto :(ElementDataIn) [1] :(t->[0.1, 0.1]*t) :(px)
#@assignto :(ElementDataIn) [1] :(t->[5.  5.]*t  ) :(pz)
#@assignto :(ElementDataIn) [1] :(t->[0., 0.]  ) :(my)
@assignto :(ElementDataIn) [2] :(t->[0.0;0.0;0.0;0.0;0.0;0.0;0.0;40.0*Int(t<=2.5)] ) :(Px)
#@assignto :(ElementDataIn) [1] :(t->[repeat([0.],1);8.0*8/10]*Int(t<=2.5)) :(Pz)
@assignto :(ElementDataIn) [2] :(t->[0.0;0.0;0.0;0.0;0.0;0.0;0.0;-160.0*Int(t<=2.5)]) :(My)


@assignto :(ElementDataIn) [1,2] :( range(-1,1,length=5) |> collect ) :(div1)
@assignto :(ElementDataIn) [1,2] :( repeat([4],4) ) :(div2)
#@assignto :(ElementDataIn) [1] :( repeat([4],30) ) :(div2)
#@assignto :(ElementDataIn) [1] :( :chebyshev2 ) :(dist)
@assignto :(ElementDataIn) [1,2] :( repeat([8],4) ) :(nInt)
@assignto :(ElementDataIn) [2] :( (true,false) ) :(relese)
#@assignto :(ElementDataIn) [1] :( true ) :(Ci)
#ElementDataIn[1].Ci = 2

#@assignto :(ElementDataIn) [1] :( re_gramschmid([[-1.,1.,0.]])) :(Ib_geom)
#@assignto :(ElementDataIn) [1] :( [2.5 -0.5] ) :(Kb)


# V O Z L I Š Č A
#@assignto :(VozDataIn) [1] :( Bool[0, 0, 1] ) :(Supp)
#@assignto :(VozDataIn) [1] :( t->[0, 0, 0] ) :(mot)

#@assignto :(VozDataIn) [2] :( Bool[1, 0, 1] ) :(Supp)

#@assignto :(VozDataIn) [1] :( pi/3. ) :(dir)

#=
macro iteration_hook()
nothing
end;
=#

#=
macro precompute_hook()
nothing
end;
=#



proc_check=true


include("../NonLinBeamRUN.jl")
#pltW = plot(time_st[1:20:end],energija(M,ElementDataIn,E,VozDataIn)[1:20:end];lc=:dodgerblue,label = latexify("EI=1000"));

plt11 = plotVar3(M,E,ElementDataIn,VozDataIn,[1;1:50:1001|>collect];init_konf=false,np=5, linecolor=:dodgerblue, label = :none,aspectratio=:equal,minorgrid=false,yflip=true);

plt12 = plotVar3(M,E,ElementDataIn,VozDataIn,[19401;19401:50:20001|>collect];init_konf=false,np=5,linecolor=:dodgerblue,label=:none,aspectratio=:equal,minorgrid=false,yflip=true);

#pltUz1=plot(time_st[1:20:end],M.uz[1,1:20:end],lc=:dodgerblue,xlabel=latexify("t"), ylabel = latexify("u_z"),yguidefontrotation=-90,label = latexify("EI=1000"));


@assignto :(ElementDataIn) [1,2] :( 1.0e6*[1. 0. 0.;0. 1. 0.; 0. 0. 0.0001] ) :(C)
include("../NonLinBeamRUN.jl")
#pltW = plot(pltW,time_st[1:20:end],energija(M,ElementDataIn,E,VozDataIn)[1:20:end];lc=:lawngreen,label = latexify("EI=100"),xlabel=latexify("t"), ylabel = latexify("W"),yguidefontrotation=-90);

plt21 = plotVar3(M,E,ElementDataIn,VozDataIn,[1;1:50:1001|>collect];init_konf=false,np=10,linecolor=:lawngreen,label = :none,aspectratio=:equal,minorgrid=false,yflip=true,xlabel=latexify("x"), ylabel = latexify("z"),yguidefontrotation=-90);

plt22 = plotVar3(M,E,ElementDataIn,VozDataIn,[19401;19401:50:20001|>collect];init_konf=false,np=10,linecolor=:lawngreen,label=:none,aspectratio=:equal,minorgrid=false,yflip=true,xlabel=latexify("x"), ylabel = latexify("z"),yguidefontrotation=-90);

#pltUz2=plot(time_st[1:20:end],M.uz[1,1:20:end],lc=:lawngreen,xlabel=latexify("t"), ylabel = latexify("u_z"),yguidefontrotation=-90,label = latexify("EI=100"));

 savefig(plt11,"out/bm5_1000_conf1.tex")
 savefig(plt12,"out/bm5_1000_conf2.tex")

 savefig(plt21,"out/bm5_100_conf1.tex")
 savefig(plt22,"out/bm5_100_conf2.tex")

 # savefig(pltUz1,"out/bm5_1000_uz.tex")
# savefig(pltW,"out/bm5_soft_W.tex")



#plt3 = plot(plt3,time_st[1:40:end],energija(M,ElementDataIn,E,VozDataIn)[1:40:end];lc=:dodgerblue,label = false,line=:dash)


#const dt::Float64 = 0.25/2
#include("../NonLinBeamRUN.jl")
#plotVar3(M,E,ElementDataIn,VozDataIn,[1;1:2:41|>collect];init_konf=false,p0=plt1,np=10,linecolor=:lawngreen,line=:dash,label=:none,aspectratio=:equal,minorgrid=false,yflip=true)

#=
plt2 = plot(vcat(map(i-> E[1].xInt[1] .+ [0.;cumsum(E[1].L)][i],1:length(E[1].indx))...),M.gamma3[:,end],label = false,linecolor = :black, xticks = 0:2.5:10|>collect,minorgrid=false)
plt3 = plot(time_st[1:40:end],energija(M,ElementDataIn,E,VozDataIn)[1:40:end];lc=:black,label = false)


@assignto :(ElementDataIn) [1] :(t->[repeat([0.],9);8.0*Int(t<=2.5)] ) :(Px)
@assignto :(ElementDataIn) [1] :(t->[repeat([0.],9);-80.0*Int(t<=2.5)]) :(My)
@assignto :(ElementDataIn) [1] :( range(-1,1,length=6) |> collect ) :(div1)
@assignto :(ElementDataIn) [1] :( [4;repeat([2],3);4] ) :(div2)
@assignto :(ElementDataIn) [1] :( repeat([11],5) ) :(nInt)
ElementDataIn[1].Ci=2
const dt::Float64        = 0.125/2

include("../NonLinBeamRUN.jl")
plotVar3(M,E,ElementDataIn,VozDataIn,[1;1:20:201|>collect];init_konf=false, p0=plt1,linecolor=:dodgerblue,line=:dash,label=:none,aspectratio=:equal,minorgrid=false,yflip=true,yrange=(-9,1))
plt2 = plot(plt2,vcat(map(i-> E[1].xInt[1] .+ [0.;cumsum(E[1].L)][i],1:length(E[1].indx))...),M.gamma3[:,end],label = false,linecolor = :dodgerblue,line=:dash, xticks =0:2.0:10|>collect,minorgrid=false)
plt3 = plot(plt3,time_st[1:40:end],energija(M,ElementDataIn,E,VozDataIn)[1:40:end];lc=:dodgerblue,label = false,line=:dash)



@assignto :(ElementDataIn) [1] :(t->[repeat([0.],7);8.0*Int(t<=2.5)] ) :(Px)
@assignto :(ElementDataIn) [1] :(t->[repeat([0.],7);-80.0*Int(t<=2.5)]) :(My)
@assignto :(ElementDataIn) [1] :( range(-1,1,length=5) |> collect ) :(div1)
@assignto :(ElementDataIn) [1] :( [4;repeat([2],2);4] ) :(div2)
@assignto :(ElementDataIn) [1] :( repeat([11],4) ) :(nInt)
const dt::Float64        = 0.125/2
ElementDataIn[1].Ci=2


include("../NonLinBeamRUN.jl")


plotVar3(M,E,ElementDataIn,VozDataIn,[1;1:20:201|>collect];init_konf=false, p0=plt1,linecolor=:lawngreen,line=:dash,label=:none,minorgrid=false,yflip=true,yrange=(-9,1),aspectratio=:equal)

plt2 = plot(plt2,vcat(map(i-> E[1].xInt[1] .+ [0.;cumsum(E[1].L)][i],1:length(E[1].indx))...),M.gamma3[:,end],label = false,linecolor = :lawngreen,line=:dash, xticks =unique([0:2.5:10|>collect;0:2.0:10|>collect]),minorgrid=false)

plt3 = plot(plt3,time_st,energija(M,ElementDataIn,E,VozDataIn);lc=:lawngreen,line=:dash,label = false)



plt1=plot(plt1,[0,0],[0,0],label = latexify("30C^0"), ylabel = latexify("z"),yguidefontrotation=-90,xlabel = latexify("x"),lc=:black)
plt1=plot(plt1,[0,0],[0,0],label = latexify("5C^2"),lc=:dodgerblue,line=:dash)
plt1=plot(plt1,[0,0],[0,0],label = latexify("4C^2"),lc=:lawngreen,line=:dash,aspectratio=:equal)

plt2=plot(plt2,[0,0],[0,0],label = latexify("30C^0"), ylabel = latexify("K_2"),yguidefontrotation=-90,xlabel = latexify("s"),lc=:black)
plt2=plot(plt2,[0,0],[0,0],label = latexify("5C^2"),lc=:dodgerblue,line=:dash)
plt2=plot(plt2,[0,0],[0,0],label = latexify("4C^2"),lc=:lawngreen,line=:dash)

plt3=plot(plt3,[0,0],[0,0],label = latexify("30C^0"), ylabel = latexify("W"),yguidefontrotation=-90,xlabel = latexify("t"),lc=:black)
plt3=plot(plt3,[0,0],[0,0],label = latexify("5C^2"),lc=:dodgerblue,line=:dash)
plt3=plot(plt3,[0,0],[0,0],label = latexify("4C^2"),lc=:lawngreen,line=:dash)

savefig(plt1,"~/0_git/dinamika_nosilcev/out/bm1_konf.tex")
savefig(plt2,"~/0_git/dinamika_nosilcev/out/bm1_K2.tex")
savefig(plt3,"~/0_git/dinamika_nosilcev/out/bm1_W.tex")
=#
