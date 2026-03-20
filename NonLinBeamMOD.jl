module NonLinBeam
	




	# E K S P O R T
	export Beam, BeamDataIn, BeamDataProcess, Node, NodeDataIn,Motion, BeamMotion, MidPoint, TimeElement,
		datainit, dataprocess, iteracija, popravek,
		R, Tan_Res, Mass,
		Integrate, QuadInt, re_gramschmid,InterpolValue,PolyValue,
		trig_re_gramschmid,TrigValue,
		plotbeams, plotmotion,plotVar,plotVar2,plotVar3,plotVar3anim,
		adj_mat,randpermute,node_permute, Cuthill_McKee
	#
	#
	#
	#
	#
	# M O D U L I
	using LinearAlgebra, Plots
    	#
	#
	#
	#
	#
	# E L E M E N T
	abstract type Beam
	end
	#
	# V N E Š E N I   P O D A T K I   N O S I L C A
	@kwdef mutable struct BeamDataIn <:Beam
		v::Vector{Int64} = [1;2] # Vozlišča - krajna
		C::Matrix{Float64} = [1. 0. 0.;0. 1. 0.;0. 0. 1.] #Materialna matrika
		M::Vector{Float64} = [0.0; 0.0] #Vektor [ρA; ρI]
		Ib_geom::Matrix{Float64} = [0.5 0.5; -0.5 0.5] #re_gramschmid(DataIn::Vector{Vector{Float64}})
		Kb::Matrix{Float64} = Array{Float64,2}(undef,(0,2))
		px::Function = t->nothing 
		pz::Function = t->nothing
		my::Function = t->nothing
		Px::Function = t->nothing
		Pz::Function = t->nothing
		My::Function = t->nothing
		div1::Array{Float64} = [-1.;1.]
		div2::Array{Int64} = [4]
		nInt::Array{Int64} = [6]
		Ci::Bool = false	#Zveznost odvodov
	end 
	#
	# P R O C E S I R A N I   P O D A T K I  N O S I L C A
	struct BeamDataProcess <:Beam
		L ::Vector{Float64}
		p0::Vector{Vector{Float64}}
		k0::Vector{Vector{Float64}}
		ksi::Vector{Vector{Float64}} # [0,xInt,L] v začetni parametrizaciji
		P ::Vector{Matrix{Float64}} #
		pb::Vector{Vector{Float64}}
		kb::Vector{Vector{Float64}}
		xInt::Vector{Vector{Float64}} # xInt v naravni parametrizaciji
		wInt::Vector{Vector{Float64}}
		indx::Vector{Vector{Int64}}
		indx_int::Vector{Vector{Int64}}
	end
	#
	#
	#
	#
	#
	# V O Z L I Š Č A
    	abstract type Node 
    	end
	#
	# V N E Š E N I   P O D A T K O I   V O Z L I Š Č A
    	@kwdef mutable struct NodeDataIn <:Node
        # te poračuna algoritem
		x::Float64 = 0.
		z::Float64 = 0.
		i::Int64 = 1
		Supp::Array{Bool} = [1, 1, 1]
		dir::Float64 = 0.
		#mot::Function t-> [0.; 0.; 0.] #Prisiljeno gibanje vozlišča    
	end
	#
	#
	#
	#
	#
	# G I B A N J E
	abstract type Motion
	end
	#
	# R E Z U L T A T I
	mutable struct BeamMotion <:Motion
		ux::Array{Float64,2}
		uz::Array{Float64,2}
		phi::Array{Float64,2}
		vx::Array{Float64,2}
		vz::Array{Float64,2}
		Omg::Array{Float64,2}
		gamma1::Array{Float64,2}
		gamma2::Array{Float64,2}
		gamma3::Array{Float64,2}
	end


	abstract type TimeIntegration
	end

	struct MidPoint <: TimeIntegration
		dt::Float64
		nt::Int64
	end

	struct TimeElement <: TimeIntegration
		dt::Float64
		nodes::Vector{Float64}
		basis::Matrix{Float64}
		DTvalue::Matrix{Float64}
		Tvalue::Matrix{Float64}
		ITvalue::Matrix{Float64}
		xInt::Vector{Float64}
		wInt::Vector{Float64}
	end

























	#Funkcija za rotacijo 3-terice za kot a okrov e3
	function R(a::Float64;n = 0)::Matrix{Float64}
		# n je stopnja odvoda
		return [0. -1. 0.;1. 0. 0.; 0. 0. 0.]^n * [cos(a) -sin(a) 0. ;sin(a) cos(a) 0. ; 0. 0. 1.]
		         #  +     -                               +            -
	end







	#Pripravi osnovne podatke
	function datainit(elem::Matrix{Int64},voz::Matrix{Float64})::Tuple{Int64,Int64,Array{BeamDataIn},Array{NodeDataIn}}
		n_elem::Int64 = length(elem[:,1])
		n_voz::Int64 = length(voz[:,1])


		element_data::Vector{BeamDataIn} = fill(BeamDataIn(),n_elem)
		voz_data::Vector{NodeDataIn} = fill(NodeDataIn(),n_voz)
		
	

		for i =1:n_elem
		    element_data[i] = BeamDataIn(v = elem[i,[1,2]]) 
		end

		for i =1:n_voz
			voz_data[i] = NodeDataIn(x=voz[i,1],z=voz[i,2],i=i)
		end

		return n_elem,n_voz,element_data,voz_data
	end #datainit

	function dataprocess(elem_dat::BeamDataIn,node_dat::Array{NodeDataIn},n_nodes::Int64,n_int_nodes::Int64;mtd::String = "gauss")::Tuple{BeamDataProcess,Int64,Int64}
		
		node1 = [node_dat[1].x node_dat[1].z]
		node2 = [node_dat[2].x node_dat[2].z]
		n_ke = length(elem_dat.div1)-1
		indx = fill(Vector{Int64}([]),n_ke)
		indx_int = fill(Vector{Int64}([]),n_ke)

		if n_ke ==1
			indx[1] = vcat([elem_dat.v[1]],collect(n_nodes+1:n_nodes-2+elem_dat.div2[1]),[elem_dat.v[2]])
			n_nodes = maximum(indx[1])
		elseif n_ke == 2
			indx[1] = vcat([elem_dat.v[1]],(n_nodes+1):(n_nodes+elem_dat.div2[1]-1+Int(elem_dat.Ci)))
			n_nodes = maximum(indx[1])
			if elem_dat.Ci
				indx[2] = vcat([n_nodes-1],collect(n_nodes+1:n_nodes-2+elem_dat.div2[2]),[elem_dat.v[2]],[n_nodes])		
				n_nodes = maximum(indx[2])
			else 
				indx[2] = vcat(collect(n_nodes:n_nodes+elem_dat.div2[2]-2),[elem_dat.v[2]])
				n_nodes = maximum(indx[2])
			end
		else

			indx[1] = vcat([elem_dat.v[1]],(n_nodes+1):(n_nodes-1+elem_dat.div2[1]+Int(elem_dat.Ci)))
			n_nodes = maximum(indx[1]) 
			if elem_dat.Ci
				for i = 2:n_ke-1
					indx[i] = vcat([n_nodes-1],collect(n_nodes+1:n_nodes+elem_dat.div2[i]-1),[n_nodes;n_nodes+elem_dat.div2[i]])
					n_nodes = maximum(indx[i])
				end
				indx[n_ke] = vcat([n_nodes-1],collect(n_nodes+1:n_nodes+elem_dat.div2[n_ke]-2),[elem_dat.v[2]],[n_nodes])
				n_nodes = maximum(indx[n_ke])

			else
				for i = 2:n_ke-1
					indx[i] = collect(n_nodes:n_nodes+elem_dat.div2[i]-1)
					n_nodes = maximum(indx[i])
				end
				indx[n_ke] = vcat(collect(n_nodes:n_nodes-2+elem_dat.div2[n_ke]),[elem_dat.v[2]])		
				n_nodes = maximum(indx[n_ke])
			end
		end

		#koeficienti razoja geometrije
		geom_koeff = vcat(node1,node2,elem_dat.Kb)
		#vektor koeficientov polinoma za x in y
		Geom_Poly = elem_dat.Ib_geom*geom_koeff
		




		#koeficienti razoja geometrije
		#
		Ki = fill(Vector{Float64}([]),n_ke)
		Pi = fill(Vector{Float64}([]),n_ke)
		ksi=copy(Pi)
		xg = fill(Vector{Float64}([]),n_ke)
		wg = fill(Vector{Float64}([]),n_ke)
		Li = fill(Float64(0.),n_ke)
		Kb = fill(Vector{Float64}([]),n_ke)
		Pb = fill(Vector{Float64}([]),n_ke)
	
		#zacetni kot v točkah za integracijo
			#odvajaj f_geom
			#izvrednoti v integracijskih tockah
			#izracunaj kot vektorja
		#zacetna ukrivljenost v točkah za integracijo
			#odvajaj 1x in 2x
			#determinanta v integracijskih tockah iz vektorja
			#norma na 3

		#dolžina krivulje k.e.
			#odvajaj polinom za interpolacijo geom
			#izvrednoti v integracijskih točkah
			#kvadriraj
			#seštej
			#koreni
			#množi z utežm
			#faktor transformacije koordinat (t1,t2)->(-1,1)


		

		#koeficienti razoja geometrije
		for i = 1:n_ke

			xg[i],wg[i] = QuadInt(elem_dat.nInt[i];mtd = mtd)
			indx_int[i] = n_int_nodes+1:n_int_nodes+length(xg[i]) |> collect
			n_int_nodes =n_int_nodes+length(xg[i])



			Li[i] = Integrate(x->norm(PolyValue(x,Geom_Poly;n=1)), [elem_dat.div1[i]], [elem_dat.div1[i+1]],(xg[i],wg[i]))

			ksi[i]=range(elem_dat.div1[i],elem_dat.div1[i+1],length=length(xg[i])+2)

			dksi = [1.0]
			x0 = [-1.;xg[i];1.]

			count=0
			xwint = QuadInt(60)

			while norm(dksi)>10^-10 && count<20
				dksi = map(j-> ((x0[j]+1.0)/2.0*Li[i] - Integrate(x->norm(PolyValue(x,Geom_Poly;n=1)),[elem_dat.div1[i]],[ksi[i][j]],xwint))/norm(PolyValue(ksi[i][j],Geom_Poly;n=1)),eachindex(ksi[i]))	
				
				ksi[i] += dksi
				count+=1
			end
			round.(ksi[i],digits=12)

			

			
			
			D1_vec = map(x->PolyValue(x,Geom_Poly;n=1),ksi[i])
			D2_vec = map(x->PolyValue(x,Geom_Poly;n=2),ksi[i]) 

			#display(D1_vec)
			#display(D2_vec)

			

			Pi[i] = -map(v-> atan(v[2],v[1]),D1_vec)
			Ki[i] = -map((v1,v2)-> (v1[1]*v2[2]-v1[2]*v2[1])/(v1[1]^2+v1[2]^2),D1_vec,D2_vec)
																	#   1



			Kb[i]= Ki[i][[1,end]]
			Pb[i] = Pi[i][[1,end]]

			Ki[i] = Ki[i][2:end-1]
			Pi[i] = Pi[i][2:end-1]
			
			wg[i] = wg[i]*Li[i]/2.0
			xg[i] = (xg[i].+1)*Li[i]/2.0
		end
	


		

		#koeficienti razoja geometrije
		# Interpolacijska baza za vsak končni element
		Ib = map(i -> (!elem_dat.Ci) ? re_gramschmid([collect(range(0.0,Li[i],length = elem_dat.div2[i]))]) : (1<i<n_ke ? re_gramschmid([collect(range(0.0,Li[i],length = elem_dat.div2[i])),[0.0,Li[i]]]) : ( i==1 ? re_gramschmid([collect(range(0.0,Li[i],length = elem_dat.div2[i])),[Li[i]]]) : re_gramschmid([collect(range(0.0,Li[i],length = elem_dat.div2[i])),[0.0]]))),1:n_ke)


		return BeamDataProcess(Li,Pi,Ki,ksi,Ib,Pb,Kb,xg,wg,indx,indx_int), n_nodes,n_int_nodes
	end














    # Funkcija za skico
	function plotbeams(EP::Array{BeamDataProcess},ED::Array{BeamDataIn},VD::Array{NodeDataIn})
		ne = length(EP)
		nv = length(VD)
		nke = map(i-> length(ED[i].div2),1:ne)

		img = plot(;title = "Konstrukcija",aspect_ratio = :equal, xticks = [], yticks = [],yflip = true,xlabel = "x",ylabel = "z")

		scatter!(; xticks = map(i->VD[i].x,1:nv),yticks = map(i->VD[i].z,1:nv))

		for i1 = 1:ne
			node1 = [VD[ED[i1].v[1]].x VD[ED[i1].v[1]].z]
			node2 = [VD[ED[i1].v[2]].x VD[ED[i1].v[2]].z]
			geom_koeff = vcat(node1,node2,ED[i1].Kb)
			Geom_Poly = ED[i1].Ib_geom*geom_koeff

			annotate!(InterpolValue(0.0,geom_koeff[:,1],ED[i1].Ib_geom),InterpolValue(0.0,geom_koeff[:,2],ED[i1].Ib_geom),("  "*string(i1),8,:black,:left))


			#POPRAVI POLOŽAJ cyan
			for i2 = 1:nke[i1]

				xl = range(0.,EP[i1].L[i2],length=ED[i1].div2[i2])
				ksi = range(ED[i1].div1[i2],ED[i1].div1[i2+1],length=ED[i1].div2[i2])

				dksi = [1.0]

				count=0
				xwint = QuadInt(60)

				while norm(dksi)>10^-10 && count<20
					dksi = map(j-> (xl[j] - Integrate(x->norm(PolyValue(x,Geom_Poly;n=1)),[ED[i1].div1[i2]],[ksi[j]],xwint))/norm(PolyValue(ksi[j],Geom_Poly;n=1)),eachindex(ksi))
					ksi = ksi + dksi
					count= count + 1
				end

				round.(ksi,digits=12)

				r0 = hcat(map(j->PolyValue(ksi[j],Geom_Poly),eachindex(ksi))...)


				#xs = range( ED[i1].div1[i2],ED[i1].div1[i2+1],length = ED[i1].div2[i2])
				#scatter!(map(x->InterpolValue(x,geom_koeff[:,1],ED[i1].Ib_geom),xs),map(x->InterpolValue(x,geom_koeff[:,2],ED[i1].Ib_geom), xs); m = :circle, markercolor = :cyan,markersize = 4,labels = :none)
				scatter!(r0[1,:],r0[2,:]; m = :circle, markercolor = :cyan,markersize = 4,labels = :none)

			end



			scatter!(map(x->InterpolValue(x,geom_koeff[:,1],ED[i1].Ib_geom), ED[i1].div1[2:end-1]),map(x->InterpolValue(x,geom_koeff[:,2],ED[i1].Ib_geom), ED[i1].div1[2:end-1]); m = :circle, markercolor = :limegreen,markersize = 6,labels = :none)
			plot!(map(x->InterpolValue(x,geom_koeff[:,1],ED[i1].Ib_geom), -1.0 : 0.05 : 1.0),map(x->InterpolValue(x,geom_koeff[:,2],ED[i1].Ib_geom), -1.0 : 0.05 : 1.0); linecolor = :black,labels = :none)
		end
		

		for i1 = 1:nv
			scatter!([VD[i1].x],[VD[i1].z]; m = :circle, markercolor = :red,markersize = 6,labels = :none,series_annotations = [("  "*string(i1),8,:red,:left)])

		end

		return img	
	end


	function plotmotion(EP::Array{BeamDataProcess},ED::Array{BeamDataIn},VD::Array{NodeDataIn},M::BeamMotion)
		#anim = plot(;aspect_ration =:equal,yflip = true, xlabel = "x", ylabel = "z")
		ne = length(EP)
		nke = map(i-> length(ED[i].div2),1:ne)
		points = 100
		
		R0  = map(i1-> map(i2-> vcat(map(x->InterpolValue(x,vcat([VD[ED[i1].v[1]].x VD[ED[i1].v[1]].z],[VD[ED[i1].v[2]].x VD[ED[i1].v[2]].z],ED[i1].Kb),ED[i1].Ib_geom) , range(ED[i1].div1[i2],ED[i1].div1[i2+1],length=points))'...),1:nke[i1]),1:ne)
		dR = map(i1-> map(i2-> vcat(map(x->InterpolValue(x,vcat([VD[ED[i1].v[1]].x VD[ED[i1].v[1]].z],[VD[ED[i1].v[2]].x VD[ED[i1].v[2]].z],ED[i1].Kb),ED[i1].Ib_geom;n=1) , range(ED[i1].div1[i2],ED[i1].div1[i2+1],length=ED[i1].div2[i1]))'...),1:nke[i1]),1:ne)
		
		
		anim = @animate for i1 = 1:size(M.ux)[2]-1
			for i2 = 1:ne
				for i3 = 1:nke[i2]

					U = [M.ux[:,i1][EP[i2].indx[i3]] M.uz[:,i1][EP[i2].indx[i3]] M.phi[:,i1][EP[i2].indx[i3]]]
					dU = vcat(map(xi-> InterpolValue(xi,U[:,1:2],EP[i2].P[i3]), range(-1,1,length=ED[i2].div2[i3]))'...)					
					
					# POpravi sin cos če je treba.
					# Še začetni kot
					
					p = -map(a -> atan(dR[i2][i3][a,2],dR[i2][i3][a,1]),1:size(dR[i2][i3])[1])
					

					P = [sin.(U[:,3]+p) cos.(U[:,3]+p)].*map(i-> norm(dU[i,:]),size(dU)[1])
					display(P)	
					bi = re_gramschmid([collect(range(-1,1,length=ED[i2].div2[i3])),collect(range(-1,1,length=ED[i2].div2[i3]))])
					A = vcat(map(xi-> InterpolValue(xi,vcat(U[:,1:2],P),bi),range(-1,1,length=points))'...)
					

					R = R0[i2][i3]+A
					plot(R[:,1],R[:,2]; label = false,linecolor = :black,yflip = true)
				end
			end
		end


		gif(anim,"in/gif.gif",fps = 1)				 
		return 	
	end
	#
	# A N I M A C I J A
	function plotVar(M::BeamMotion,EP::BeamDataProcess,ED::BeamDataIn,VD::Array{NodeDataIn},var::Symbol = :uz)

		plot_dat = getfield(M,var)

		anim = @animate for i1 = 1:size(M.ux)[2]
			plot(;yflip = true)
			for i2 = 1:length(EP.P)
				xrange = sum(EP.L[1:i2-1]) .+ (0:0.05:1)*EP.L[i2] |> collect
				xl = (0:0.05:1)*EP.L[i2] |> collect
				plot!(xrange,PolyValue(xl,EP.P[i2]*plot_dat[EP.indx[i2],i1]))
			end
		end
		return anim
	end


	function plotVar2(M::BeamMotion,EP::BeamDataProcess,ED::BeamDataIn,VD::Array{NodeDataIn},var1::Symbol = :ux,var2::Symbol = :uz;ratio::Bool = false,tstep::Int64 = 1,clr::Bool = true)

		plot_dat1 = getfield(M,var1)
		plot_dat2 = getfield(M,var2)
		plot(;yflip = true)
		anim = @animate for i1 = 1:tstep:size(M.ux)[2]

		(clr ? plot(;yflip = true) : nothing)
		for i2 = 1:length(EP.P)
			xrange = sum(EP.L[1:i2-1]) .+ (0:0.05:1)*EP.L[i2] |> collect
			xl = (0:0.05:1)*EP.L[i2] |> collect
			plot!(xrange + PolyValue(xl,EP.P[i2]*plot_dat1[EP.indx[i2],i1]),PolyValue(xl,EP.P[i2]*plot_dat2[EP.indx[i2],i1]))
			plot!(xrange , zeros(length(xrange));linecolor = :black, aspect_ratio = (ratio ? :equal : :auto ),legends = :none)
		end
		end
	return anim
	end


	function plotVar3(M::BeamMotion,EP::Array{BeamDataProcess},ED::Array{BeamDataIn},VD::Array{NodeDataIn},tstep::Int64)

	p = plot(;yflip = true)

	for i1 = 1:length(EP)
		for i2 = 1:length(EP[i1].P)

			xl = (0:0.05:1)*EP[i1].L[i2] |> collect

			node1 = [VD[ED[i1].v[1]].x VD[ED[i1].v[1]].z]
			node2 = [VD[ED[i1].v[2]].x VD[ED[i1].v[2]].z]

			ksi::Array{Float64} = range(ED[i1].div1[i2],ED[i1].div1[i2+1],length=length(xl))
			#Pi = copy(ksi)
			#Ki = copy(ksi)

			geom_koeff = vcat(node1,node2,ED[i1].Kb)
			Geom_Poly = ED[i1].Ib_geom*geom_koeff

			dksi = [1.0]#copy(ksi)*0.0

			count=0
			xwint = QuadInt(60)

			while norm(dksi)>10^-7 && count<20
				dksi = map(j-> (xl[j] - Integrate(x->norm(PolyValue(x,Geom_Poly;n=1)),[ED[i1].div1[i2]],[ksi[j]],xwint))/norm(PolyValue(ksi[j],Geom_Poly;n=1)),eachindex(ksi))
			ksi += dksi
			count+=1
			end

			round.(ksi,digits=12)

			r0 = hcat(map(j->PolyValue(ksi[j],Geom_Poly),eachindex(ksi))...)

			plot!(r0[1,:] + PolyValue(xl,EP[i1].P[i2]*M.ux[EP[i1].indx[i2],tstep]),r0[2,:]+PolyValue(xl,EP[i1].P[i2]*M.uz[EP[i1].indx[i2],tstep]))
			plot!(r0[1,:],r0[2,:];linecolor = :black, aspect_ratio = :equal,legends = :none)
		end
	end

	return p
	end

	function plotVar3anim(M::BeamMotion,EP::Array{BeamDataProcess},ED::Array{BeamDataIn},VD::Array{NodeDataIn},indx_step::Int64)

		anim = @animate for it = 1:indx_step:size(M.ux)[2]
			plotVar3(M,EP,ED,VD,it)
		end
		return anim
	end



	#function MassM()
	#end



	# Na KE
	function Tan_Res(xInt::Vector{Float64},wInt::Vector{Float64},ux::Matrix{Float64},uz::Matrix{Float64},phi::Matrix{Float64},vx::Matrix{Float64},vz::Matrix{Float64},omg::Matrix{Float64},gamma1::Vector{Float64},gamma2::Vector{Float64},gamma3::Vector{Float64},Pval::Matrix{Float64},dPval::Matrix{Float64},Ib::Matrix{Float64},p0::Vector{Float64},k0::Vector{Float64},C::Matrix{Float64},M::Vector{Float64},Fpx::Vector{Float64},Fpz::Vector{Float64} ,Fmy::Vector{Float64},Px::Vector{Float64},Pz::Vector{Float64},My::Vector{Float64},dt::MidPoint,pb::Vector{Float64},kb::Vector{Float64},L::Float64,g::Vector{Float64})

		#dt = tInt.dt

		ux1 = ux[:,1]; uz1 = uz[:,1]; phi1 = phi[:,1]; ux2 = ux[:,2]; uz2 = uz[:,2]; phi2 = phi[:,2];
		vx1 = vx[:,1]; vz1 = vz[:,1]; omg1 = omg[:,1]; vx2 = vx[:,2]; vz2 = vz[:,2]; omg2 = omg[:,2]
		#ux,uz,... so vrednosti v interpolacijskih točkah



		F = fill(Vector{Float64}([0.0;0.0;0.0]),length(ux1))
		dlF = fill(zeros(Float64,(3,3)),(length(ux1),length(ux1)))
		gamma1plus1 = copy(gamma1)
		gamma2plus1 = copy(gamma2)
		gamma3plus1 = copy(gamma3)

		indx2 = CartesianIndex.((1:size(Pval)[2]),(1:size(Pval)[2])')

		D = [0.;1.;0.;;-1.0;0.0;0.0;;0.0;0.0;0.0]
			   #-      +
		Id = [1.0;0.0;0.0;;0.0;1.0;0.0;;0.0;0.0;1.0]
		Mai = [M[1] 0.0 0.0; 0.0 M[1] 0.0; 0.0 0.0 M[2]]

		for i1 = eachindex(xInt)

			e0 = [-1.; 0.; -k0[i1]]

			#   H I T R O S T I
			V1 = map(vi->InterpolValue(xInt[i1],vi,Ib),[vx1,vz1,omg1])
			V2 = map(vi->InterpolValue(xInt[i1],vi,Ib),[vx2,vz2,omg2])
			dV1 = map(vi->InterpolValue(xInt[i1],vi,Ib;n=1),[vx1,vz1,omg1])#*2/L
			dV2 = map(vi->InterpolValue(xInt[i1],vi,Ib;n=1),[vx2,vz2,omg2])#*2/L
			
			V = (V1 +V2) /2.0
			dV= (dV1+dV2)/2.0


			#   P O M I K I
			U1 = map(ui->InterpolValue(xInt[i1],ui,Ib),[ux1,uz1,phi1])
			dU1 = map(ui->InterpolValue(xInt[i1],ui,Ib;n=1), [ux1,uz1,phi1])#*2/L

			U = U1+V*dt/2.0
			dU = dU1 +dV*dt/2.0

			#	D E F O R M A C I J E   F I K S N A   B A Z A
			dq1 = [cos(p0[i1]);-sin(p0[i1]);k0[i1]]+dU1
			dq = dq1 + dV*dt/2.0


			Rm = R(U[3]+p0[i1])

			#	D E F O R M A C I J E
			E1  = [gamma1[i1],gamma2[i1],gamma3[i1]]
			E_e = Rm*dq
			E   = E_e+e0
			E2  = E1 +dt*(Rm*dV + V[3]*D*E_e)
			gamma1plus1[i1] = E2[1]
			gamma2plus1[i1] = E2[2]
			gamma3plus1[i1] = E2[3]

			#	N O T R A N J E   S I L E
			N = C*(E1+E2)/2

			#	N O T R A N J E   S I L E   F I K S N A
			Re = Rm'*N

			#	V E K T O R S K I   P R O D U K T
			X = -[0.0;0.0;1.0]*dot(N,D,E_e)
			   #-

			# 	L I N E A R I Z A C I J A
			LR = (Id+dt/2*V[3]*D)



			#dlE2 = C*dt*D*(e0)
			dlRe= (dt/2*Rm'*(-D*N + C*D*(dt/2*Rm*dV + LR*E_e)), dt/2*Rm'*C*LR*Rm)
							#-   +
			dlX = ([0.0;0.0;dt/2.0]*(E_e'*D*(-D*N + C*D*(dt/2*Rm*dV + LR*E_e)))',
											#-
					[0.0;0.0;dt/2.0]*(D*Re + (E_e'*D*C*LR*Rm)' )')
					#                  +

			#OBTEŽBA  !!!POPRAVI INTERPOLACIJO!!!
			p = map(pj -> PolyValue(xInt[i1],[1. 0.;-1.0/L  1.0/L]*reshape(pj,(2))),[Fpx,Fpz,-Fmy])
						#                                                                    +

			#POSPEŠEK

			tv = Mai*((V2-V1) + dt*[g;0.])


			F .+= map(i2 ->  dPval[i1,i2]*(-Re*dt) + Pval[i1,i2]*(p*dt - tv + X*dt),1:length(Ib[:,1])) * wInt[i1]


			dlF .+= map( ij ->
				-dPval[i1,ij[1]]*(
					Pval[i1,ij[2]]*[0.0;0.0;0.0;;0.0;0.0;0.0;;dlRe[1]*dt]
					+dPval[i1,ij[2]]*dlRe[2]*dt)
				+Pval[i1,ij[1]]*(
					Pval[i1,ij[2]]*([0.0;0.0;0.0;;0.0;0.0;0.0;;dlX[1]*dt])
					+dPval[i1,ij[2]]*dlX[2]*dt),
				indx2)*wInt[i1]

		end
		P = [Px[1];Pz[1];-My[1]]
		#    			 +
		F .+= map(i2 -> dt*P*PolyValue(0.0,Ib[:,i2]) ,1:length(Ib[:,1]))
		P = [Px[2];Pz[2];-My[2]]
		#    			 +
		F .+= map(i2 -> dt*P*PolyValue(L,Ib[:,i2]) ,1:length(Ib[:,1]))

		return dlF, F, gamma1plus1,gamma2plus1,gamma3plus1
	end

	# Na KE
	function Tan_Res(xInt::Vector{Float64}, wInt::Vector{Float64}, ux::Matrix{Float64}, uz::Matrix{Float64}, phi::Matrix{Float64}, vx::Matrix{Float64}, vz::Matrix{Float64}, omg::Matrix{Float64}, gamma1::Vector{Float64}, gamma2::Vector{Float64}, gamma3::Vector{Float64}, Pval::Matrix{Float64}, dPval::Matrix{Float64}, Ib::Matrix{Float64}, p0::Vector{Float64}, k0::Vector{Float64}, C::Matrix{Float64}, M::Vector{Float64}, Fpx::Matrix{Float64}, Fpz::Matrix{Float64} , Fmy::Matrix{Float64}, Px::Matrix{Float64}, Pz::Matrix{Float64}, My::Matrix{Float64}, tInt::TimeElement, pb::Vector{Float64}, kb::Vector{Float64}, L::Float64, g::Vector{Float64})

	#dt = tInt.dt
	gamma1plus1 = copy(gamma1)
	gamma2plus1 = copy(gamma2)
	gamma3plus1 = copy(gamma3)

	F = fill(Vector{Float64}([0.0;0.0;0.0]),size(vx))
	dlF = fill(zeros(Float64,(3,3)),(length(vx),length(vx)))


	indxF = CartesianIndex.((1:size(Pval)[2]),(1:size(tInt.Tvalue)[2])')
	indx2 = CartesianIndex.((1:size(Pval)[2]),(1:size(Pval)[2])')

	D = [0.;1.;0.;;-1.0;0.0;0.0;;0.0;0.0;0.0]
	#-      +
	Id = [1.0;0.0;0.0;;0.0;1.0;0.0;;0.0;0.0;1.0]
	Mai = [M[1] 0.0 0.0; 0.0 M[1] 0.0; 0.0 0.0 M[2]]

	#=
	struct TimeElement <: TimeIntegration
	dt::Float64
	nodes::Vector{Float64}
	basis::Matrix{Float64}
	DTvalue::Matrix{Float64}
	Tvalue::Matrix{Float64}
	ITvalue::Matrix{Float64}
	xInt::Vector{Float64}
	wInt::Vector{Float64}
	end
	=#

	for i1 = eachindex(xInt)
		e0 = [-1.; 0.; -k0[i1]]

		#   H I T R O S T I
		V = vcat(map(vi->PolyValue(xInt[i1],Ib*vi),[vx1,vz1,omg1])...)
		#  V = vcat(map(vi-> Pval[i1,:]*vi ,[vx1,vz1,omg1])...)
		dV = vcat(map(vi->PolyValue(xInt[i1],Ib*vi;n=1),[vx1,vz1,omg1])...)#*2/L
		# dV = vcat(map(vi-> dPval[i1,:]*vi ,[vx1,vz1,omg1])...)

		#   P O M I K I
		Un = map(ui->InterpolValue(xInt[i1],ui,Ib),[ux1[:,1],uz1[:,1],phi1[:,1]])
		#  Un = map(ui-> Pval[i1,:]*ui ,[ux1[:,1],uz1[:,1],phi1[:,1]])
		dUn = map(ui->InterpolValue(xInt[i1],ui,Ib;n=1), [ux1[:,1],uz1[:,1],phi1[:,1]])#*2/L
		# dUn = map(ui-> dPval[i1,:]*ui ,[ux1[:,1],uz1[:,1],phi1[:,1]])

		dqn = [cos(p0[i1]);-sin(p0[i1]);k0[i1]]+dUn

		R0 = R(p0[i1])

		for t1 = eachindex(tInt.wInt)

			U = Un + V*tInt.ITvalue[t1,:]'
			dU = dUn + dV*tInt.ITvalue[t1,:]'

		#	D E F O R M A C I J E   F I K S N A   B A Z A

			dq = dqn + dV*tInt.ITvalue[t1,:]'

			Rst = R(U[3])*R0
			E_e = Rst*dq

			Est = E_e + e0
			Nst = C*Est
			Rst = Rst'*Nst
			Xst = -[0.;0.;1.]*Nst'*D*E_e

			tV = Mai*(V*tInt.DTvalue[t1,:]' + [g;0.])




			#dlE2 = C*dt*D*(e0)

			#-   +
			dlX = ([0.0;0.0;dt/2.0]*(E_e'*D*(-D*N + C*D*(dt/2*Rm*dV + LR*E_e)))',
				#-
				[0.0;0.0;dt/2.0]*(D*Re + (E_e'*D*C*LR*Rm)' )')
			#                  +


			p = vcat(map(pj -> PolyValue(xInt[i1],[1. 0.;-1.0/L  1.0/L]*pj),[Fpx[:,t1],Fpz[:,t1],-Fmy[:,t1]])...)
			#                                                                                    +


			F .+= map(it2 ->  tInt.Tvalue[t1,it2[2]]*(-(dPval[i1,it2[1]]*Id + Pval[i1,it2[1]]*[0.;0.;1.0]*dq'*D)*Rst + Pval[i1,i2]*(p - tV ))   ,indxF) * wInt[i1]*tInt.wInt[t1]


			dlF .+= map( ij ->
					-dPval[i1,ij[1]]*(
						Pval[i1,ij[2]]*[0.0;0.0;0.0;;0.0;0.0;0.0;;dlRe[1]*dt]
						+dPval[i1,ij[2]]*dlRe[2]*dt)
					+Pval[i1,ij[1]]*(
						Pval[i1,ij[2]]*([0.0;0.0;0.0;;0.0;0.0;0.0;;dlX[1]*dt])
						+dPval[i1,ij[2]]*dlX[2]*dt),
					indx2)*wInt[i1]

		end
	end
	P = vcat([Px[1,:];Pz[1,:];-My[1,:]]...)
	#    			 +
	for t1 = eachindex(tInt.wInt)
		F .+= map(it2 -> P[:,t1]*PolyValue(0.0,Ib[:,it2[1]])*tInt.Tvalue[t1,it2[2]] ,indxF)*tInt.wInt[t1]
	end
	P = vcat([Px[2,:];Pz[2,:];-My[2,:]]...)
	#    			 +
	for t1 = eachindex(tInt.wInt)
		F .+= map(it2 -> P[:,t1]*PolyValue(L,Ib[:,it2[1]])*tInt.Tvalue[t1,it2[2]] ,indxF)*tInt.wInt[t1]
	end

	return dlF, F, gamma1plus1,gamma2plus1,gamma3plus1

	end

	function iteracija(ux::Vector{Float64}, uz::Vector{Float64}, phi::Vector{Float64}, vx::Vector{Float64}, vz::Vector{Float64}, Omg::Vector{Float64}, Dvx::Vector{Float64}, Dvz::Vector{Float64}, DvOmg::Vector{Float64}, tInt::MidPoint)
		dt = tInt.dt

		vx += Dvx*2
		vz += Dvz*2
		Omg += DvOmg*2

		ux += Dvx*tInt.dt
		uz += Dvz*tInt.dt
		phi += DvOmg*tInt.dt

	return vx,vz,Omg,ux,uz,phi
	end

	function prediktor(ux::Matrix{Float64}, uz::Matrix{Float64}, phi::Matrix{Float64}, vx::Matrix{Float64}, vz::Matrix{Float64}, Omg::Matrix{Float64}, tInt::MidPoint)
		return vx[:,1],vz[:,1],Omg[:,1],ux[:,1]+vx[:,1]*dt, uz[:,1]+vz[:,1]*dt,phi[:,1]+Omg[:,1]*dt
	end


	function Mass(xInt::Vector{Float64}, wInt::Vector{Float64}, Pval::Matrix{Float64}, dPval::Matrix{Float64}, M::Vector{Float64},tInt::MidPoint)
		Mass = fill(zeros(Float64,(3,3)),(size(Pval)[2],size(Pval)[2]))
		indx2 = CartesianIndex.((1:size(Pval)[2]),(1:size(Pval)[2])')
		Mai = diagm(M[[1,1,2]])
		for i1 = eachindex(xInt)
			Mass .+= map( ij -> Pval[i1,ij[1]]*Pval[i1,ij[2]]*(-2.0*Mai), indx2)*wInt[i1]
		end
		return Mass
	end


#	 A L G O R I T M I:
#
#
#		- I N T E R P O L A C I J S K A   	B A Z A
#
#		- N U M E R I Č N A   			I N T E G R A C I J A
#
#		- I Z V R E D N O T E N J E   		I N T E R P O L A C I J E
#
#		- P E R M U T A C I J E 		V O Z L I Š Č			G R A F A
	
	
	
	
	
	# V E K T O R   V R E D N O S T I   B A Z E
	function DotP(a::Array{Float64},DataIn::Vector{Vector{Float64}})
		m = length(DataIn)
		n = sum(length.(DataIn))-1
	#Diferencialni operator
		Dp = diagm(1=>1:n)

	#Seznam funkcijskih vrednosti polinoma glede na DataIn
		b = vcat(map( (A,i) -> (((Dp^i)*a)'* (A'.^(0:n)) )' ,DataIn,0:(m-1))...)
		return b
	end
	#
	# S K A L A R N I   P R O D U K T
	function DotP(a::Array{Float64},b::Array{Float64},DataIn::Vector{Vector{Float64}})
		#Skalarni produkt vektorjev iz DotP(a,A)
		return DotP(a::Array{Float64},DataIn::Vector{Vector{Float64}})'*DotP(b::Array{Float64},DataIn::Vector{Vector{Float64}})
	end
	#
	# O R T O G O N A L I Z A C I J A
	function gramschmid(DataIn::Vector{Vector{Float64}}; B0::Union{Array{Float64},Int64} = 1)
        	n::Int64 = sum(length.(DataIn))-1
        	
		#Standardna baza        
       	 	Id = Matrix{Float64}(I,n+1,n+1);

		#Re-ortogoanalizacija
        	if typeof(B0) == Int64 
		    B0::Array{Float64} = copy(Id)
		end

		#Diferencialni operator
		Dp::Array{Float64} = diagm(1=>Float64.(1:n))

		#Skalarni produkti standardne baze + normiranje
		ei = B0./sqrt.(DotP(B0,B0,DataIn)[CartesianIndex.(1:n+1,1:n+1)])

		#Ortogonalizacija
		bi = copy(ei);
		for i1 = 2:n+1
		    c1 =ei[:,i1] - sum(DotP(ei[:,i1],bi[:,1:i1-1],DataIn) .*bi[:,1:i1-1],dims=2)[:,1]
		    bi[:,i1] = c1/sqrt(DotP(c1,c1,DataIn)[1])
		end
		
		#Interpolacija
		Ib = hcat(map(i-> round.(sum(DotP(bi,DataIn)[i,:]'.*bi,dims=2),digits = 13) ,1:n+1)...)

		return Ib,n,Dp,bi
	end
	#
	# D V O J N A   O R T O G O N A L I Z A C I J A
	function re_gramschmid(DataIn::Vector{Vector{Float64}})
		bi = gramschmid(DataIn)[4]
		bi = gramschmid(DataIn;B0 = bi)[1]
		if all(bi[end,:].<10^(-13))
			error("Interpolacijska baza ne obstaja")
		end
		return bi
	end
	#
	#
	#
	#
	#
	# I N T E G R A C I J A   Z   V O Z L I Š Č I   I N   U T E Ž M I
	function QuadInt(n::Int64;mtd::String = "gauss")

		if mtd == "gauss"
			# 2n-1
			b = map( i-> (i+1)/(((2*i+1)*(2*i+3))^0.5 ),0:n-2)

			E = eigen(diagm(1=>b,-1=>b))

			xg = E.values
			wg = E.vectors[1,:].^2*2

			return xg,wg

		elseif mtd == "lobatto"
			# 2n-3
			# 2(n+1)-3 = 2n-1
			n += 1
			if n == 2
				xg = [-1.0,1.0]
				wg = [1.0,1.0]
			else
				n = n-2
				b = map( i1->  i1/(2*i1+3)*(i1+2)/(2*i1+1), 1:n-1).^0.5
				E = eigen(diagm(1=>b,-1=>b))
				xg = E.values

				# Poiskusi optimizirat
				P = [[1.0],[0.0,1.0]]
				bI = map( i-> (i+1)^2.0/(((2*i+1)*(2*i+3)) ),0:n-1)
				for i = 1:n
					P = [P[2],vcat([0],P[2])-bI[i]*vcat(P[1],[0.0;0.0])]
				end
				wg = 2.0./((n+2)*(n+1)*map(xi->sum(P[2].*(xi.^(0:(length(P[2])-1))))^2,xg))
				wg = vcat(2.0/((n+1)*(n+2)),wg,2/((n+1)*(n+2)))
				wg[2:end-1] = (2-2*wg[1])/sum(wg[2:end-1])*wg[2:end-1]
				xg = vcat(-1.0,xg,1.0)
				return xg,wg
			end
		end
		return xg,wg
	end
	#
	# I N T E G R A C I J A
	function Integrate(f::Function,a0::Float64,a1::Float64;n::Int64 = 30,mtd::String = "gauss")
		xw = QuadInt(n;mtd = mtd)
		I = Integrate(f::Function,a0::Float64,a1::Float64,xw::Tuple{Array{Float64,1},Array{Float64,1}})
		return I
	end
	#
	#
	function Integrate(f::Function,a0::Float64,a1::Float64,xw::Tuple{Array{Float64,1},Array{Float64,1}})
		I = xw[2]'*f.(xw[1]*(a1-a0)/2 .+(a1+a0)/2)*(a1-a0)/2
		return I
	end
	#
	#
	function Integrate(f::Function,a0::Vector{Float64},a1::Vector{Float64};n=30,mtd = "gauss")
		xw = QuadInt(n)
		xtform = map(i -> (xw[1]*(a1[i]-a0[i])/2 .+ (a1[i]+a0[i])/2) , eachindex(a0))
		wtform = map(i -> (xw[2]*(a1[i]-a0[i])/2),eachindex(a0))
		I = Integrate(f::Function,a0::Vector{Float64},a1::Vector{Float64},repeat([xw],length(a0))...)
		return I
	end
	#
	#
	function Integrate(f::Function,xw::Tuple{Array{Float64,1},Array{Float64,1}}...)

		xtform = map(i->xw[i][1],eachindex(xw))
		wtform = map(i->xw[i][2],eachindex(xw))
		indx = Tuple.(CartesianIndex.(map(i-> reshape(1:length(xtform[i]),vcat(ones(Int64,i-1),[length(xtform[i])])...), eachindex(xw))...))
		I = map(i -> prod(getindex.(wtform,i))*f(getindex.(xtform,i)...),indx) |> sum
		I = round.(I,digits = 12)
		return I
	end

	function Integrate(f::Function,a0::Vector{Float64},a1::Vector{Float64},xw::Tuple{Array{Float64,1},Array{Float64,1}}...)
		xtform = map(i -> (xw[i][1]*(a1[i]-a0[i])/2 .+ (a1[i]+a0[i])/2) , eachindex(a0))
		wtform = map(i -> (xw[i][2]*(a1[i]-a0[i])/2),eachindex(a0))

		indx = Tuple.(CartesianIndex.(map(i-> reshape(1:length(xtform[i]),vcat(ones(Int64,i-1),[length(xtform[i])])...), eachindex(a0))...))
		I = map(i -> prod(getindex.(wtform,i))*f(getindex.(xtform,i)...),indx) |> sum
		I = round.(I,digits = 12)
		return I
	end
	#
	#
	function Integrate(P::Array{Float64,2})
		return vcat(zeros(1,size(P)[2]),P.*(1.0./(1:size(P)[1])))
	end
	#
	#
	#
	#
	# I Z V R E D N O T E N J E   I N T E R P O L A C I J E
	function InterpolValue(x::Float64,Kb::Array{Float64},Ib::Matrix{Float64};n::Int64=0)
		f = PolyValue(x::Float64,Ib*Kb::Array{Float64};n = n)
		return f
	end
	#
	# V R E D N O S T   P O L I N O M A 
	function PolyValue(x::Float64,Ki::Array{Float64,1};n::Int64 = 0)
		dP =reverse( (diagm(1 => 1. : size(Ki)[1]-1.)^n)[1:end-n,:]*Ki)
		f = 0.0
		for i in eachindex(dP)
			f *= x
			f += dP[i]
		end
		#f = dP'*x.^(0:size(Ki)[1]-1)
		return f
	end
	function PolyValue(x::Float64,Ki::Array{Float64,2};n::Int64 = 0)
		dP =reverse( (diagm(1 => 1. : size(Ki)[1]-1.)^n)[1:end-n,:]*Ki,dims=1)
		f = zeros(Float64,size(dP)[2])
		for i in 1:size(dP)[1]
			f *= x
			f += dP[i,:]
		end
		#f = dP'*x.^(0:size(Ki)[1]-1)
		return f
	end
	function PolyValue(x::Array{Float64},Ki::Array{Float64,1};n::Int64 = 0)
		dP =reverse( (diagm(1 => 1. : size(Ki)[1]-1.)^n)[1:end-n,:]*Ki)
		f = copy(x)*0
		for i1 in eachindex(f)
			for i2 in eachindex(dP)
				f[i1] *= x[i1]
				f[i1] += dP[i2]
			end
		end
		#f = dP'*x.^(0:size(Ki)[1]-1)
		return f
	end

	function PolyValue(x::Vector{Float64},Ki::Matrix{Float64};n::Int64 = 0)
	dP =reverse( (diagm(1 => 1. : size(Ki)[1]-1.)^n)*Ki,dims=1)
	f = zeros(length(x),size(dP)[2])
	for i1 in eachindex(x)
		for i2 in 1:size(dP)[2]
			f[i1,:] *= x[i1]
			f[i1,:] += dP[i2,:]
		end
	end
	#f = dP'*x.^(0:size(Ki)[1]-1)
	return f
	end







	function trig_DotP(a::Array{Float64},b::Array{Float64},DataIn::Vector{Vector{Float64}})
	#Skalarni produkt vektorjev iz DotP(a,A)
	return TrigValue(DataIn,a)'*TrigValue(DataIn,b)
	end

	#
	# O R T O G O N A L I Z A C I J A
	function trig_re_gramschmid(DataIn::Vector{Vector{Float64}}; B0::Array{Float64} = zeros(Float64,2,2), rep::Int64 = 0)
	na = sum(length.(DataIn))
	n::Int64 = Int(ceil((sum(length.(DataIn)))/2))*2
	nb = size(B0)[2]

	#Standardna baza
	Id = Matrix{Float64}(I,n,n);

	#Re-ortogoanalizacija
	if B0 == zeros(Float64,2,2)
	B0::Array{Float64} = copy(Id)
	nb = n
		end

		#Diferencialni operator
		#Dp::Array{Float64} = diagm(1=>Float64.(1:n))

		#Skalarni produkti standardne baze + normiranje
		ei = B0./sqrt.(trig_DotP(B0,B0,DataIn)[CartesianIndex.(1:nb,1:nb)])'
		bi = copy(ei);

		for i1 = 2:nb

		c1 =ei[:,i1] - sum(trig_DotP(ei[:,i1],bi[:,1:i1-1],DataIn) .*bi[:,1:i1-1],dims=2)[:,1]

		if (a=sqrt(trig_DotP(c1,c1,DataIn)[1]))<10^-10 || any(isnan.(bi[:,i1]))
		bi[:,i1] .= 0.0
		else
		bi[:,i1] = c1/a
		end
		end

		#Interpolacija
		Ib = hcat(map(i-> round.(sum(TrigValue(DataIn,bi)[i,:]'.*bi,dims=2),digits = 13) ,1:na)...)

		#Ib = bi*PolyValue(DataIn::Vector{Vector{Float64}},bi::Matrix{Float64})

		if rep > 0
			Ib,bi = trig_re_gramschmid(DataIn;B0 = bi, rep = rep-1)
		end

		return Ib,bi
		end

	function TrigValue(x::Float64,ai::Matrix{Float64};n::Int64 = 0)
		p = zeros(Float64,1,size(ai)[2])
		if n>0
		da = hcat([[i;0.] for i = 1:Int(size(ai)[1]/2)]...)
			da = reshape(da,length(da))
			Dp = diagm(1=>da,-1=>-da)[1:end-1,1:end-1]
			ai = Dp^n * ai
		end
		for i in 1:2:size(ai)[1]
		p +=  [sin(i*x);; cos(i*x)] * ai[[i,i+1],:]
		end
		return p
	end
	function TrigValue(x::Float64,ai::Vector{Float64};n::Int64 = 0)
		p = 0.
		if n>0
		da = hcat([[i;0.] for i = 1:Int(size(ai)[1]/2)]...)
			da = reshape(da,length(da))
			Dp = diagm(1=>da,-1=>-da)[1:end-1,1:end-1]
			ai = Dp^n * ai
		end
		for i in 1:2:size(ai)[1]
		p +=  ([sin(i*x);; cos(i*x)] * ai[[i,i+1]] )[1]
		end
		return p
	end
	function TrigValue(x::Vector{Float64},ai::Matrix{Float64};n::Int64 = 0)
		return vcat(map(i-> TrigValue(i::Float64,ai::Matrix{Float64};n=n),x)...)
	end
	function TrigValue(x::Vector{Float64},ai::Vector{Float64};n::Int64 = 0)
		return vcat(map(i-> TrigValue(i::Float64,ai::Vector{Float64};n=n),x)...)
	end
	function TrigValue(x::Vector{Vector{Float64}},ai::Matrix{Float64})
		return vcat(map(i-> TrigValue(x[i]::Vector{Float64},ai::Matrix{Float64};n=i-1),eachindex(x))...)
	end
	function TrigValue(x::Vector{Vector{Float64}},ai::Vector{Float64})
		return vcat(map(i-> TrigValue(x[i]::Vector{Float64},ai::Vector{Float64};n=i-1),eachindex(x))...)
	end
	#
	#
	#
	#
	# S O S E D N O S T N A   M A T R I K A
    	function adj_mat(conn::Matrix{Int64},n::Int64)::Matrix{Int64}
		m = size(conn)[1]

		A = zeros(Int64,(n,n))
			
		for i=1:m	
			A[conn[i,1],conn[i,2]] = A[conn[i,2],conn[i,1]] = 1
		end
		return A
	end
	#
	# N A K L J U Č N A   P E R M U T A C I J A
	function randpermute(n::Int64)
		indx  = collect(1:n)
		P = Matrix(I,(n,n))
		P2 = copy(P)

		for i= eachindex(indx)
			c = rand(indx)
			P[:,i] = P2[:,c]
			popat!(indx,findfirst(indx .== c))
		end
		return P
	end
	#
	# M A N J Š A N J E   D I A G O N A L N E G A   P A S U
	function node_permute(conn::Matrix{Int64},nodes::Matrix{Float64};rep = 10)
		m = size(nodes)[1]
		n = size(conn)[1]

		A = adj_mat(conn,m)
		P = Matrix(I,(m,m))

		indx = collect(1:m)
		n_diag = maximum(map(ij-> A[ij] == 0 ? 0 : abs(ij[1]-ij[2]),CartesianIndex.((1:m)',1:m)))

		c = 1
	  	while c < rep
			P2 = randpermute(m)
			A2 = adj_mat((P2*indx)[conn],m)
			n_diag2 = maximum(map(ij -> A2[ij] == 0 ? 0 : abs(ij[1]-ij[2]), CartesianIndex.((1:m)',1:m)))

			if n_diag2 < n_diag
				A = A2
				n_diag = n_diag2
				P = P2
				c = 1
			else
				c += 1
			end
			if n_diag == 1
				break
			end
		end

		nodes = P'*nodes
		conn = (P*indx)[conn]

		return conn,nodes
	end

	function Cuthill_McKee(conn::Array{Int64,2},n::Int64)
		indx_permute = collect(1:n)
		indx = copy(indx_permute)

	#=
		Adj = adj_mat(conn,n)
		deg = sum(Adj;dims = 2)[:,1]

		node_deg = getindex.(findall(deg .== collect(1:n)'),[1 2]) # Min degree > 0

		#indx_permute[1] = node_deg[1,1]

		#Re indeksirane
		R = [node_deg[1][1]]
		re = 1
		#Prvi sosedje
		Q = findall(R .== findall(Adj[:,R].==1))
		n_node = 2

		while

			Q = Q[getindex.(findall(Q[:,2] .== sort(Q[:,w])'),1),:]
			#indx_permute[n_node.+1:size(Q)[1]] = Q[:,1]
			push!(R,Q[:,1]	)

		end

	=#


		return indx_permute
	end
		

end # module
