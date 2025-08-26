module NonLinBeam
    
	export Beam, BeamDataIn, BeamDataProcess, Node, NodeDataIn,Motion, BeamMotion,
		datainit, readdata, dataprocess,
		R, GaussInt, re_gramschmid, Tan_Res, VarsAtX,InterpolValue,PolyValue


    using LinearAlgebra


    # Struktura podatkov za nosilec
    abstract type Beam end

    @kwdef mutable struct BeamDataIn <:Beam
        v::Vector{Int64} = [1;2] # Vozlišča - krajna
        
	C::Matrix{Float64} = [1. 0. 0.;0. 1. 0.;0. 0. 1.] #Materialna matrika
        M::Matrix{Float64} = [1.0 1.0] #Vektor [ρA; ρI]

	Ib_geom::Matrix{Float64} = [0.5 0.5; -0.5 0.5] #re_gramschmid(DataIn::Vector{Vector{Float64}})
	Kb::Matrix{Float64} = Array{Float64,2}(undef,(0,2)) 
	
        px::Function = t->[0. 0.] 
	pz::Function = t->[0. 0.]
        my::Function = t->[0. 0.]

        div1::Array{Float64} = [-1.,1.]
	#Naj bo vedno med -1 in 1

        div2::Array{Int64} = [4]
	#Številow vozlišč v sekundarni delitvi. Dve sta robni. Int za vsak element.

        nInt::Array{Int64} = [20]
	#Število integracijskih točk v posameznem div1. Smiselno je približno 3*div2. Ponovno Int za vsak element
	

	Ci::Bool = false
	#Ib::Array{Matrix{Float64}} = ... 
	#Zveznost odvodov
	#Bi shranil tudi baze
    end 

	#Li,Pi,Ki,Ib,xi,wi
    struct BeamDataProcess <:Beam
       # C::Array{Float64} #
       # M::Array{Float64} #
        L::Vector{Float64} 
	p0::Vector{Vector{Float64}}
	k0::Vector{Vector{Float64}}
        P::Vector{Matrix{Float64}} #
	pb::Vector{Vector{Float64}}
	kb::Vector{Vector{Float64}}
	# dP::Array{Array{Array{Float64}}} #
       # A1::Array{Array{Array{Float64}}} #
       # A2::Array{Array{Array{Float64}}} #
       # A3::Array{Array{Array{Float64}}} #
       # A4::Array{Array{Array{Float64}}} #
       # A5::Array{Array{Array{Float64}}} #
       # A6::Array{Array{Array{Float64}}} #
        xInt::Vector{Vector{Float64}}
        wInt::Vector{Vector{Float64}} # Jih dejansko rabim s sabo če so v A1,...? Lahko jih poračunam samo lokalno.
       # v::Array{Array{Int64}} #
	indx::Vector{Vector{Int64}}
    end
    


    abstract type Node 
    end

    @kwdef mutable struct NodeDataIn <:Node
        # te poračuna algoritem
	x::Float64 = 0.
        z::Float64 = 0.
        i::Int64 = 1

	# prilagodi tako, da lahko nastavis drsno pod kotom.
	# sprostitev v smeri vektorja [x,y] pomeni vezno enačbo
	# y*ux-x*uy = 0 [ux,uy] je pravokoten na [y,-x]
        Supp::Array{Bool} = [1, 1, 1]
	dir::Float64 = 0.
	# lahko nastaviš številko za rotacijo globalnih koordinat v kateri je ta sprostitev jasna in je potem supp podan za ta rotiran koordinatni sistem. smiselno je le za 0-pi/2
	# dokler ne urediš vsega naj bo vedno 0.
    end



    abstract type Motion end

    mutable struct BeamMotion <:Motion
        ux::Array{Float64}
        uz::Array{Float64}
        phi::Array{Float64}
        vx::Array{Float64}
        vz::Array{Float64}
        Omg::Array{Float64}
    end





    #Funkcija za rotacijo 3-terice za kot a okrov e3
    function R(a;n = 0)::Matrix{Float64}
	    # n je stopnja odvoda
	    return [0. 1. 0.;-1. 0. 0.; 0. 0. 0.]^n * [cos(a) sin(a) 0. ; -sin(a) cos(a) 0. ; 0. 0. 1.]
    end




    #Funkcija za določitev uteži wi in koordinat xi za gaussovo integracijo
    function GaussInt(n::Int64)
        if n<2
            error("Število integracijskih točn naj bo večje od 2. (Vnešena $n)")
        end

        b = map( i-> (i+1)/(((2*i+1)*(2*i+3))^0.5 ),0:n-2)
        K = diagm(1=>b,-1=>b)
        E = eigen(K)

        wg = E.vectors[1,:].^2*2
        xg = E.values

        return xg,wg
    end
    function GaussInt(n::Array{Int64})
        if any(n.<2)
            error("Število integracijskih točn naj bo večje od 2. (Vnešena $n)")
        end

        b = map.( i-> (i+1)/(((2*i+1)*(2*i+3))^0.5 ), collect.(range.(0,n.-2)))
        K = map(k ->diagm(1=>k,-1=>k),b)
        E = eigen.(K)

        wg = map(k->E[k].vectors[1,:].^2*2,1:length(n))
        xg = map(k->E[k].values,1:length(n))

        return xg,wg
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

	function dataprocess(elem_dat::BeamDataIn,node_dat::Array{NodeDataIn},n_nodes::Int64)::Tuple{BeamDataProcess,Int64}
		node1 = [node_dat[1].x node_dat[1].z]
		node2 = [node_dat[2].x node_dat[2].z]
		n_ke = length(elem_dat.div1)-1
		indx = fill(Vector{Int64}([]),n_ke)
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
					indx[i] = vcat([n_nodes-1],collect(n_nodes+1:n_nodes+elem_dat.div[i]-1),[n_nodes;n_nodes+elem_dat.div2[i]])
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



		f_geom = sum(map( i -> elem_dat.Ib_geom[:,i].*geom_koeff[i,:]',1:size(elem_dat.Ib_geom)[1]))
		#diff operator za geometrijo



		
		Df_geom = diagm(1=>1:length(geom_koeff[:,1])-1)



		

		#koeficienti razoja geometrije
		Ki = fill(Vector{Float64}([]),n_ke)
		Pi = fill(Vector{Float64}([]),n_ke)
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
			#množi z utežmi
			#faktor transformacije koordinat (t1,t2)->(-1,1)


		

		#koeficienti razoja geometrije
		for i = 1:n_ke

			xg[i],wg[i] = GaussInt(elem_dat.nInt[i])

			x_trans = vcat([elem_dat.div1[i] ; elem_dat.div1[i+1]],xg[i]/2. *(elem_dat.div1[i+1]-elem_dat.div1[i]).+(elem_dat.div1[i+1]+elem_dat.div1[i])/2.)

			D1_vec = map(x->sum((Df_geom*f_geom).*(x.^(0:length(geom_koeff[:,1])-1)),dims = 1),x_trans)
			D2_vec =  map(x->sum(((Df_geom^2)*f_geom).*(x.^(0:length(geom_koeff[:,1])-1)),dims = 1),x_trans) 
			
			Pi[i] = -map(v-> atan(v[2],v[1]),D1_vec)
			Ki[i] = map((v1,v2)-> abs(det([v1;v2]))/norm(v1)^3,D1_vec,D2_vec)		
			
			Kb[i]= Ki[i][1:2]
			Pb[i] = Pi[i][1:2]

			Ki[i] = Ki[i][3:end]
			Pi[i] = Pi[i][3:end]
			Li[i] = (sum(norm.(D1_vec[3:end]).*wg[i])*(elem_dat.div1[i+1]-elem_dat.div1[i])/2.)

		end
	


		

		#koeficienti razoja geometrije
		# Interpolacijska baza za vsak končni element
		Ib = map(i -> (!elem_dat.Ci) ? re_gramschmid([collect(range(-1.,1.,length = elem_dat.div2[i]))]) : (1<i<n_ke ? re_gramschmid([collect(range(-1.,1.,length = elem_dat.div2[i])),[-1.,1.]]) : ( i==1 ? re_gramschmid([collect(range(-1.,1.,length = elem_dat.div2[i])),[-1.]]) : re_gramschmid([collect(range(-1.,1.,length = elem_dat.div2[i])),[1.]]))),1:n_ke)


		return BeamDataProcess(Li,Pi,Ki,Ib,Pb,Kb,xg,wg,indx), n_nodes	
	end


	# Funkcija za branje datoteka z podatki
	function readdata()::Tuple{String,String,Array{String}}
		print("Pot do datoteke z podatki: ")
		file = readline()
		
		if isempty(findall(".txt",file))
			file = file*".txt"
		end

		data = readlines(file)
		data = data[setdiff(1:length(data),findall(sizeof.(data).==0))]
		data = data[findall(isempty.(findall.("#",replace.(data,"\t"=>""))))]


		data1 = data[1:findfirst(isempty.(findall.("elementi::Array",data)).==0) - 1]
		data2 = data[(length(data1)+1):findfirst(isempty.(findall.("@assignto",data)).==0)-1]
		data1 = prod(data1)
		data2 = prod(data2)
		data3 = data[setdiff(1:length(data),findall(isempty.(findall.("@assignto",data))))]
		return data1,data2,data3
    	end #readdata



	# Funkcije za interpolacijsko bazo	
    	function DotP(a::Array{Float64},DataIn::Vector{Vector{Float64}})
        	m = length(DataIn)
       	 	n = sum(length.(DataIn))-1
		#Diferencialni operator
        	Dp = diagm(1=>1:n)
		
		#Seznam funkcijskih vrednosti polinoma glede na DataIn
        	b = vcat(map( (A,i) -> (((Dp^i)*a)'* (A'.^(0:n)) )' ,DataIn,0:(m-1))...)
        	return b
    	end
	function DotP(a::Array{Float64},b::Array{Float64},DataIn::Vector{Vector{Float64}})
		#Skalarni produkt vektorjev iz DotP(a,A)
		return DotP(a::Array{Float64},DataIn::Vector{Vector{Float64}})'*DotP(b::Array{Float64},DataIn::Vector{Vector{Float64}})
	end
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
	function re_gramschmid(DataIn::Vector{Vector{Float64}})
		bi = gramschmid(DataIn)[4]
		bi = gramschmid(DataIn;B0 = bi)[1]
		if all(bi[end,:].<10^(-13))
			error("Interpolacijska baza ne obstaja")
		end
		return bi
	end






	#Funkcije za izvrednotenje linearne kombinacije ali polinom sam
	function InterpolValue(x::Float64,Kb::Vector{Float64},Ib::Matrix{Float64};n::Int64=0)::Float64
		Df = diagm(1 => 1. :length(Kb)-1.)^n
		f = (Df*Ib*Kb)'*x.^(0:length(Kb)-1)
		return f
	end
	function PolyValue(x::Float64,Ki::Vector{Float64};n::Int64 = 0)
		m = length(Ki)
		a = n
		f = InterpolValue(x,Ki,Matrix{Float64}(I,(m,m)); n = a)
		return f
	end
	




	#Funkcija za račun količin v integracijskih točkah nosilca
	function VarsAtX(x::Float64,ux::Array{Float64},uz::Array{Float64},phi::Array{Float64},vx::Array{Float64},vz::Array{Float64},omg::Array{Float64},Ib::Matrix{Float64},p0::Float64,k0::Float64,C::Matrix{Float64})

		U = map(qi -> InterpolValue(x,qi,Ib),[ux,uz,phi])	
		V = map(qi -> InterpolValue(x,qi,Ib),[vx,vz,omg])
		
		dU = map(qi -> InterpolValue(x,qi,Ib;n=1),[ux,uz,phi])	
		dV = map(qi -> InterpolValue(x,qi,Ib;n=1),[vx,vz,omg])
		

		D = dU+[cos(p0); sin(p0); 0.0]
		E = R(U[3]+p0)*D+[-1.0;0.0;-k0]
		Re = R(U[3]+p0)'* C* E
		
		# (dlPhi , dlD)
		dlRe = ( R(U[3]+p0;n=1)'*C*E + R(U[3]+p0)'*C*R(U[3]+p0;n=1)*D, R(U[3]+p0)'*C*R(U[3]+p0) )
		return U,V,dU,dV,D,Re,dlRe
		
	end
	

	# Na KE
	function Tan_Res(xInt::Vector{Float64},wInt::Vector{Float64},ux::Matrix{Float64},uz::Matrix{Float64},phi::Matrix{Float64},vx::Matrix{Float64},vz::Matrix{Float64},omg::Matrix{Float64},Ib::Matrix{Float64},p0::Vector{Float64},k0::Vector{Float64},C::Matrix{Float64},M::Vector{Float64},Fpx::Vector{Float64},Fpz::Vector{Float64} ,Fmy::Vector{Float64},dt::Float64,pb::Vector{Float64},kb::Vector{Float64},L::Float64,g::Float64)

		ux1 = ux[:,1]; uz1 = uz[:,1]; phi1 = phi[:,1];
		vx1 = vx[:,1]; vz1 = vz[:,1]; omg1 = omg[:,1]; vx2 = vx[:,2]; vz2 = vz[:,2]; omg2 = omg[:,2]
		#ux,uz,... so vrednosti v interpolacijskih točkah
		
		F = fill(Vector{Float64}([0.0;0.0;0.0]),length(ux1))
		dlF = fill(zeros(Float64,(3,3)),(length(ux1),length(ux1)))

		indx2 = CartesianIndex.((1:length(Ib[:,1]))',1:length(Ib[:,1]))
		for i1 = eachindex(xInt)
			#Poračunaj količine v xg
				#za čas tn in tn+1
			
			V1 = map(vi->InterpolValue(xInt[i1],vi,Ib),[vx1,vz1,omg1])
			V2 = map(vi->InterpolValue(xInt[i1],vi,Ib),[vx2,vz2,omg2])
			dV1 = map(vi->InterpolValue(xInt[i1],vi,Ib;n=1),[vx1,vz1,omg1])
			dV2 = map(vi->InterpolValue(xInt[i1],vi,Ib;n=1),[vx2,vz2,omg2])

			V = (V1+V2)/2
			dV = (dV1+dV2)/2


			U1 = map(ui->InterpolValue(xInt[i1],ui,Ib),[ux1,uz1,phi1])
			dU1 = map(ui->InterpolValue(xInt[i1],ui,Ib;n=1), [ux1,uz1,phi1])

			U = U1+V*dt/2
			U2 = U1+V*dt
			
			D1 = dU1+[cos(p0[i1]);sin(p0[i1]);0.0]
			D = D1+dt/2*dV
			D2 = D1+dt*dV


			E = R(U[3]+p0[i1])*D + [-1.0; 0.0; -k0[i1]]

			#E1 = R(U1[3])*D1 + [-1.0;0.0;-k0[i1]]
			#E2 = R(U2[3])*D2 + [-1.0;0.0;-k0[i1]] 
			
			Re = R(U[3]+p0[i1])'*C*E
			#Re1 = R(U1[3])*C*E1
			#Re2 = R(U2[3])*C*E2
			#Re = (Re1+Re2)/2
		
			N = C*E
			
			dlE = (R(U[3]+p0[1i1];n=1)*D , R(U[3]+p0[i1]))
			dlRe =(R(U[3]+p0[i1];n=1)'*C*E + R(U[3]+p0[i1])'*C*R(U[3]+p0[i1];n=1)*D , R(U[3]+p0[i1])'*C*R(U[3]+p0[i1]))
			dlN = (C*R(U[3]+p0[i1];n=1)*E , C*R(U[3]+p0[i1]))
			

			#=
			U1,V1,dU1,dV1,D1,Re1,dlRe1 = VarsAtX(xInt[i1],ux1,uz1,phi1,vx1,vz1,omg1,Ib,p0[i1],k0[i1],C)
			   V2,    dV2		   = VarsAtX(xInt[i1],ux2,uz2,phi2,vx2,vz2,omg2,Ib,p0[i1],k0[i1],C)[[2,4]]
			

			#Količine v času tn+1/2
			V = (V1 + V2)/2.
			dV = (dV1+dV2)/2.
			
			U = U1 + V*dt/2.
			dU = dU1 + dV*dt/2.
			D = D1 + dV*dt/2.


			Re = (Re1 + Re2)/2.
			dlRe = (dlRe1 .+ dlRe2)./2.
			
			E = R(U[3]+p0[i1])*D+[-1.0;0.0;-k0[i1]]
			N = C*E
			dlE = ( R(U[3]+p0[i1];n=1)*D, R(U[3]+p0[i1]) )
			dlN = (C*dlE[1], C*dlE[2])
			=#
			
			# Obtežba v xInt za vmesni čas
			p = map(pj -> PolyValue(xInt[i1],[0.5 0.5;-0.5 0.5]*reshape(pj,(2))),[Fpx,Fpz,Fmy])
			# Rabim se vektor pospeskov v xInt
			tv =(V2-V1)/dt+ [0.0; g; 0.0]
			

			F .+= map(i2 -> -Re*PolyValue(xInt[i1],Ib[:,i2];n=1) + (p + [0.0;0.0;dot(N,[E[2];-(1.0 +E[1]);0.0])] - tv.*[M[1];M[1];M[2]] )*PolyValue(xInt[i1],Ib[:,i2]) ,1:length(Ib[:,1])) * wInt[i1]	
	
			# dlRe
			dlF .+= map( ij -> PolyValue(xInt[i1],Ib[:,ij[1]];n=1) * ([[0.0;0.0;0.0] [0.0;0.0;0.0] dlRe[1]]*PolyValue(xInt[i1],Ib[:,ij[2]]) + dlRe[2]*PolyValue(xInt[i1],Ib[:,ij[2]];n=1)) ,indx2) * wInt[i1]*dt/2.
			# dlN[1] dlE[1]
			dlF .+= map( ij ->  PolyValue(xInt[i1],Ib[:,ij[1]])*([[0.0 0.0 0.0; 0.0 0.0 0.0];[0.0 0.0 (E[2]*dlN[1][1] - (1+E[1])*dlN[1][2] + N[1]*dlE[1][2] - N[2]*dlE[2][1])]])*PolyValue(xInt[i1],Ib[:,ij[2]]), indx2 )*wInt[i1]*dt/2.
			# dlN[2] dlE[2]
			dlF .+= map( ij ->  PolyValue(xInt[i1],Ib[:,ij[1]])*PolyValue(xInt[i1],Ib[:,ij[2]];n=1)*[[0.0 0.0 0.0; 0.0 0.0 0.0]; (E[2]*dlN[2][1,:] - (1+E[1])*dlN[2][2,:] - N[2]*dlE[2][1,:] + N[1]*dlE[2][2,:])'] ,indx2)*wInt[i1]*dt/2.
		end

		
		xb = [-1.,1.]
		for i in 1:2 
			#Poračunaj količine v x
				#za čas tn in tn+1
			#Treba je zračunat kot in ukrivljenost v robnih točkah		
	
			
			V1 = map(vi->InterpolValue(xb[i],vi,Ib),[vx1,vz1,omg1])
			V2 = map(vi->InterpolValue(xb[i],vi,Ib),[vx2,vz2,omg2])
			dV1 = map(vi->InterpolValue(xb[i],vi,Ib;n=1),[vx1,vz1,omg1])
			dV2 = map(vi->InterpolValue(xb[i],vi,Ib;n=1),[vx2,vz2,omg2])

			V = (V1+V2)/2
			dV = (dV1+dV2)/2


			U1 = map(ui->InterpolValue(xb[i],ui,Ib),[ux1,uz1,phi1])
			dU1 = map(ui->InterpolValue(xb[i],ui,Ib;n=1), [ux1,uz1,phi1])

			U = U1+V*dt/2
			U2 = U1+V*dt
			
			D1 = dU1+[cos(pb[i]);sin(pb[i]);0.0]
			D = D1+dt/2*dV
			D2 = D1+dt*dV


			E = R(U[3]+pb[i])*D + [-1.0; 0.0; -kb[i]]

			#E1 = R(U1[3])*D1 + [-1.0;0.0;-k0[i1]]
			#E2 = R(U2[3])*D2 + [-1.0;0.0;-k0[i1]] 
			
			Re = R(U[3]+pb[i])'*C*E		
			#=
			U1,V1,dU1,dV1,D1,Re1,dlRe1 = VarsAtX(xb[i],ux1,uz1,phi1,vx1,vz1,omg1,Ib,pb[i],kb[i],C)
			U2,V2,dU2,dV2,D2,Re2,dlRe2 = VarsAtX(xb[i],ux2,uz2,phi2,vx2,vz2,omg2,Ib,pb[i],kb[i],C)
	
			#Količine v času tn+1/2
			Re = (Re1 + Re2)/2.
			=#
			F .+= map(i2 -> Re*PolyValue(xb[i],Ib[:,i2]) ,1:length(Ib[:,1])) 
		end

		F.*=L/2.
		dlF.*= L/2.

		#display(dlF[1])
		#display(F[1])
		return dlF,F

	end


end # module
