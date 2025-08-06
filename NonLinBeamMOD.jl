module NonLinBeam
    
	export Beam, BeamDataIn, BeamDataProcess, Node, NodeDataIn,Motion, BeamMotion,
		datainit, readdata, dataprocess,
		R2, R, GaussInt, InterpolKoeff, Interpolated, JacKoeff, BalanceEq, evaluateKoefficient, re_gramschmid


    using LinearAlgebra


    # Struktura podatkov za nosilec
    abstract type Beam end

    @kwdef mutable struct BeamDataIn <:Beam
        v::Vector{Int64} = [1;2] # Vozlišča - krajna
        
	C::Matrix{Float64} = [1 0 0;0 1 0;0 0 1] #Materialna matrika
        M::Vector{Float64} = [1; 1] #Vektor [ρA; ρI]

	Ib_geom::Matrix{Float64} = [0.5 0.5; -0.5 0.5] #re_gramschmid(DataIn::Vector{Vector{Float64}})
	Kb::Matrix{Float64} = Array{Float64,2}(undef,(0,2)) 
	
        px::Function = t->[0. 0.] 
	pz::Function = t->[0. 0.]
        my::Function = t->[0. 0.]

        div1::Array{Float64} = [-1.,1.]
	#Naj bo vedno med -1 in 1

        div2::Array{Int64} = [4]
	#Število vozlišč v sekundarni delitvi. Dve sta robni. Int za vsak element.

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
    end
    


    abstract type Node 
    end

    @kwdef mutable struct NodeDataIn <:Node
        # te poračuna algoritem
	x::Float64 = 0.
        y::Float64 = 0.
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




    
    #Funkcija za rotacijo 2-terice za kot
    function R2(a)::Matrix{Float64} 
	    return [cos(a) -sin(a); sin(a) cos(a)]
    end
    #Funkcija za rotacijo 3-terice za kot a okrov e3
    function R(a;n = 0)::Matrix{Float64}
	    # n je stopnja odvoda
	    return [0. -1. 0.;1. 0. 0.; 0. 0. 0.]^n * [cos(a) -sin(a) 0. ; sin(a) cos(a) 0. ; 0. 0. 1.]
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
	function datainit(elem::Matrix{Int64},voz::Matrix{Float64})
		n_elem::Int64 = length(elem[:,1])
		n_voz::Int64 = length(voz[:,1])


		element_data::Vector{BeamDataIn} = fill(BeamDataIn(),n_elem)
		voz_data::Vector{NodeDataIn} = fill(NodeDataIn(),n_voz)


		for i =1:n_elem
		    element_data[i] = BeamDataIn(v = elem[i,[1,2]]) 
		end

		for i =1:n_voz
		    voz_data[i] = NodeDataIn(x=voz[i,1],y=voz[i,2],i=i)
		end

		return n_elem,n_voz,element_data,voz_data
	end #datainit

	function dataprocess(elem_dat::BeamDataIn,node_dat::Array{NodeDataIn})::BeamDataProcess
		node1 = [node_dat[1].x node_dat[1].y]
		node2 = [node_dat[2].x node_dat[2].y]
		n_ke = length(elem_dat.div1)-1
	

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

			x_trans = xg[i]/2. *(elem_dat.div1[i+1]-elem_dat.div1[i]).+(elem_dat.div1[i+1]+elem_dat.div1[i])/2.

			D1_vec = map(x->sum((Df_geom*f_geom).*(x.^(0:length(geom_koeff[:,1])-1)),dims = 1),x_trans)
			D2_vec =  map(x->sum(((Df_geom^2)*f_geom).*(x.^(0:length(geom_koeff[:,1])-1)),dims = 1),x_trans) 
			
			Pi[i] = map(v-> atan(v[2],v[1]),D1_vec)
			Ki[i] = map((v1,v2)-> abs(det([v1;v2]))/norm(v1)^3,D1_vec,D2_vec)		

			Li[i] = sum(norm.(D1_vec).*wg[i])*(elem_dat.div1[i+1]-elem_dat.div1[i])/2.

		end
	


		

		#koeficienti razoja geometrije
		# Interpolacijska baza za vsak končni element
		Ib = map(i -> (!elem_dat.Ci) ? re_gramschmid([collect(range(-1.,1.,length = elem_dat.div2[i]))]) : (1<i<n_ke ? re_gramschmid([collect(range(-1.,1.,length = elem_dat.div2[i])),[-1.,1.]]) : ( i==1 ? re_gramschmid([collect(range(-1.,1.,length = elem_dat.div2[i])),[-1.]]) : re_gramschmid([collect(range(-1.,1.,length = elem_dat.div2[i])),[1.]]))),1:n_ke)


		return BeamDataProcess(Li,Pi,Ki,Ib,xg,wg)	
	end



	function readdata()::Tuple{String,String,Array{String}}
		print("Pot do datoteke z podatki:")
		file = readline()
		
		if isempty(findall(".txt",file))
			file = file*".txt"
		end

		data = readlines(file)
		data = data[setdiff(1:length(data),findall(sizeof.(data).==0))]
		data = data[findall(isempty.(findall.("#",replace.(data,"\t"=>""))))]
		print(data)

		data1 = data[1:findfirst(isempty.(findall.("elementi::Array",data)).==0) - 1]
		data2 = data[(length(data1)+1):findfirst(isempty.(findall.("@assignto",data)).==0)-1]
		data1 = prod(data1)
		data2 = prod(data2)
		data3 = data[setdiff(1:length(data),findall(isempty.(findall.("@assignto",data))))]

		return data1,data2,data3
    	end #readdata



    #Funkcija za določitev koeficientov standardne baze za interpolacijo v danih točkah
    #
    #Tole zamenjaj za gramshchmita
    #Funkcija za račun funkcijske vrednosti v x za interpolacijsko bazo Ib in koeficiente razvoja Kn

    #tole malo polepšaj
    

    Interpolated(x::Float64,Kn::Array{Float64},Ib::Vector{Array{Float64}}) = ((Kn' *Ib)'*(x.^(0:(length(Ib)-1))) )
    Interpolated(x::Float64,Kn::Float64,Ib::Array{Float64}) = (Kn*Ib)' * (x.^(0:(length(Ib)-1))) 
    Interpolated(x::Array{Float64},Kn::Array{Float64},Ib::Vector{Array{Float64}}) = map(x->Interpolated(x::Float64,Kn::Array{Float64},Ib::Vector{Array{Float64}}),x)
    Interpolated(x::Array{Float64},Kn::Float64,Ib::Array{Float64}) = map(x->Interpolated(x::Float64,Kn::Float64,Ib::Array{Float64}),x)





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
		return bi
	end


	function InterpolValue(x::Float64,Kb::Vector{Float64},Ib::Matrix{Float64};n::Int64=0)::Float64
		Df = diagm(1=>1. :length(Kb)-1.)^n
		f = (Df*Ib*Kb)'*x.^(0:length(Kb)-1)
		return f
	end

	# TOLE JE NASTY |
	# 		ˇ


	function VarsAtX(x::Float64,ux::Array{Float64},uz::Array{Float64},phi::Array{Float64},vx::Array{Float64},vz::Array{Float64},omg::Array{Float64},Ib::Matrix{Float64},p0::Float64,k0::Float64,C::Matrix{Float64})

		U = map(qi -> InterpolaValue(x,qi,Ib),[ux,uz,phi])	
		V = map(qi -> InterpolaValue(x,qi,Ib),[vx,vz,omg])
		
		dU = map(qi -> InterpolaValue(x,qi,Ib;n=1),[ux,uz,phi])	
		dV = map(qi -> InterpolaValue(x,qi,Ib;n=1),[vx,vz,omg])
		
		D = U+[cos(p0); sin(p0); 0.0]
		
		dlD = Matrix(I,(3,3))

		E = R(U[3])*D + [-1.0 ; 0.0; k0]
		R = R(U[3])*C*E
		
		dlE = ( R(U[3];n=1)*D , R(U[3])*dlD )
		dlR = ( R(U[3];n=1)*C*E + R(U[3])*C*R(U[3];n=1)*D, R(U[3])*dlD)

		
		
		
	end


    function evaluateKoefficient(ux,uz,phi,vx,vz,Omg,phi0,L,C,px,pz,my,g,xInt,wInt,dt,P,dP,A1,A2,A3,A4,A5,A6)
        #=  Popravi
            ‖_ Hitrosti so povprečene pred klicom funkcije
            ‖_ Preveri kaj lahko narediš z matričnim računom in rotacijo R(ϕ) 
            ‖_ 
        =#


        x_Int = vcat(-1.,xInt,1.)
        
        #pomiki v zečetnem času
        uxₙ₁ = Interpolated(x_Int,ux,P)
        uzₙ₁ = Interpolated(x_Int,uz,P)
        ϕₙ₁ = Interpolated(x_Int,phi,P)
        duxₙ₁ = Interpolated(x_Int,ux,dP)
        duzₙ₁ = Interpolated(x_Int,uz,dP)
        dϕₙ₁ = Interpolated(x_Int,phi,dP)


        #povprečne hitrosti
        v̄x = Interpolated(x_Int,vx,P)
        v̄z = Interpolated(x_Int,vz,P)
        Ω = Interpolated(x_Int,Omg,P)
        dv̄x = Interpolated(x_Int,vx,dP)
        dv̄z = Interpolated(x_Int,vz,dP)
        dΩ = Interpolated(x_Int,Omg,dP)
        

        #defomacije v fiksni bazi 
        Dxₙ₁ = cos(phi0) .+ duxₙ₁
        Dzₙ₁ = sin(phi0) .+ duzₙ₁
        Dxₙ₂ = Dxₙ₁ + dt/2*dv̄x
        Dzₙ₂ = Dzₙ₁ + dt/2*dv̄z
        Dxₙ₃ = Dxₙ₁ + dt*dv̄x
        Dzₙ₃ = Dzₙ₁ + dt*dv̄z


        
        Einit = [-1.;0.;0.]
        # Začetni časp 
        Rot0 = R.(ϕₙ₁)
        dRot0 = dR.(ϕₙ₁)
        D0 = map( k -> [Dxₙ₁[k]; Dzₙ₁[k];dϕₙ₁[k]], eachindex(x_Int))
        
        # Vmesni čas
        RotM = R.(ϕₙ₁ + dt/2*Ω)
        dRotM = dR.(ϕₙ₁ + dt/2*Ω)
        DM = map( k -> [Dxₙ₂[k]; Dzₙ₂[k];dϕₙ₁[k]+dt/2*dΩ[k]], eachindex(x_Int))
        
        # Končni čas
        Rot1 = R.(ϕₙ₁ + dt*Ω)
        dRot1 = dR.(ϕₙ₁ + dt*Ω)
        D1 = map( k -> [Dxₙ₃[k]; Dzₙ₃[k];dϕₙ₁[k]+dt*dΩ[k]], eachindex(x_Int))

        #Vstaviko količine iz vmesenga časa
        #[ δΩ ; δvx' ; δvz' ; δΩ']

        EM = RotM.*DM .+ [Einit]
        RiM = RotM .* map(k->C*k,EM)
        NiM = map(k-> C*k ,EM)
        δEM = map(k-> k*[dt/2. 0. 0. 0.] ,dRotM .*DM ) + map(k->k*[[0. dt/2 0. 0.];[0. 0. dt/2 0.];[0. 0. 0. dt/2]],RotM)
        δRiM =  map(j -> j*[dt/2. 0. 0. 0.] , dRotM.*map(i-> C*i,RotM.*DM ) + RotM.*map(i-> C*i,dRotM.*DM ) + map(i->i*C*Einit,dRotM))  + map(k->k*[[0. dt/2 0. 0.];[0. 0. dt/2 0.];[0. 0. 0. dt/2]],RotM.*map(i-> C*i,RotM ))
            
        #= Print
            println("δRiM[1] = ")  
            display(δRiM[1])
            println("RiM[1] = ")  
            display(RiM[1])
            println("NiM[1] = ")  
            display(NiM[1])
            println("δEM[1] = ")  
            display(δEM[1])
            println("EM[1] = ")  
            display(EM[1])
        =#

        RxM = vcat(map(i->-RiM[i][1],eachindex(x_Int))...)
        RzM = vcat(map(i->-RiM[i][2],eachindex(x_Int))...)
        MyM = vcat(map(i->-RiM[i][3],eachindex(x_Int))...)


        # Ni še inercijskih sil
        fxi = map(j->RxM[2:end-1]'*A6[j] ,eachindex(A6)) + map(j-> ((px[1]+px[2] .+xInt*(px[2]-px[1]))/2)'*A5[j] ,eachindex(A5)) + map(i-> RxM[[1,end]]' *Interpolated([-1.,1.],1.,P[i])  , eachindex(P))
        fzi = map(j->RzM[2:end-1]'*A6[j] ,eachindex(A6)) + map(j-> ((pz[1]+pz[2] .+xInt*(pz[2]-pz[1]))/2)'*A5[j] ,eachindex(A5)) + map(i-> RzM[[1,end]]' *Interpolated([-1.,1.],1.,P[i])  , eachindex(P))
        fmi = map(j->MyM[2:end-1]'*A6[j] ,eachindex(A6)) + map(j-> (RxM[2:end-1].*(duzₙ₁+dt/2*̄dv̄z) - (1. .+duxₙ₁+dt/2*dv̄x).*RzM[2:end-1] + ((my[1]+my[2] .+xInt*(my[2]-my[1]))/2))'*A6[j] ,eachindex(A6)) + map(i-> MyM[[1,end]]' *Interpolated([-1.,1.],1.,P[i])  , eachindex(P))

        
        δfxi = map(j->RxM[2:end-1]'*A6[j] ,eachindex(A6)) 
        δfzi = map(j->RzM[2:end-1]'*A6[j] ,eachindex(A6)) 
        δfmi = map(j->MyM[2:end-1]'*A6[j] ,eachindex(A6)) 


        #=
            #deformacije ob vmesnem času
            ϵₙ₂ =  Dxₙ₂.*cos.(ϕₙ₁ + dt/2*Ω) + Dzₙ₂.*sin.(ϕₙ₁ + dt/2*Ω) .-1
            γₙ₂ = -Dxₙ₂.*sin.(ϕₙ₁ + dt/2*Ω) + Dzₙ₂.*cos.(ϕₙ₁ + dt/2*Ω)
            κₙ₂ =  dϕₙ₁ + dt/2*dΩ


            #rezultante napetosti
            Rxₙ₂ = cos.(ϕₙ₁ + dt/2*Ω).*(C[1,1]*ϵₙ₂ + C[1,2]*γₙ₂ + C[1,3]*κₙ₂) + sin.(ϕₙ₁ + dt/2*Ω).*(C[2,1]*ϵₙ₂ + C[2,2]*γₙ₂ + C[2,3]*κₙ₂)   
            Rzₙ₂ = -sin.(ϕₙ₁ + dt/2*Ω).*(C[1,1]*ϵₙ₂ + C[1,2]*γₙ₂ + C[1,3]*κₙ₂) + cos.(ϕₙ₁ + dt/2*Ω).*(C[2,1]*ϵₙ₂ + C[2,2]*γₙ₂ + C[2,3]*κₙ₂)   
            Mcₙ₂ = C[3,1]*ϵₙ₂ + C[3,2]*γₙ₂ + C[3,3]*κₙ₂


            

            δRxₙ₂ = [[dt/2*cos.(ϕₙ₁ + dt/2*Ω).*(C[1,1]*cos.(ϕₙ₁ + dt/2*Ω)-C[1,2]*sin.(ϕₙ₁ + dt/2*Ω)) + sin.(ϕₙ₁ + dt/2*Ω).*(C[2,1]*cos.(ϕₙ₁ + dt/2*Ω)-C[2,2]*sin.(ϕₙ₁ + dt/2*Ω))];
                    [dt/2*cos.(ϕₙ₁ + dt/2*Ω).*(C[1,1]*sin.(ϕₙ₁ + dt/2*Ω)+C[1,2]*cos.(ϕₙ₁ + dt/2*Ω)) + sin.(ϕₙ₁ + dt/2*Ω).*(C[2,1]*sin.(ϕₙ₁ + dt/2*Ω)+C[2,2]*cos.(ϕₙ₁ + dt/2*Ω))];
                    [dt/2*((C[1,1]*(cos.(ϕₙ₁ + dt/2*Ω).*(sin(phi0) .+ duzₙ₁ + dt/2*dv̄z) - sin.(ϕₙ₁ + dt/2*Ω).*(cos(phi0).+duxₙ₁+dt/2*dv̄x)) + C[1,2]*(-sin.(ϕₙ₁ + dt/2*Ω).*(sin(phi0) .+ duzₙ₁ + dt/2*dv̄z) + cos.(ϕₙ₁ + dt/2*Ω).*(cos(phi0).+duxₙ₁+dt/2*dv̄x) )).*cos.(ϕₙ₁ + dt/2*Ω) + (C[2,1]*(cos.(ϕₙ₁ + dt/2*Ω).*(sin(phi0) .+ duzₙ₁ + dt/2*dv̄z) - sin.(ϕₙ₁ + dt/2*Ω).*(cos(phi0).+duxₙ₁+dt/2*dv̄x)) + C[2,2]*(-sin.(ϕₙ₁ + dt/2*Ω).*(sin(phi0) .+ duzₙ₁ + dt/2*dv̄z) + cos.(ϕₙ₁ + dt/2*Ω).*(cos(phi0).+duxₙ₁+dt/2*dv̄x) )).*sin.(ϕₙ₁ + dt/2*Ω))];
                    [dt/2*(C[1,3]*cos.(ϕₙ₁ + dt/2*Ω)+C[2,3]*sin.(ϕₙ₁ + dt/2*Ω))]
            ]
            
            δRzₙ₂ = [[dt/2*(-sin.(ϕₙ₁ + dt/2*Ω).*(C[1,1]*cos.(ϕₙ₁ + dt/2*Ω)-C[1,2]*sin.(ϕₙ₁ + dt/2*Ω)) + cos.(ϕₙ₁ + dt/2*Ω).*(C[2,1]*cos.(ϕₙ₁ + dt/2*Ω)-C[2,2]*sin.(ϕₙ₁ + dt/2*Ω)))];
                    [dt/2*(-sin.(ϕₙ₁ + dt/2*Ω).*(C[1,1]*sin.(ϕₙ₁ + dt/2*Ω)+C[1,2]*cos.(ϕₙ₁ + dt/2*Ω)) + cos.(ϕₙ₁ + dt/2*Ω).*(C[2,1]*sin.(ϕₙ₁ + dt/2*Ω)+C[2,2]*cos.(ϕₙ₁ + dt/2*Ω)))];
                    [-dt/2*((C[1,1]*(cos.(ϕₙ₁ + dt/2*Ω).*(sin(phi0) .+ duzₙ₁ + dt/2*dv̄z) - sin.(ϕₙ₁ + dt/2*Ω).*(cos(phi0).+duxₙ₁+dt/2*dv̄x)) + C[1,2]*(-sin.(ϕₙ₁ + dt/2*Ω).*(sin(phi0) .+ duzₙ₁ + dt/2*dv̄z) + cos.(ϕₙ₁ + dt/2*Ω).*(cos(phi0).+duxₙ₁+dt/2*dv̄x) )).*sin.(ϕₙ₁ + dt/2*Ω) + (C[2,1]*(cos.(ϕₙ₁ + dt/2*Ω).*(sin(phi0) .+ duzₙ₁ + dt/2*dv̄z) - sin.(ϕₙ₁ + dt/2*Ω).*(cos(phi0).+duxₙ₁+dt/2*dv̄x)) + C[2,2]*(-sin.(ϕₙ₁ + dt/2*Ω).*(sin(phi0) .+ duzₙ₁ + dt/2*dv̄z) + cos.(ϕₙ₁ + dt/2*Ω).*(cos(phi0).+duxₙ₁+dt/2*dv̄x) )).*cos.(ϕₙ₁ + dt/2*Ω))];
                    [dt/2*(-C[1,3]*sin.(ϕₙ₁ + dt/2*Ω)+C[2,3]*cos.(ϕₙ₁ + dt/2*Ω))]
            ]

            δMcₙ₂ = [[dt/2*(C[3,1]*cos.(ϕₙ₁ + dt/2*Ω) - C[3,2]*sin.(ϕₙ₁ + dt/2*Ω))] ;
                    [dt/2*(C[3,1]*sin.(ϕₙ₁ + dt/2*Ω) + C[3,2]*cos.(ϕₙ₁ + dt/2*Ω))];
                    [dt/2*((C[3,1]*(cos.(ϕₙ₁ + dt/2*Ω).*(sin(phi0) .+ duzₙ₁ + dt/2*dv̄z) - sin.(ϕₙ₁ + dt/2*Ω).*(cos(phi0).+duxₙ₁+dt/2*dv̄x)) + C[3,2]*(-sin.(ϕₙ₁ + dt/2*Ω).*(sin(phi0) .+ duzₙ₁ + dt/2*dv̄z) + cos.(ϕₙ₁ + dt/2*Ω).*(cos(phi0).+duxₙ₁+dt/2*dv̄x) )))];
                    dt/2*C[3,3]
            ]


            
            #variacije rezultant
            
            #println(A3)

            Jxx = map(a -> -δRxₙ₂[1][2:end-1]'* A3[a],   CartesianIndex.(1:length(P),(1:length(P))'))
            Jxz = map(a -> -δRxₙ₂[2][2:end-1]'* A3[a],   CartesianIndex.(1:length(P),(1:length(P))'))
            JxO = map(a -> -δRxₙ₂[3][2:end-1]'* A1[a] - δRxₙ₂[4][2:end-1]'* A3[a],   CartesianIndex.(1:length(P),(1:length(P))'))
            
            Jzx = map(a -> -δRzₙ₂[1][2:end-1]'* A3[a],   CartesianIndex.(1:length(P),(1:length(P))'))
            Jzz = map(a -> -δRzₙ₂[2][2:end-1]'* A3[a],   CartesianIndex.(1:length(P),(1:length(P))'))
            JzO = map(a -> -δRzₙ₂[3][2:end-1]'* A1[a] - δRxₙ₂[4][2:end-1]'* A3[a],   CartesianIndex.(1:length(P),(1:length(P))'))

            #println(-δMcₙ₂)
            JOx = map(a -> -δMcₙ₂[1][2:end-1]'* A3[a] .+ (δRxₙ₂[1].*(duzₙ₁ + dt/2*dv̄z)-δRzₙ₂[1].*(1. .+duxₙ₁+dt/2*dv̄x)-Rzₙ₂)[2:end-1]'* A2[a] ,   CartesianIndex.(1:length(P),(1:length(P))'))
            #println(JOx)
            JOz = map(a -> -δMcₙ₂[2][2:end-1]'* A3[a] +(δRxₙ₂[2].*(duzₙ₁ + dt/2*dv̄z)-δRzₙ₂[2].*(1. .+duxₙ₁+dt/2*dv̄x)+Rxₙ₂)[2:end-1]'* A2[a] ,   CartesianIndex.(1:length(P),(1:length(P))'))
            JOO = map(a -> -δMcₙ₂[3][2:end-1]'* A1[a] - sum(δMcₙ₂[4]* A3[a]) + (δRxₙ₂[3].*(duzₙ₁+dt/2*dv̄z)-δRzₙ₂[3].*(1. .+duxₙ₁+dt/2*dv̄x))[2:end-1]'* A2[a] + (δRxₙ₂[4].*(duzₙ₁+dt/2*dv̄z)-δRzₙ₂[4].*(1. .+duxₙ₁+dt/2*dv̄x))[2:end-1]'* A4[a],   CartesianIndex.(1:length(P),(1:length(P))'))
            
            J= [Jxx Jxz JxO; Jzx Jzz JzO; JOx JOz JOO]
            

            #fxi = map( i -> wInt'*( -Rxₙ₂.*Interpolated(xInt,1., dP[i]) +  ((px[1]+px[2] .+xInt*(px[2]-px[1]))/2 - (vxInt[2:end-1] -Interpolated(xInt,vx0,P))/dt*M[1])  .* Interpolated(xInt,1.,P[i]) )  + Rx0' *Interpolated([-1.,1.],1.,P[i])  , eachindex(P))
            #fzi = map( i -> wInt'*( -Rzₙ₂.*Interpolated(xInt,1., dP[i]) +  ((pz[1]+pz[2] .+xInt*(pz[2]-pz[1]))/2 - ((vzInt[2:end-1]-Interpolated(xInt,vz0,P))/dt.-g)*M[1])  .*Interpolated(xInt,1.,P[i]) )  + Rz0' *Interpolated([-1.,1.],1.,P[i]) , eachindex(P))
            #fpi = map( i -> wInt'*( -Mcₙ₂.*Interpolated(xInt,1., dP[i]) +  (Rx.*(duzInt[2:end-1]+dt/2*dvzInt[2:end-1]) -Rz.*(1 .+duxInt[2:end-1]+dt/2*dvzInt[2:end-1])+  (my[1]+my[2] .+xInt*(my[2]-my[1]))/2 -  M[2]*(OmgInt[2:end-1]-Interpolated(xInt,Omg0,P))/dt  ).*Interpolated(xInt,1.,P[i]) )  +Mc0' *Interpolated([-1.,1.],1.,P[i]) , eachindex(P))

            #Mankajo še čelni obtežbe in pospeškov
            Fx = map(a -> -Rxₙ₂[2:end-1]' * A6[a] + ((px[1]+px[2] .+xInt*(px[2]-px[1]))/2)[2:end-1]'*A5[a] + Rxₙ₂[[1,end]]'*Interpolated([-1.,1.],1.,P[a]),   1:length(P))*L/2
            Fz = map(a -> -Rzₙ₂[2:end-1]' * A6[a] + ((pz[1]+pz[2] .+xInt*(pz[2]-pz[1]))/2 .- g)[2:end-1]'*A5[a] + Rzₙ₂[[1,end]]'*Interpolated([-1.,1.],1.,P[a]),   1:length(P))*L/2
            FO = map(a -> -Mcₙ₂[2:end-1]' * A6[a] + (Rxₙ₂.*(duzₙ₁+dt/2*dv̄z) - Rzₙ₂[1].*(1. .+duxₙ₁+dt/2*dv̄x) + (my[1]+my[2] .+xInt*(my[2]-my[1]))/2)[2:end-1]'*A5[a] + Mcₙ₂[[1,end]]'*Interpolated([-1.,1.],1.,P[a]),  1:length(P))*L/2

            F=[Fx;Fz;FO]
            
        
            return hcat(J,F)
        =#
    end

end # module
