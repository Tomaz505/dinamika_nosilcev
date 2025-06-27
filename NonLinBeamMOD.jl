module NonLinBeam
    
    export Beam,BeamDataIn,BeamDataProcess,
    Motion, BeamMotion,
    CMat,DiaMat,SymMat,AnyMat,CMat_class,
    Node, NodeDataIn,
    R2, GaussInt,InterpolKoeff,Interpolated,JacKoeff,BalanceEq,evaluateKoefficient


    using LinearAlgebra


    # Materialna matrika
    abstract type CMat 
    end
    struct DiaMat<:CMat
        Cd::Array{Float64}
    end
    struct SymMat<:CMat
        Cd::DiaMat
        Cs::Array{Float64}
    end
    struct AnyMat<:CMat
        #Cd::DiaMat
        Cs::SymMat
        Ca::Array{Float64}
    end
    function CMat_class(C::Array)
        C_diag =  C.*Matrix(I,3,3)
        C_sym = (C + C')/2 - C_diag
        C_ant_sym = (C - C')/2

        if all(C_sym.==0)
            CM = DiaMat(C_diag)
        elseif all(C_ant_sym.==0)
            CM = SymMat(DiaMat(C_diag),C_sym)
        else
            CM = AnyMat( SymMat( DiaMat(C_diag), C_sym),C_ant_sym)
        end
        return CM
    end



    # Struktura podatkov za nosilec
    abstract type Beam 
    end
    @kwdef mutable struct BeamDataIn <:Beam
        v::Array{Int64} = [1,1] # Vozlišča - krajna
        C::Array{Float64} = [1 0 0;0 1 0;0 0 1] #Materialna matrika
        M::Array{Float64} = [1; 1] #Vektor [ρA; ρI]
        f::Function = x->[0,x] #Funkcija krivulje
        df::Function = x->[0,1] #Odvod funkcije krivulje
        px::Function = t->[0. 0.] #Linearna distribucija obtežbe med robnima vrednostima [p1 p2] na poddelitvah
        pz::Function = t->[0. 0.]
        my::Function = t->[0. 0.]
        div1::Array{Float64} = [-1,1] #Relativne koordinate vozlišč primarne delitve.
        div2::Array{Int64} = [4] #Število vozlišč v sekundarni delitvi. Dve sta robni
        nInt::Array{Int64} = [100] #Število integracijskih točk v posameznem div1
    end 
    struct BeamDataProcess <:Beam
        C::Array{Float64}
        M::Array{Float64} 
        L::Array{Float64}
        ϕ0::Array{Float64}
        P::Array{Array{Array{Float64}}}
        dP::Array{Array{Array{Float64}}}
        A1::Array{Array{Array{Float64}}}
        A2::Array{Array{Array{Float64}}}
        A3::Array{Array{Array{Float64}}}
        A4::Array{Array{Array{Float64}}}
        A5::Array{Array{Array{Float64}}}
        A6::Array{Array{Array{Float64}}}
        xInt::Array{Array{Float64}}
        wInt::Array{Array{Float64}} # Jih dejansko rabim s sabo če so v A1,...? Lahko jih poračunam samo lokalno.
        v::Array{Array{Int64}}
    end
    


    abstract type Node 
    end
    @kwdef mutable struct NodeDataIn <:Node
        x::Float64 = 0.
        y::Float64 = 0.
        i::Int64 = 1
        Supp::Array{Bool} = [1 1 1]
    end



    abstract type Motion
    end
    mutable struct BeamMotion <:Motion
        ux::Array{Float64}
        uz::Array{Float64}
        phi::Array{Float64}
        vx::Array{Float64}
        vz::Array{Float64}
        Omg::Array{Float64}
    end




    # Funkcije
    #Funkcija za rotacijo 2-terice za kot a
    R2(a) = [cos(a) -sin(a); sin(a) cos(a)]
    R(a) = [cos(a) -sin(a) 0 ; sin(a) cos(a) 0 ; 0 0 1]
    dR(a) = [-sin(a) -cos(a) 0 ; cos(a) -sin(a) 0 ; 0 0 0]
    ddR(a) = [-cos(a) sin(a) 0 ; -sin(a) -cos(a) 0 ; 0 0 0]


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

    #Funkcija za določitev koeficientov standardne baze za interpolacijo v danih točkah
    function InterpolKoeff(IterpPoint::Union{Array,StepRangeLen,StepRange})
        n = length(IterpPoint)-1
        IterpPoint=reshape(collect(IterpPoint),(n+1,1))
        Dp = diagm(1=>1:n)
        DotP(a::Array)= (a' * (IterpPoint').^(0:n))'
        DotP(a::Array,b::Array) = (DotP(a))'*DotP(b)
    
    
        ei = map(k->Matrix(I,n+1,n+1)[k,(1:n+1)],1:n+1)
        Edots = hcat(map.(k->DotP(ei[k],ei),1:n+1)...)
        ei = map(k-> ei[k]./sqrt(Edots[k,k]),1:n+1)
    
        map!(k  ->  (ei[k] - (k==1 ? 0 .*ei[k] : sum( map(t->DotP(ei[k],ei[t]),1:k-1) .*ei[1:k-1]) )) / sqrt(DotP((ei[k] - (k==1 ? 0 .*ei[k] : sum( map(t->DotP(ei[k],ei[t]),1:k-1) .*ei[1:k-1]))) ,(ei[k] - (k==1 ? 0 .*ei[k] : sum( map(t->DotP(ei[k],ei[t]),1:k-1) .*ei[1:k-1]))) )) ,ei , 1:n+1)
    
        sp = length(DotP(ei[1]))
        Iunit = Matrix(I,sp,sp)
        UnitV = map(k->Iunit[k,(1:sp)],1:sp)
    
        Ib = map( t-> round.(sum( map(k->  ((DotP(ei[k]))'* UnitV[t]).*ei[k] ,  1:n+1)),digits=12) ,1:sp )
        dIb = map(k->Dp*k,Ib)
        return Ib,dIb
    end

    #Funkcija za račun funkcijske vrednosti v x za interpolacijsko bazo Ib in koeficiente razvoja Kn
    Interpolated(x::Float64,Kn::Array{Float64},Ib::Vector{Array{Float64}}) = ((Kn' *Ib)'*(x.^(0:(length(Ib)-1))) )
    Interpolated(x::Float64,Kn::Float64,Ib::Array{Float64}) = (Kn*Ib)' * (x.^(0:(length(Ib)-1))) 
    Interpolated(x::Array{Float64},Kn::Array{Float64},Ib::Vector{Array{Float64}}) = map(x->Interpolated(x::Float64,Kn::Array{Float64},Ib::Vector{Array{Float64}}),x)
    Interpolated(x::Array{Float64},Kn::Float64,Ib::Array{Float64}) = map(x->Interpolated(x::Float64,Kn::Float64,Ib::Array{Float64}),x)

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
        # Začetni čas
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

    #Metode funkcij za račun koeficientov pred integracijo koefcientov v Lineariziranih ravnotežnih enačbah
    function JacKoeff(C::DiaMat,ux::Array{Float64},uz::Array{Float64},phi::Array{Float64},vx::Array{Float64},vz::Array{Float64},Omg::Array{Float64},phi0::Float64,xInt::Array{Float64},dt::Float64,P::Vector{Array{Float64}},dP::Vector{Array{Float64}})
        C=C.Cd[CartesianIndex.(1:3,1:3)]
        
        

        #Mogoče bi lahko prej izračunal sin in cos...

        #uxInt = Interpolated(xInt,ux,P)
        #uzInt = Interpolated(xInt,uz,P)
        phiInt = Interpolated(xInt,phi,P)

        duxInt = Interpolated(xInt,ux,dP)
        duzInt = Interpolated(xInt,uz,dP)
        #dphiInt = Interpolated(xInt,phi,dP)

        #vxInt = Interpolated(xInt,vx,P)
        #vzInt = Interpolated(xInt,vz,P)
        OmgInt = Interpolated(xInt,Omg,P)
        
        dvxInt = Interpolated(xInt,vx,dP)
        dvzInt = Interpolated(xInt,vz,dP)
        #dOmgInt = Interpolated(xInt,Omg,dP)


        #Popravi!! Vrjetno mam vse količine zgoraj.
        #Preveri k_ph_Om !! Zadnji C[2] ?!?!
        epsInt = (cos(phi0) .+duxInt).*cos.(phiInt) + (sin(phi0) .+duzInt).*sin.(phiInt) .-1
        gamInt = -(cos(phi0) .+duxInt).*sin.(phiInt) + (sin(phi0) .+duzInt).*cos.(phiInt)


        k_x_Om::Array{Float64} =    (2*(C[2]-C[1])*(duzInt + dt/2*dvzInt .+ sin(phi0)) 
                                    + 4*(C[1]*epsInt.*sin.(phiInt + dt/2*OmgInt) - C[2]*gamInt.*cos.(phiInt+dt/2*OmgInt)) 
                                    + (C[1]+C[2])*( 2*(1+dt)*sin.(-phi0.+2*phiInt+dt*OmgInt) + 2*dt*OmgInt.*cos.(phi0.+2*phiInt+dt*OmgInt) + (2*duxInt +3*dt*dvxInt -dt*OmgInt.*(2*duzInt + dt*dvzInt)).*sin.(2*phiInt+dt*OmgInt) - (2*duzInt +3*dt*dvzInt - dt*OmgInt.*(2*duxInt + dt*dvxInt)).*cos.(2*phiInt+dt*OmgInt))) 
        k_x_dvx::Array{Float64} =   ( 2(C[2]-C[1]) .+ (C[1]+C[2])*( -2*cos.(2*phiInt +dt*OmgInt) + dt*OmgInt.*sin.(2*phiInt +dt*OmgInt)))
        k_x_dvz::Array{Float64} =   ( dt*(C[2]-C[1])*OmgInt - (C[1]+C[2])*( 2*sin.(2*phiInt +dt*OmgInt) + dt*OmgInt.*cos.(2*phiInt +dt*OmgInt)))

        k_z_Om::Array{Float64} =    (2*(C[2]-C[1])*(duxInt + dt/2*dvxInt .+ cos(phi0)) 
                                    + 4*(C[1]*epsInt.*cos.(phiInt + dt/2*OmgInt) + C[2]*gamInt.*sin.(phiInt+dt/2*OmgInt)) 
                                    + (C[1]+C[2])*( 2*(1+dt)*cos.(phi0.+2*phiInt+dt*OmgInt) - 2*dt*OmgInt.*sin.(phi0.+2*phiInt+dt*OmgInt) + (2*duxInt +3*dt*dvxInt -dt*OmgInt.*(2*duzInt + dt*dvzInt)).*cos.(2*phiInt+dt*OmgInt) + (2*duzInt +3*dt*dvzInt -dt*OmgInt.*(2*duxInt + dt*dvxInt)).*sin.(2*phiInt+dt*OmgInt))) 
        k_z_dvx::Array{Float64} =   ( dt*(C[2]+C[1])*OmgInt.*cos.(2*phiInt +dt*OmgInt) +    (C[2]+C[1])*2*sin.(2*phiInt +dt*OmgInt) + dt*OmgInt*(C[2]-C[1]))
        k_z_dvz::Array{Float64} =   ( 2(C[2]-C[1]) .+ (C[1]+C[2])*( -2*cos.(2*phiInt +dt*OmgInt) + dt*OmgInt.*sin.(2*phiInt +dt*OmgInt)))     #Možno da ima napako

        k_ph_Om::Array{Float64}  =  (4*((2*duzInt+dt*dvzInt).*(C[2]*gamInt.*cos.(phiInt+dt*OmgInt)-C[1]*epsInt.*sin.(phiInt+dt/2*OmgInt))+(2 .+2*duxInt+dt*dvxInt).*(C[1]*epsInt.*cos.(phiInt+dt*OmgInt)+C[2]*gamInt.*sin.(phiInt+dt/2*OmgInt)))
                                    +(C[1]+C[2])*(-2*(2*duzInt+dt*dvzInt).*(dt*OmgInt.*cos.(phi0.+ 2*phiInt +dt*OmgInt)+(1+dt)*sin.(-phi0.+2*phiInt+dt*OmgInt)) - 2*(2 .+2*duxInt+dt*dvxInt).*(dt*OmgInt.*sin.(phi0.+ 2*phiInt +dt*OmgInt)-(1+dt)*cos.(-phi0.+2*phiInt+dt*OmgInt)) + cos.(2*phiInt +dt*OmgInt).*(4*(duxInt + duxInt.^2 +duzInt.^2) + 2*dt*(3*dvxInt+4*duxInt.*dvxInt+4*duzInt.*dvzInt-4*OmgInt.*duzInt-8*OmgInt.*duxInt.*duzInt)+ dt^2 *(3*dvxInt.^2 + 3*dvzInt.^2-4*OmgInt.*duzInt.*dvxInt-2*dvzInt.*OmgInt-4*duxInt.*dvzInt.*OmgInt) -dt^3*dvxInt.*dvzInt.*OmgInt)+sin.(2*phiInt+dt*OmgInt).*(4*duzInt+dt*(-4*duzInt.*dvxInt+6*dvzInt+4*duxInt.*dvzInt-4*duxInt.*OmgInt) +dt^2*(-4*duxInt+4*duzInt-2*dvxInt-4*duxInt.*dvxInt+4*duzInt.*dvzInt).*OmgInt+dt^3*(dvzInt.^2-dvxInt.^2).*OmgInt))#             (6. .+4*duxInt).*dvzInt*dt + 4*dvzInt.^2 *dt.*OmgInt + dt^3 *dvzInt.^2 .*OmgInt - dt*(2*duxInt +dt*dvxInt).*(2. .+ 2*duxInt +dt*dvxInt).*OmgInt + duzInt.*(4. .- 4*dvxInt*dt +4*dt^2 *dvzInt.*OmgInt)))
                                    +(C[2]-C[1])*(4*(duxInt+duxInt.^2+dvzInt.^2)+2*dt*(dvxInt+2*duxInt.*dvxInt-2*duzInt.*dvzInt)+dt^2 *(dvxInt.^2-dvzInt.^2)+2*(2. .+2*duxInt+dt*dvxInt)*cos(phi0) -2*(2*duzInt+dt*dvzInt)*sin(phi0)))/2
        k_ph_dOm::Float64 =  -4*C[3]
        k_ph_dvx::Array{Float64} =  (4*C[1]*epsInt.*sin.(phiInt + dt/2*OmgInt)-4*C[2]*gamInt.*cos.(phiInt+dt/2*OmgInt)
                                    +(C[1]-C[2])*(dt*(sin(phi0).-OmgInt*cos(phi0)) + 2*(duzInt +dt*dvzInt) -dt*(1. .+2*duxInt+dt*dvxInt).*OmgInt)
                                    +(C[1]+C[2])*((2*duzInt+dt*(1. .+ 2*duxInt+dt*dvxInt)).*cos.(2*phiInt +dt*OmgInt) + (2*(1. .+duxInt+dt*dvxInt) -dt*(2*duzInt +dt*dvzInt)).*sin.(2*phiInt+dt*OmgInt)+dt*(sin.(-phi0 .+2*phiInt +dt*OmgInt) +OmgInt.*cos.(phi0 .+2*phiInt +dt*OmgInt))))/2
        k_ph_dvz::Array{Float64} =  (4*C[1]*epsInt.*cos.(phiInt + dt/2*OmgInt)+4*C[2]*gamInt.*sin.(phiInt+dt/2*OmgInt)
                                    +(C[1]-C[2])*(2*(1. .+duxInt +dt*dvxInt) +dt*(cos(phi0) .- OmgInt*sin(phi0)))#   (-4 .-4*duxInt -2*dt*(dvxInt+dvzInt) +dt*OmgInt.*(2*duxInt-2*duzInt+dt*(dvxInt-dvzInt))+2*dt*(OmgInt*cos(phi0) .- sin(phi0)))
                                    +(C[1]+C[2])*(-2*(1+dt)*cos.(2*phiInt+dt*OmgInt) +(2*duzInt +2*dvzInt +dt*OmgInt).*sin.(2*phiInt+dt*OmgInt) + dt*(cos.(-phi0 .+2*phiInt+dt*OmgInt)-sin.(phi0 .+2*phiInt+dt*OmgInt)) ))/2   #((-4 .-4*duxInt -2*dt*(dvxInt+dvzInt).*(1 .-dt/2*OmgInt)+2*dt*OmgInt.*(duxInt+duzInt)).*cos.(2*phiInt+dt*OmgInt) + (4*duzInt+2*dt*(dvxInt+dvzInt)+dt^2*OmgInt.*(dvxInt-dvzInt) + 2*dt*OmgInt.*(1 .+duxInt-duzInt)).*sin.(2*phiInt+dt*OmgInt) + 2*dt*(sin.(-phi0 .+2*phiInt+dt*OmgInt)+OmgInt.*cos.(phi0 .+2*phiInt+dt*OmgInt))))/2
        
        k = dt.*[k_x_Om, k_x_dvx, k_x_dvz, k_z_Om,k_z_dvx,k_z_dvz,k_ph_Om,k_ph_dOm,k_ph_dvx,k_ph_dvz]./8
        
        return k
    end
    function JacKoeff(C::SymMat,ux::Array{Float64},uz::Array{Float64},phi::Array{Float64},vx::Array{Float64},vz::Array{Float64},Omg::Array{Float64},phi0::Float64,xInt::Array{Float64},dt::Float64,P::Array{Float64},dP::Array{Float64})
        k = JacKoeff(C.Cd::DiaMat,ux::Array{Float64},uz::Array{Float64},phi::Array{Float64},vx::Array{Float64},vz::Array{Float64},Omg::Array{Float64},phi0::Float64,xInt::Array{Float64},dt::Float64,P::Array{Float64},dP::Array{Float64})

        g = 1

        return k,g
    end
    function JacKoeff(C::AnyMat,ux::Array{Float64},uz::Array{Float64},phi::Array{Float64},vx::Array{Float64},vz::Array{Float64},Omg::Array{Float64},phi0::Float64,xInt::Array{Float64},dt::Float64,P::Array{Float64},dP::Array{Float64})
        k,g = JacKoeff(C.Cs::SymMat,ux::Array{Float64},uz::Array{Float64},phi::Array{Float64},vx::Array{Float64},vz::Array{Float64},Omg::Array{Float64},phi0::Float64,xInt::Array{Float64},dt::Float64,P::Array{Float64},dP::Array{Float64})
        #k = JacKoeff(C.Cs.Cd::DiaMat,ux::Array{Float64},uz::Array{Float64},phi::Array{Float64},vx::Array{Float64},vz::Array{Float64},Omg::Array{Float64},phi0::Float64,xInt::Array{Float64},dt::Float64,P::Array{Float64},dP::Array{Float64})

        c=1

        return k,g,c
    end



    #Funkcija za račun ravnotežnih enačb... Želimo: fxi = 0, fzi = 0, fpi = 0
    function BalanceEq(C::Union{DiaMat,SymMat,AnyMat},M::Array{Float64},g::Float64,ux::Array{Float64},uz::Array{Float64},phi::Array{Float64},vx::Array{Float64},vz::Array{Float64},Omg::Array{Float64},phi0::Float64,xInt::Array{Float64},wInt::Array{Float64},dt::Float64,P::Vector{Array{Float64}},dP::Vector{Array{Float64}},px::Vector{Float64},pz::Vector{Float64},my::Vector{Float64})
        if C isa DiaMat
            C = C.Cd
        elseif C isa SymMat
            C = C.Cs + C.Cd.Cd
        else
            C = C.Ca + C.Cs.Cs + C.Cs.Cd.Cd
        end
        
        #lahko bi vx, vz, Omg vnesel kot dva zaporedna stolpca in potem vzel za pospešek (vx[2]-vx[1])/dt ali nekaj podobnega.

        xInt_ = vcat(-1,xInt,1)

        #uxInt = Interpolated(xInt_,ux,P)
        #uzInt = Interpolated(xInt_,uz,P)
        Omg0=Omg[:,1]
        Omg = Omg[:,2]
        vx0 = vx[:,1]
        vx = vx[:,2]
        vz0 = vz[:,1]
        vz = vz[:,2]

        phiInt = Interpolated(xInt_,phi,P)

        duxInt = Interpolated(xInt_,ux,dP)
        duzInt = Interpolated(xInt_,uz,dP)
        dphiInt = Interpolated(xInt_,phi,dP)

        vxInt = Interpolated(xInt_,vx,P)
        vzInt = Interpolated(xInt_,vz,P)
        OmgInt = Interpolated(xInt_,Omg,P)
        
        dvxInt = Interpolated(xInt_,vx,dP)
        dvzInt = Interpolated(xInt_,vz,dP)
        dOmgInt = Interpolated(xInt_,Omg,dP)

        epsInt = (cos(phi0) .+duxInt).*cos.(phiInt) + (sin(phi0) .+duzInt).*sin.(phiInt) .-1
        gamInt = -(cos(phi0) .+duxInt).*sin.(phiInt) + (sin(phi0) .+duzInt).*cos.(phiInt)
        #kapInt = dphiInt

        tepsInt = sin.(phiInt+dt/2*OmgInt).*(sin.(phi0 .+ dvzInt - OmgInt.*(cos(phi0) .+ duxInt +dt/2*dvxInt))) + cos.(phiInt+dt/2*OmgInt).*(cos(phi0) .+ dvxInt + OmgInt.*(sin(phi0) .+ duzInt +dt/2*dvzInt))
        tgamInt = cos.(phiInt+dt/2*OmgInt).*(sin.(phi0 .+ dvzInt - OmgInt.*(cos(phi0) .+ duxInt +dt/2*dvxInt))) - sin.(phiInt+dt/2*OmgInt).*(cos(phi0) .+ dvxInt + OmgInt.*(sin(phi0) .+ duzInt +dt/2*dvzInt))
        #tkapInt = dOmgInt



        Rx =    (C[1,1]*(epsInt +dt/2*tepsInt) +C[1,2]*(gamInt +dt/2*tgamInt)+C[1,3]*(dphiInt +dt/2*dOmgInt)).*cos.(phiInt + dt/2*OmgInt)
                +(C[2,2]*(gamInt +dt/2*tgamInt)+C[2,1]*(epsInt +dt/2*tepsInt) + C[2,3]*(dphiInt +dt/2*dOmgInt)).*sin.(phiInt + dt/2*OmgInt)
        Rz =    -(C[1,1]*(epsInt +dt/2*tepsInt)+C[1,2]*(gamInt +dt/2*tgamInt)+C[1,3]*(dphiInt +dt/2*dOmgInt)).*sin.(phiInt + dt/2*OmgInt)
                +(C[2,2]*(gamInt +dt/2*tgamInt)+C[2,1]*(epsInt +dt/2*tepsInt) + C[2,3]*(dphiInt +dt/2*dOmgInt)).*cos.(phiInt + dt/2*OmgInt)
        Mc =    C[3,3]*(dphiInt + dt/2*dOmgInt) + C[3,1]*(epsInt + dt/2*tepsInt) + C[3,2]*(gamInt + dt/2*tgamInt)


        Rx0 = Rx[[1,end]]
        Rz0 = Rz[[1,end]]
        Mc0 = Mc[[1,end]]


        Rx = Rx[2:end-1]
        Rz = Rz[2:end-1]
        Mc = Mc[2:end-1]

        


        #Popravi pospeške... Predznak, gravitacija,...
        fxi = map( i -> wInt'*( -Rx.*Interpolated(xInt,1., dP[i]) +  ((px[1]+px[2] .+xInt*(px[2]-px[1]))/2 - (vxInt[2:end-1] -Interpolated(xInt,vx0,P))/dt*M[1])  .* Interpolated(xInt,1.,P[i]) )  + Rx0' *Interpolated([-1.,1.],1.,P[i])  , eachindex(P))
        fzi = map( i -> wInt'*( -Rz.*Interpolated(xInt,1., dP[i]) +  ((pz[1]+pz[2] .+xInt*(pz[2]-pz[1]))/2 - ((vzInt[2:end-1]-Interpolated(xInt,vz0,P))/dt.-g)*M[1])  .*Interpolated(xInt,1.,P[i]) )  + Rz0' *Interpolated([-1.,1.],1.,P[i]) , eachindex(P))
        fpi = map( i -> wInt'*( -Mc.*Interpolated(xInt,1., dP[i]) +  (Rx.*(duzInt[2:end-1]+dt/2*dvzInt[2:end-1]) -Rz.*(1 .+duxInt[2:end-1]+dt/2*dvzInt[2:end-1])+  (my[1]+my[2] .+xInt*(my[2]-my[1]))/2 -  M[2]*(OmgInt[2:end-1]-Interpolated(xInt,Omg0,P))/dt  ).*Interpolated(xInt,1.,P[i]) )  +Mc0' *Interpolated([-1.,1.],1.,P[i]) , eachindex(P))

        #round.(fxi,digits = 14)
        #round.(fzi,digits = 14)
        #round.(fpi,digits = 14)

        return fxi,fzi,fpi
    end


    struct Polynomial
        ai :: Array{<:Real}
    end

end