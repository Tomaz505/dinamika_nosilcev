#   B R A N J E   P O D A T K O V   I Z   D A T O T E K E

data1,data2,data3 = readdata();
eval(Meta.parse(data1));
eval(Meta.parse(data2));
n_elem,n_voz,ElementDataIn,VozDataIn = datainit(elementi,vozlisca);
eval.(Meta.parse.(data3));



#   P R E D P O C E S I R A N J E
E = Array{BeamDataProcess}(undef,n_elem)
for i = 1:n_elem
	E[i] = dataprocess(ElementDataIn[i],VozDataIn[ElementDataIn[i].v])
end





#=







#=   P R I P R A V A   P O D A T K O V
      Podatki se shranijo v seznamu elementov tipa BeamDataProcess

      E[i].A1[j][k] -seznam za i-ti element, z vrednostmi produkta interpolacijskih polinomov j in k. A2,A3 ... vsebujejo produkte odvodov. 
=#

begin
    E = Array{BeamDataProcess}(undef,n_elem,1)
    for i in eachindex(ElementDataIn)
        el = ElementDataIn[i]


        g = el.f.(el.div1)
        g = map(k-> k-g[1],g)
        g = g/norm(g[end])
        g = round.(map(k->(g[k][1] .+ im*g[k][2])/(g[end][1] .+ im*g[end][2]),eachindex(g)),digits=14)
        g =  map(k-> g[k]-g[k-1],setdiff(eachindex(g),[1]))
        ϕ0 = angle.(g)
        L = abs.(g)
        z = im*(VozDataIn[el.v[2]].y-VozDataIn[el.v[1]].y) + (VozDataIn[el.v[2]].x-VozDataIn[el.v[1]].x)
        ϕ0 = ϕ0.+ angle(z)
        L = abs(z)*L

        P = InterpolKoeff.(map(k->range(-1,1,length=k),el.div2))
        dP = vcat(map(k->[P[k][2]],eachindex(P))...)
        P = vcat(map(k->[P[k][1]],eachindex(P))...)
        xInt,wInt = GaussInt(el.nInt::Array{Int64})

        A1= Array{Array{Array{Float64}}}(undef,(length(P),1))
        A2 = copy(A1)
        A3 = copy(A1)
        A4 = copy(A1)
        A5 = copy(A1)
        A6 = copy(A1)
        for j in eachindex(P)
            A1[j] = L[j]/2*map(k-> ( P[j][k[1]]' * (xInt[j]' .^(0:length( P[j][k[1]])-1)))' .* ( dP[j][k[2]]' * (xInt[j]' .^(0:length( dP[j][k[2]])-1)))' .* wInt[j],   CartesianIndex.(eachindex(P[j]),eachindex(P[j])'))
            A2[j] = L[j]/2*map(k-> (dP[j][k[1]]' * (xInt[j]' .^(0:length(dP[j][k[1]])-1)))' .* ( P[j][k[2]]' * (xInt[j]' .^(0:length( P[j][k[2]])-1)))' .* wInt[j],   CartesianIndex.(eachindex(P[j]),eachindex(P[j])'))
            A3[j] = L[j]/2*map(k-> (dP[j][k[1]]' * (xInt[j]' .^(0:length(dP[j][k[1]])-1)))' .* (dP[j][k[2]]' * (xInt[j]' .^(0:length(dP[j][k[2]])-1)))' .* wInt[j],   CartesianIndex.(eachindex(P[j]),eachindex(P[j])'))
            A4[j] = L[j]/2*map(k-> ( P[j][k[1]]' * (xInt[j]' .^(0:length( P[j][k[1]])-1)))' .* ( P[j][k[2]]' * (xInt[j]' .^(0:length( P[j][k[2]])-1)))' .* wInt[j],   CartesianIndex.(eachindex(P[j]),eachindex(P[j])'))
            A5[j] = L[j]/2*map(k-> ( P[j][k]' * (xInt[j]' .^(0:length( P[j][k])-1)))' .* wInt[j],   eachindex(P[j]))
            A6[j] = L[j]/2*map(k-> ( dP[j][k]' * (xInt[j]' .^(0:length( dP[j][k])-1)))' .* wInt[j],   eachindex(P[j]))

        end

        #k = sum(el.div2[1:length(el.div1)-1])-length(el.div1)+2 
        v = map(j-> collect(0:j-1),el.div2)
        k = sum(el.div2[1:length(el.div1)-1])-length(el.div1)
        for i in eachindex(v)
            v[i] = v[i].+(n_voz)
            n_voz = n_voz+length(v[i])-1
        end
        n_voz = Int(length(vozlisca)/2)+k
        v[1][1]  = el.v[1]
        v[end][end] = el.v[2]
        

        E[i] = BeamDataProcess(ElementDataIn[i].C,el.M,L,ϕ0,P,dP,A1,A2,A3,A4,A5,A6,xInt,wInt,v)
    end
    n_voz = Int(length(vozlisca)/2)
    n_v = maximum(E[end].v[end])
end;
=#






#=

#   P O D A T K I   I T R A C I J S K E G A   P O S T O P K A
begin
    ti = 0.
    tf = 1.
    dt = 10^-3
    
    g = 0.

    #   METODA
    begin

        t = collect(ti:dt:tf)
        j = length(t)
        global M =  BeamMotion(zeros(n_v,j),zeros(n_v,j),zeros(n_v,j),zeros(n_v,j),zeros(n_v,j),zeros(n_v,j) )


        indx = setdiff(vcat(map(j->(VozDataIn[j].i .+[0;n_v;2*n_v]).*(-VozDataIn[j].Supp'.+1), 1:n_voz )...),[0])
        indx = setdiff(1:3*n_v,indx)

        indxX = indx[findall(indx.<= n_v)]
        indxZ = indx[findall(n_v+1 .<= indx .<= 2*n_v)].-n_v
        indxP = indx[findall(2*n_v+1 .<=indx.<= 3*n_v)].-2*n_v

        #J=spzeros(3*n_v,3*n_v)
        #R=zeros(3*n_v,1)

        JR=spzeros(3*n_v,3*n_v+1)
        Dv=zeros(3*n_v,1)

        for i=1:1 # in eachindex(t)
            for iter=1:1
                #display(iter)

                JR[:,:] .=0

                for j in eachindex(E)
                    for l in eachindex(E[j].P)
                        #=
                            vx = (M.vx[E[j].v[l],i] + M.vx[E[j].v[l],i-1])/2
                            vz = (M.vx[E[j].v[l],i] + M.vx[E[j].v[l],i-1])/2
                            Ω = (M.vx[E[j].v[l],i] + M.vx[E[j].v[l],i-1])/2


                            
                            k = JacKoeff(E[j].C, M.ux[E[j].v[l],i] ,M.uz[E[j].v[l],i],M.phi[E[j].v[l],i],vx,vz,Ω,E[j].ϕ0[l],E[j].xInt[l],dt,E[j].P[l],E[j].dP[l])
                            fxi,fzi,fpi = BalanceEq(E[j].C,E[j].M,g,M.ux[E[j].v[l],i] ,M.uz[E[j].v[l],i],M.phi[E[j].v[l],i],M.vx[E[j].v[l],i:i+1],M.vz[E[j].v[l],i:i+1],M.Omg[E[j].v[l],i:i+1],E[j].ϕ0[l],E[j].xInt[l],E[j].wInt[l],dt,E[j].P[l],E[j].dP[l],ElementDataIn[j].px(t[i])[l,:],ElementDataIn[j].pz(t[i])[l,:],ElementDataIn[j].my(t[i])[l,:])

                            dfxOmg = map(a->k[1]'* E[j].A2[l][a],   CartesianIndex.(1:length(E[j].P[l]),(1:length(E[j].P[l]))'))
                            dfxdvx = map(a->k[2]'* E[j].A3[l][a],   CartesianIndex.(1:length(E[j].P[l]),(1:length(E[j].P[l]))'))
                            dfxdvz = map(a->k[3]'* E[j].A3[l][a],   CartesianIndex.(1:length(E[j].P[l]),(1:length(E[j].P[l]))'))
                            
                            dfzOmg = map(a->k[4]'* E[j].A2[l][a],   CartesianIndex.(1:length(E[j].P[l]),(1:length(E[j].P[l]))'))
                            dfzdvx = map(a->k[5]'* E[j].A3[l][a],   CartesianIndex.(1:length(E[j].P[l]),(1:length(E[j].P[l]))'))
                            dfzdvz = map(a->k[6]'* E[j].A3[l][a],   CartesianIndex.(1:length(E[j].P[l]),(1:length(E[j].P[l]))'))
                            
                            dfpOmg = map(a-> k[7]'* E[j].A2[l][a] + sum(k[8]* E[j].A3[l][a]),   CartesianIndex.(1:length(E[j].P[l]),(1:length(E[j].P[l]))'))
                            dfpdvx = map(a-> k[9]'* E[j].A2[l][a[2],a[1]],   CartesianIndex.(1:length(E[j].P[l]),(1:length(E[j].P[l]))'))
                            dfpdvz = map(a-> k[10]'* E[j].A2[l][a[2],a[1]],   CartesianIndex.(1:length(E[j].P[l]),(1:length(E[j].P[l]))'))
                            
                            J[E[j].v[l] , E[j].v[l] ] = dfxdvx
                            J[E[j].v[l] , E[j].v[l].+n_v ] = dfxdvz
                            J[E[j].v[l] , E[j].v[l].+n_v*2 ] = dfxOmg

                            J[E[j].v[l].+n_v , E[j].v[l] ] = dfzdvx
                            J[E[j].v[l].+n_v , E[j].v[l].+n_v ] = dfzdvz
                            J[E[j].v[l].+n_v , E[j].v[l].+n_v*2 ] = dfzOmg

                            J[E[j].v[l].+n_v*2 , E[j].v[l] ] = dfpdvx
                            J[E[j].v[l].+n_v*2 , E[j].v[l].+n_v ] = dfpdvz
                            J[E[j].v[l].+n_v*2 , E[j].v[l].+n_v*2 ] = dfpOmg

                            R[E[j].v[l]] = fxi*E[j].L[l]/2
                            R[E[j].v[l].+n_v] = fzi*E[j].L[l]/2
                            R[E[j].v[l].+n_v*2] = fpi*E[j].L[l]/2
                        =#

                        #println(vcat([E[j].v[l],E[j].v[l].+n_v,E[j].v[l].+n_v*2]...),vcat([E[j].v[l],E[j].v[l].+n_v,E[j].v[l].+n_v*2,n_v*3+1]...))
                        #V enem kosu prišteje komponente Jakobija in desnih strani da omogoči račun z več Threadsi
                        #JR[vcat([E[j].v[l],E[j].v[l].+n_v,E[j].v[l].+n_v*2]...),vcat([E[j].v[l],E[j].v[l].+n_v,E[j].v[l].+n_v*2,n_v*3+1]...)] += 
                        evaluateKoefficient(M.ux[E[j].v[l],i] ,M.uz[E[j].v[l],i], M.phi[E[j].v[l],i], (M.vx[E[j].v[l],i+1] + M.vx[E[j].v[l],i])/2, (M.vz[E[j].v[l],i+1] + M.vz[E[j].v[l],i])/2, (M.Omg[E[j].v[l],i+1] + M.Omg[E[j].v[l],i])/2, E[j].ϕ0[l], E[j].L[l],E[j].C,ElementDataIn[j].px(t[i])[l,:],ElementDataIn[j].pz(t[i])[l,:],ElementDataIn[j].my(t[i])[l,:],g,E[j].xInt[l],E[j].wInt[l],dt,E[j].P[l],E[j].dP[l],E[j].A1[l],E[j].A2[l],E[j].A3[l],E[j].A4[l],E[j].A5[l],E[j].A6[l])
                        #display(vcat([E[j].v[l],E[j].v[l].+n_v,E[j].v[l].+n_v*2]))
                    end
                end #Račun komponent jakobija in desnih strani
                #display(JR)
                #Dv =  -J[indx,indx]\R[indx]
                #Dv[indx] = -JR[indx,indx]\Vector(JR[indx,end])
                
                #M.vx[indxX,i+1] += Dv[indxX]
                #M.vz[indxZ,i+1] += Dv[indxZ.+n_v] 
                #M.Omg[indxP,i+1] += Dv[indxP.+n_v*2] 

                #println("i = $iter, ||Dv|| = $(norm(Dv))")
            end #while

            #M.ux[indxX,i+1] += M.ux[indxX,i] + dt*(M.vx[indxX,i+1]+M.vx[indxX,i])/2
            #M.uz[indxZ,i+1] += M.uz[indxZ,i] + dt*(M.vz[indxZ,i+1]+M.vz[indxZ,i])/2
            #M.phi[indxP,i+1] += M.phi[indxP,i] + dt*(M.Omg[indxP,i+1]+M.Omg[indxP,i])/2

            #M.vx[indxX,i+2] = M.vx[indxX,i+1]
            #M.vz[indxZ,i+2] = M.vz[indxX,i+1] 
            #M.Omg[indxP,i+2] = M.Omg[indxX,i+1] 
        end #for i in t
    end
end=#

