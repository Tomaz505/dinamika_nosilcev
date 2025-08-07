#   B R A N J E   P O D A T K O V   I Z   D A T O T E K E
data1,data2,data3 = readdata();
eval(Meta.parse(data1));
eval(Meta.parse(data2));
n_elem,n_voz,ElementDataIn,VozDataIn = datainit(elementi,vozlisca);
eval.(Meta.parse.(data3));
println("\n[  Ok  ]  Vnos podatkov")



#   P R E D P O C E S I R A N J E
E = Array{BeamDataProcess}(undef,n_elem)
n_nodes = n_voz
for i = 1:n_elem
	global n_nodes 
	E[i], n_nodes = dataprocess(ElementDataIn[i],VozDataIn[ElementDataIn[i].v],n_nodes)
end
println("[  Ok  ]  Procesiranje podatkov")



#   P O D A T K I   I T R A C I J S K E G A   P O S T O P K A
# Tole spravi v data.txt
ti = 0.
tf = 1.
dt = 10.0^-3
g = 0.


t = collect(ti:dt:tf)
n_time = length(t)

#Komponente za razvoj interpolacije (vključuje odvod če je interpoliran)
#n_nodes = n_voz + sum(map(i -> .+(size.(E[i].P)...)[1]-2-(length(ElementDataIn[i].div2)-1)*(1+Int(ElementDataIn[i].Ci)),1:n_elem))

M =  BeamMotion(zeros(n_nodes,n_time),zeros(n_nodes,n_time),zeros(n_nodes,n_time),zeros(n_nodes,n_time),zeros(n_nodes,n_time),zeros(n_nodes,n_time))


indx = setdiff(vcat(map(j->(VozDataIn[j].i .+[0;n_nodes;2*n_nodes]).*(-VozDataIn[j].Supp'.+1), 1:n_voz )...),[0])
indx = setdiff(1:3*n_nodes,indx)

indxX = indx[findall(indx.<= n_nodes)]
indxZ = indx[findall(n_nodes+1 .<= indx .<= 2*n_nodes)].-n_nodes
indxP = indx[findall(2*n_nodes+1 .<=indx.<= 3*n_nodes)].-2*n_nodes

# Tangentna in rezidual v eni matriki
JR=spzeros(Float64,3*n_nodes,3*n_nodes+1)
Dv=zeros(Float64,3*n_nodes,1)


println("[  Ok  ]  Priprava na račun")


i_el = 1
i_time = 2
i_ke = 1
i_indx = 1:4

F = Tan_Res(E[i_el].xInt[i_ke],E[i_el].wInt[i_ke],M.ux[i_indx,[i_time,i_time-1]],M.uz[i_indx,[i_time,i_time-1]],M.phi[i_indx,[i_time,i_time-1]],M.vx[i_indx,[i_time,i_time-1]],M.vz[i_indx,[i_time,i_time-1]],M.Omg[i_indx,[i_time,i_time-1]],E[i_el].P[i_ke],E[i_el].p0[i_ke],E[i_el].k0[i_ke],ElementDataIn[i_el].C,ElementDataIn[i_el].M[:,i_ke],ElementDataIn[i_el].px(t[i_time])[i_ke,:],ElementDataIn[i_el].pz(t[i_time])[i_ke,:],ElementDataIn[i_el].my(t[i_time])[i_ke,:],dt,E[i_el].pb[i_ke],E[i_el].kb[i_ke],E[i_el].L[i_ke],g)

#=
for i_time in 2:n_time
	while norm(Dv) > 10.0^-6
		JR[:,:] .= 0
      		for i_el in eachindex(E)
               		for i_ke in eachindex(E[j].P)
               		end # i_ke
              	end # i_el
		Dv = J\R
		v[indx] -= Dv
		u[indx] = u[indx-1] + dt*v[indx]
       	end # while norm(Re)
	v[indx+1] = 2.0*v[indx] - v[indx-1]
end # i_time
=#
println("[  Ok  ]  Račun")

