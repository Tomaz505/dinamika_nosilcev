#   B R A N J E   P O D A T K O V   I Z   D A T O T E K E


println("Pot do datoteke z podatki")
file = readline()
include(file*".jl")
#data1,data2,data3 = readdata(file);
#eval(Meta.parse(data1));
#eval(Meta.parse(data2));
#n_elem,n_voz,ElementDataIn,VozDataIn = datainit(elementi,vozlisca);
#eval.(Meta.parse.(data3));
println("\n[  Ok  ]  Vnos podatkov")









#   P R E D P O C E S I R A N J E
E = Array{BeamDataProcess}(undef,n_elem)
n_nodes = n_voz
for i = 1:n_elem
	global n_nodes 
	E[i], n_nodes = dataprocess(ElementDataIn[i],VozDataIn[ElementDataIn[i].v],n_nodes)
end
println("[  Ok  ]  Procesiranje podatkov")









#   P L O T
plotbeams(E,ElementDataIn,VozDataIn)
println("[  Ok  ]  Risnaje Konstrikcije")



println("\n Kako nadaljujem?")
println("0\t\t-> preklici postopek")
println("Enter\t\t-> nadaljuj račun")

canc = readline()
if canc == "0"
	error("Preklic")
end










#   P O D A T K I   I T R A C I J S K E G A   P O S T O P K A
# Tole spravi v data.txt
#

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

indx_solve = sort(vcat(3*indxX.-2,3*indxZ.-1,3*indxP))
println("[  Ok  ]  Priprava na račun")
#=
i_time = 2
i_el=1
i_ke=1

dlF,F = Tan_Res(E[i_el].xInt[i_ke],E[i_el].wInt[i_ke],M.ux[E[i_el].indx[i_ke],[i_time-1,i_time]],M.uz[E[i_el].indx[i_ke],[i_time-1,i_time]],M.phi[E[i_el].indx[i_ke],[i_time-1,i_time]],M.vx[E[i_el].indx[i_ke],[i_time-1,i_time]],M.vz[E[i_el].indx[i_ke],[i_time-1,i_time]],M.Omg[E[i_el].indx[i_ke],[i_time-1,i_time]],E[i_el].P[i_ke],E[i_el].p0[i_ke],E[i_el].k0[i_ke],ElementDataIn[i_el].C,reshape(ElementDataIn[i_el].M[i_ke,:],2),ElementDataIn[i_el].px(t[i_time])[i_ke,:],ElementDataIn[i_el].pz(t[i_time])[i_ke,:],ElementDataIn[i_el].my(t[i_time])[i_ke,:],dt,E[i_el].pb[i_ke],E[i_el].kb[i_ke],E[i_el].L[i_ke],g)
=#
for i_time in 2:n_time
	Dv = [1.;1.]

	println("i_time = ",i_time,"\t\tt = ",t[i_time])
	println()

	while norm(Dv) > 10.0^-7
		# Tangentna in rezidual v eni matriki
		Ja = spzeros(Float64,3*n_nodes,3*n_nodes)
		Re = zeros(Float64,3*n_nodes,1)


			for i_el in eachindex(E)
				for i_ke in eachindex(E[i_el].P)

					dlF,F = Tan_Res(E[i_el].xInt[i_ke],E[i_el].wInt[i_ke],M.ux[E[i_el].indx[i_ke],[i_time-1,i_time]],M.uz[E[i_el].indx[i_ke],[i_time-1,i_time]],M.phi[E[i_el].indx[i_ke],[i_time-1,i_time]],M.vx[E[i_el].indx[i_ke],[i_time-1,i_time]],M.vz[E[i_el].indx[i_ke],[i_time-1,i_time]],M.Omg[E[i_el].indx[i_ke],[i_time-1,i_time]],E[i_el].P[i_ke],E[i_el].p0[i_ke],E[i_el].k0[i_ke],ElementDataIn[i_el].C,reshape(ElementDataIn[i_el].M[i_ke,:],2),ElementDataIn[i_el].px(t[i_time])[i_ke,:],ElementDataIn[i_el].pz(t[i_time])[i_ke,:],ElementDataIn[i_el].my(t[i_time])[i_ke,:],dt,E[i_el].pb[i_ke],E[i_el].kb[i_ke],E[i_el].L[i_ke],g)


					#indx_inv_perm =  vcat(map(i->findall(E[i_el].indx[i_ke].==i),sort(E[i_el].indx[i_ke]))...)
					#indx_inv = ((E[i_el].indx[i_ke])[indx_inv_perm])[indx_inv_perm]


					indx_dof = reshape(hcat(E[i_el].indx[i_ke]*3 .-2,E[i_el].indx[i_ke]*3 .-1, E[i_el].indx[i_ke]*3)',(3*length(E[i_el].indx[i_ke])))
					#indx_dof = vcat(E[i_el].indx[i_ke]*3 .-2,E[i_el].indx[i_ke]*3 .-1, E[i_el].indx[i_ke]*3)

					Ja[indx_dof,indx_dof] +=hvcat(length(F),dlF...)
					Re[indx_dof] += vcat(F...)
				end # i_ke
			end # i_el
		# inverzni indeksi
		# c = vcat(map(i->findall(sort(arr) .== i), arr)...)
		# inv_arr_indx = invpermute(arr,c)
		#
		Dv = round.(-Ja[indx_solve,indx_solve]\Re[indx_solve],digits = 10)

		# *2 je zaradi povprečenja za vmesno ??? je to ok ???
		M.vx[indxX,i_time] += Dv[1:3:end]*2
		M.vz[indxZ,i_time] += Dv[2:3:end]*2
		M.Omg[indxP,i_time] += Dv[3:3:end]*2

		println("\e[2K\e[1G\t|Dv| = ",norm(Dv))
	
       	end # while norm(Dv) > x


	M.vx[indxX,i_time+1] = M.vx[indxX,i_time]
	M.vz[indxZ,i_time+1] = M.vz[indxZ,i_time]
	M.Omg[indxP,i_time+1] = M.Omg[indxP,i_time]

	M.ux[indxX,i_time] = M.ux[indxX,i_time-1] + dt*(M.vx[indxX,i_time]+M.vx[indxX,i_time-1])/2.
	M.uz[indxZ,i_time] = M.uz[indxZ,i_time-1] + dt*(M.vz[indxZ,i_time]+M.vz[indxZ,i_time-1])/2.
	M.phi[indxP,i_time] = M.phi[indxP,i_time-1] + dt*(M.Omg[indxP,i_time]+M.Omg[indxP,i_time-1])/2.
		

end # i_time

println("[  Ok  ]  Račun")

