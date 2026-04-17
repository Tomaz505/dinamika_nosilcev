#   B R A N J E   P O D A T K O V   I Z   D A T O T E K E
#
begin
println("Ime datoteke v ./in")
file = readline()
include("in/"*file*".jl")
println()
@info "Vnos podatkov\n\t\t[  Ok  ]"
end


#   P R E D P O C E S I R A N J E
#
begin
E = Array{BeamDataProcess}(undef,n_elem)
n_nodes = n_voz
n_int_nodes = 0
for i = 1:n_elem
	global n_nodes 
	global n_int_nodes
	E[i], n_nodes,n_int_nodes = dataprocess(ElementDataIn[i],VozDataIn[ElementDataIn[i].v],n_nodes,n_int_nodes;mtd = Integracija)
end
@info "Procesiranje podatkov\n\t\t[  Ok  ]"
end


#   P L O T
#
begin
	println("\n Kako nadaljujem?")
	println("0 → prekini postopek")
	println("1 → nariši konstrukcijo in nadaljuj")
	println("2 → nariši konstrukcijo in prekini")
	println("↲ → nadaljuj račun")
	local elt = time()
	local test = readline()

	if (test == "0")
		error("Preklic")
	elseif (test == "1")
		konstr_img = plotbeams(E,ElementDataIn,VozDataIn)
		display(konstr_img)
		@info "Risnaje konstrukcije\n\t\t[  Ok  ]"
	elseif (test == "2")
	konstr_img = plotbeams(E,ElementDataIn,VozDataIn)
	display(konstr_img)
	@info "Risnaje konstrukcije\n\t\t[  Ok  ]"
	error("Preklic")
	end
end





time_st = collect(ti:dt:tf)
n_time = length(time_st)






begin




	#Tip casovne integracije -> Multiple dispatch
	if contains(metoda_t_integracije, "timeelement")

		tInt = if contains(metoda_t_integracije, "T")
				Ibtime = trig_re_gramschmid([tnodes])
				xt,wt = QuadInt(it)
				# !!
				xt = (xt.+1)*dt/2*(nt-1)
				wt = wt*dt/2*(nt-1)
				DTvals = TrigValue(xt,Ibtime;n=1)
				Tvals = TrigValue(xt,Ibtime)
				ITvals = TrigValue(xt,Ibtime;n=-1) .- TrigValue(0.0,Ibtime;n=-1)

				TimeElement(dt,tnodes,Ibtime,DTvals,Tvals,ITvals,xt,wt,nt)
			elseif contains(metoda_t_integracije, "P")
				Ibtime = re_gramschmid([tnodes])
				xt,wt = QuadInt(it)
				# !!
				xt = (xt.+1)*dt/2*(nt-1)
				wt = wt*dt/2*(nt-1)
				DTvals = PolyValue(xt,Ibtime;n=1)
				Tvals = PolyValue(xt,Ibtime)
				ITvals = PolyValue(xt,Ibtime;n=-1)

				TimeElement(dt,tnodes,Ibtime,DTvals,Tvals,ITvals,xt,wt,nt)
			end
	elseif metoda_t_integracije == "midpoint"
		tInt = MidPoint(dt,1)
	end


	#Komponente za razvoj interpolacije (vključuje odvod če je interpoliran)
	#n_nodes = n_voz + sum(map(i -> .+(size.(E[i].P)...)[1]-2-(length(ElementDataIn[i].div2)-1)*(1+Int(ElementDataIn[i].Ci)),1:n_elem))


	indx = setdiff(vcat(map(j->(VozDataIn[j].i .+[0;n_nodes;2*n_nodes]).*(-VozDataIn[j].Supp'.+1), 1:n_voz )...),[0])
	indx = setdiff(1:3*n_nodes,indx)

	indxX = indx[findall(indx.<= n_nodes)]
	indxZ = indx[findall(n_nodes+1 .<= indx .<= 2*n_nodes)].-n_nodes
	indxP = indx[findall(2*n_nodes+1 .<=indx.<= 3*n_nodes)].-2*n_nodes
	for i = eachindex(VozDataIn)
		if VozDataIn[i].Supp[1]==1
			global indxX = [indxX;i]
		end
		if VozDataIn[i].Supp[2]==1
			global indxZ = [indxZ;i]
		end
		if VozDataIn[i].Supp[3]==1
			global indxP = [indxP;i]
		end
	end

	indxX = sort(unique(indxX))
	indxZ = sort(unique(indxZ))
	indxP = sort(unique(indxP))


	if metoda_t_integracije == "timeelement"
		indx_solve = vcat(map(i-> sort(vcat(3*indxX.-2,3*indxZ.-1,3*indxP)) .+n_nodes*(i-1)*3, 2:nt)...)
	elseif metoda_t_integracije == "midpoint"
		indx_solve =  sort(vcat(3*indxX.-2,3*indxZ.-1,3*indxP))
	end

	indxX_solve = map(i->findfirst(indx_solve.==(indxX[i]*3 .-2)),eachindex(indxX))
	indxZ_solve = map(i->findfirst(indx_solve.==(indxZ[i]*3 .-1)),eachindex(indxZ))
	indxP_solve = map(i->findfirst(indx_solve.==(indxP[i]*3)),eachindex(indxP))

	n_nodes2 = n_nodes*3
	n_pnodes = maximum(indxP_solve)
	dp_nodes = 0
	for i1 = eachindex(E)
		if ElementDataIn[i1].relese[1]
			global n_nodes2 += 1
			global n_pnodes +=1
			global dp_nodes +=1
			global indx_solve = [indx_solve; n_nodes2]
			global indxP_solve = [indxP_solve;n_pnodes]
			global E[i1].indxP[1][1] = n_nodes+dp_nodes
		end
		if ElementDataIn[i1].relese[2]
			global n_nodes2 +=1
			global n_pnodes +=1
			global dp_nodes +=1
			global indx_solve = [indx_solve; n_nodes2]
			global indxP_solve = [indxP_solve;n_pnodes]
			global E[i1].indxP[end][end] = n_nodes+dp_nodes
		end
	end



	if metoda_t_integracije == "timeelement"
		local Ja = zeros(Float64,3*n_nodes*nt,3*n_nodes*nt)
		local Re0 = ones(Float64,3*n_nodes*nt,1)
	elseif metoda_t_integracije == "midpoint"
		local Ja = zeros(Float64,n_nodes2,n_nodes2)
		local Re0 = ones(Float64,n_nodes2)
	end

	local Dv::Array{Float64}


	#vrednosti polinomov v integracijskih točkah
	local Pvalues = map(i1-> map(i2 -> PolyValue(E[i1].xInt[i2],E[i1].P[i2]), eachindex(E[i1].L)),eachindex(E))
	local dPvalues = map(i1-> map(i2 -> PolyValue(E[i1].xInt[i2],E[i1].P[i2];n=1), eachindex(E[i1].L)),eachindex(E))

	n_nodes2 = n_nodes*3
	indx_dof =  map(i_el->
				 map(i_ke->
					begin

						local indx1 = E[i_el].indx[i_ke]*3 .-2
						local indx2 = E[i_el].indx[i_ke]*3 .-1
						local indx3 = E[i_el].indx[i_ke]*3

						if length(E[i_el].L) == 1
							if ElementDataIn[i_el].relese[1] && ElementDataIn[i_el].relese[2]
								global n_nodes2 +=2
								indx3 = [n_nodes2-1;indx3[2:end-1];n_nodes2]
							elseif ElementDataIn[i_el].relese[1]
								global n_nodes2 +=1
								indx3 = [n_nodes2;indx3[2:end]]
							elseif ElementDataIn[i_el].relese[2]
								global n_nodes2 +=1
								indx3 = [indx3[1:end-1];n_nodes2]
							end

						else
							if i_ke == 1 && ElementDataIn[i_el].relese[1]
								global n_nodes2 +=1
								indx3 = [n_nodes2;indx3[2:end]]
							elseif i_ke == length(E[i_el].indx) && ElementDataIn[i_el].relese[2]
								global n_nodes2 +=1
								indx3 = [indx3[1:end-1];n_nodes2]
							end
						end
						#display(indx1)
						#display(indx2)
						#display(indx3)

						return reshape(hcat(indx1,indx2,indx3)',3*length(indx1))
					end,
				 eachindex(E[i_el].P)),
				eachindex(E))

	if dp_nodes >0
		push!(indxP,(1:dp_nodes).+maximum(indxP) ...)
	end



	mass_m = copy(Ja)
	for i1 = eachindex(E)
		for i2 = eachindex(E[i1].P)
			mi = Mass(E[i1].xInt[i2],E[i1].wInt[i2],Pvalues[i1][i2],dPvalues[i1][i2],ElementDataIn[i1].M,tInt)
			mass_m[indx_dof[i1][i2],indx_dof[i1][i2]] = mass_m[indx_dof[i1][i2],indx_dof[i1][i2]]+hvcat(size(mi)[1],mi...)
		end
	end

	M =  BeamMotion(zeros(n_nodes,n_time),zeros(n_nodes,n_time),zeros(n_nodes+dp_nodes,n_time),zeros(n_nodes,n_time),zeros(n_nodes,n_time),zeros(n_nodes+dp_nodes,n_time),zeros(n_int_nodes,n_time),zeros(n_int_nodes,n_time),zeros(n_int_nodes,n_time))



	#Vsiljeno gibanje podprte prostostne stopnje
	for i1 in eachindex(VozDataIn)
		local V = (hcat(VozDataIn[i1].mot.(time_st)...)-hcat(VozDataIn[i1].mot.(time_st.-dt)...))/dt
		local U = hcat(VozDataIn[i1].mot.(time_st)...)

		M.vx[i1,:] = V[1,:]
		M.vz[i1,:] = V[2,:]
		M.Omg[i1,:] = V[3,:]
		M.ux[i1,:] = U[1,:]
		M.uz[i1,:] = U[2,:]
		M.phi[i1,:] = U[3,:]

	#=for i2 = 2:length(time_st)
	M.vx[i1,i2] = -M.vx[i1,i2-1]+2*V[1,i2-1]
	M.vz[i1,i2] = -M.vz[i1,i2-1]+2*V[2,i2-1]
	M.Omg[i1,i2] = -M.Omg[i1,i2-1]+2*V[3,i2-1]
	end=#
	end

	@info "Priprava na račun\n\t\t[  Ok  ] "

	@time for i_time in 2:nt-1:n_time
		Dv = [1.;1.]
		Re = Re0

		# Prediktor
		M.vx[indxX,i_time],M.vz[indxZ,i_time],M.Omg[indxP,i_time], M.ux[indxX,i_time],M.uz[indxZ,i_time],M.phi[indxP,i_time] = prediktor(M.ux[indxX,i_time.+(-1:0)], M.uz[indxZ,i_time.+(-1:0)], M.phi[indxP,i_time.+(-1:0)], M.vx[indxX,i_time.+(-1:0)], M.vz[indxZ,i_time.+(-1:0)], M.Omg[indxP,i_time.+(-1:0)],tInt)




		count::Int64 = 0
		
		while (norm(Re[indx_solve]) > 10.0^dv_norm_tol_exp || norm(Dv) > 10.0^dv_norm_tol_exp )&& count < nwt_iter_max_count

			count += 1
			Ja *=0.0
			Re *=0.0

			for i_el in eachindex(E)
				for i_ke in eachindex(E[i_el].P)
					
						# Obtežba ov vemsnem času
						px = (isnothing(ElementDataIn[i_el].px(0.0)) ? [0.0;0.0] : ElementDataIn[i_el].px((time_st[i_time]+time_st[i_time-1])/2)[i_ke,:])
						pz = (isnothing(ElementDataIn[i_el].pz(0.0)) ? [0.0;0.0] : ElementDataIn[i_el].pz((time_st[i_time]+time_st[i_time-1])/2)[i_ke,:])
						my = (isnothing(ElementDataIn[i_el].my(0.0)) ? [0.0;0.0] : ElementDataIn[i_el].my((time_st[i_time]+time_st[i_time-1])/2)[i_ke,:])


						Px = (isnothing(ElementDataIn[i_el].Px(0.0)) ? [0.0;0.0] : ElementDataIn[i_el].Px((time_st[i_time]+time_st[i_time-1])/2)[[2*i_ke-1,2*i_ke]])
						Pz = (isnothing(ElementDataIn[i_el].Pz(0.0)) ? [0.0;0.0] : ElementDataIn[i_el].Pz((time_st[i_time]+time_st[i_time-1])/2)[[2*i_ke-1,2*i_ke]])
						My = (isnothing(ElementDataIn[i_el].My(0.0)) ? [0.0;0.0] : ElementDataIn[i_el].My((time_st[i_time]+time_st[i_time-1])/2)[[2*i_ke-1,2*i_ke]])


						dlF,F,M.gamma1[E[i_el].indx_int[i_ke],i_time],M.gamma2[E[i_el].indx_int[i_ke],i_time],M.gamma3[E[i_el].indx_int[i_ke],i_time] = Tan_Res(E[i_el].xInt[i_ke], E[i_el].wInt[i_ke], M.ux[E[i_el].indx[i_ke],[i_time-1,i_time]], M.uz[E[i_el].indx[i_ke],[i_time-1,i_time]], M.phi[E[i_el].indxP[i_ke],[i_time-1,i_time]], M.vx[E[i_el].indx[i_ke],[i_time-1,i_time]], M.vz[E[i_el].indx[i_ke],[i_time-1,i_time]], M.Omg[E[i_el].indxP[i_ke],[i_time-1,i_time]], M.gamma1[E[i_el].indx_int[i_ke],i_time-1], M.gamma2[E[i_el].indx_int[i_ke],i_time-1], M.gamma3[E[i_el].indx_int[i_ke],i_time-1], Pvalues[i_el][i_ke], dPvalues[i_el][i_ke], E[i_el].P[i_ke], E[i_el].p0[i_ke], E[i_el].k0[i_ke], ElementDataIn[i_el].C, ElementDataIn[i_el].M, px, pz, my, Px, Pz, My, tInt, E[i_el].pb[i_ke], E[i_el].kb[i_ke], E[i_el].L[i_ke], g, ElementDataIn[i_el].beta)
						
						
						Ja[indx_dof[i_el][i_ke],indx_dof[i_el][i_ke]] += hvcat(length(F),dlF...)
						Re[indx_dof[i_el][i_ke]] += vcat(F...)


				end # i_ke
			end # i_el

			Ja = Ja+mass_m

			#display(sparse(Ja[indx_solve,indx_solve]))

			Dv = -(Ja[indx_solve,indx_solve]\Re[indx_solve])	
			Dv = reshape(Dv,(Int(length(Dv)/(nt-1)),nt-1))
			#display(Dv)
			#display(Re[indx_solve])
			#display(Dv)
			#display(norm(Dv))
			#display(norm(Re[indx_solve]))
			#println()

			#println("\tDv->",norm(Dv))
			#println("\tDv->",norm(Dv),"    \tRe->",norm(Re[indx_solve]))
			#print(".")
			#Dv = round.(Dv/dt,digits=12)


			# Popavek hitrost

			M.vx[indxX,i_time-1:i_time], M.vz[indxZ,i_time-1:i_time], M.Omg[indxP,i_time-1:i_time], M.ux[indxX,i_time-1:i_time], M.uz[indxZ,i_time-1:i_time], M.phi[indxP,i_time-1:i_time] = iteracija(M.ux[indxX,i_time-1:i_time], M.uz[indxZ,i_time-1:i_time], M.phi[indxP,i_time-1:i_time], M.vx[indxX,i_time-1:i_time], M.vz[indxZ,i_time-1:i_time], M.Omg[indxP,i_time-1:i_time], Dv[indxX_solve],Dv[indxZ_solve],Dv[indxP_solve],tInt)



		end # while norm(Dv) > x
		print("it = ",i_time,"    \tt = ",time_st[i_time])
		print("    \titr.cnt.=",count,"\n")
	end # i_time
end
@info "Račun\n\t\t[  Ok  ]"

