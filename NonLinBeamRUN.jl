#   B R A N J E   P O D A T K O V   I Z   D A T O T E K E

println("Pot do datoteke z podatki")
include("in/"*readline()*".jl")
println()
@info "Vnos podatkov\n\t\t[  Ok  ]"









#   P R E D P O C E S I R A N J E
E = Array{BeamDataProcess}(undef,n_elem)
n_nodes = n_voz
n_int_nodes = 0
for i = 1:n_elem
	global n_nodes 
	global n_int_nodes
	E[i], n_nodes,n_int_nodes = dataprocess(ElementDataIn[i],VozDataIn[ElementDataIn[i].v],n_nodes,n_int_nodes;mtd = Integracija)
end
@info "Procesiranje podatkov\n\t\t[  Ok  ]"









#   P L O T
#

#konstr_img = plotbeams(E,ElementDataIn,VozDataIn)
#display(konstr_img)

#@info "Risnaje konstrukcije\n\t\t[  Ok  ]"

#println("\n Kako nadaljujem?")
#println("0 -> prekini postopek")
#println("Enter\t\t-> nadaljuj račun")
#if (readline() == "0")
	#error("Preklic")
#end








time = collect(ti:dt:tf)
n_time = length(time)

M =  BeamMotion(zeros(n_nodes,n_time),zeros(n_nodes,n_time),zeros(n_nodes,n_time),zeros(n_nodes,n_time),zeros(n_nodes,n_time),zeros(n_nodes,n_time),zeros(n_int_nodes,n_time),zeros(n_int_nodes,n_time),zeros(n_int_nodes,n_time))





begin
#=
	if metoda_t_integracije == "timeelement"
		Ibtime = trig_re_gramschmid([tnodes])
		xt,wt = QuadInt(it)
		xt = (xt.+1)*dt/2*(nt-1)
		wt = wt*dt/2*(nt-1)
		Tvals = TrigValue
		tInt = TimeElement(dt,tnodes,Ibtime,)
	elseif metoda_t_integracije == "midpoint"
		tInt = MidPoint(dt)
	end

=#
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



	if metoda_t_integracije == "timeelement"
		local Ja = zeros(Float64,3*n_nodes*nt,3*n_nodes*nt)
		local Re0 = ones(Float64,3*n_nodes*nt,1)
	elseif metoda_t_integracije == "midpoint"
		local Ja = zeros(Float64,3*n_nodes,3*n_nodes)
		local Re0 = ones(Float64,3*n_nodes)
	end

	local Dv::Array{Float64}


	#vrednosti polinomov v integracijskih točkah
	local Pvalues = map(i1-> map(i2 -> PolyValue(E[i1].xInt[i2],E[i1].P[i2]), eachindex(E[i1].L)),eachindex(E))
	local dPvalues = map(i1-> map(i2 -> PolyValue(E[i1].xInt[i2],E[i1].P[i2];n=1), eachindex(E[i1].L)),eachindex(E))

	indx_dof =  map(i_el->
				 map(i_ke->
					reshape(hcat(E[i_el].indx[i_ke]*3 .-2,E[i_el].indx[i_ke]*3 .-1, E[i_el].indx[i_ke]*3)',
					(3*length(E[i_el].indx[i_ke]))),
				 eachindex(E[i_el].P)),
				eachindex(E))


	mass_m = copy(Ja)
	for i1 = eachindex(E)
		for i2 = eachindex(E[i1].P)
			mi = Mass(E[i1].xInt[i2],E[i1].wInt[i2],Pvalues[i1][i2],dPvalues[i1][i2],ElementDataIn[i1].M)
			mass_m[indx_dof[i1][i2],indx_dof[i1][i2]] = mass_m[indx_dof[i1][i2],indx_dof[i1][i2]]+hvcat(size(mi)[1],mi...)
		end
	end

	@info "Priprava na račun\n\t\t[  Ok  ] "

	@time for i_time in 2:nt-1:n_time
		Dv = [1.;1.]
		Re = Re0

		#Re = ones(Float64,3*n_nodes*nt,1)

		print("\nit = ",i_time,"    \tt = ",time[i_time])


		# Prediktor
		begin
			M.vx[indxX,i_time:i_time+nt-2] = repeat(M.vx[indxX,i_time-1],inner=(1,nt-1))
			M.vz[indxZ,i_time:i_time+nt-2] =  repeat(M.vz[indxZ,i_time-1],inner=(1,nt-1))
			M.Omg[indxP,i_time:i_time+nt-2]=  repeat(M.Omg[indxP,i_time-1],inner=(1,nt-1))

			M.ux[indxX,i_time:i_time+nt-2] = repeat(M.ux[indxX,i_time-1] + M.vx[indxX,i_time:i_time+nt-2]*dt,inner=(1,nt-1))
			M.uz[indxZ,i_time:i_time+nt-2] =  repeat(M.uz[indxZ,i_time-1] + M.vz[indxZ,i_time:i_time+nt-2]*dt,inner=(1,nt-1))
			M.phi[indxP,i_time:i_time+nt-2]=  repeat(M.phi[indxP,i_time-1] + M.Omg[indxP,i_time:i_time+nt-2]*dt,inner=(1,nt-1))
		end



		count::Int64 = 0
		
		while (norm(Re[indx_solve]) > 10.0^dv_norm_tol_exp || norm(Dv) > 10.0^dv_norm_tol_exp )&& count < nwt_iter_max_count

			count += 1
			Ja *=0.0
			Re *=0.0

			for i_el in eachindex(E)
				for i_ke in eachindex(E[i_el].P)
					
						# Obtežba ov vemsnem času
						px = (isnothing(ElementDataIn[i_el].px(0.0)) ? [0.0;0.0] : ElementDataIn[i_el].px((time[i_time]+time[i_time-1])/2)[i_ke,:]) 	
						pz = (isnothing(ElementDataIn[i_el].pz(0.0)) ? [0.0;0.0] : ElementDataIn[i_el].pz((time[i_time]+time[i_time-1])/2)[i_ke,:]) 
						my = (isnothing(ElementDataIn[i_el].my(0.0)) ? [0.0;0.0] : ElementDataIn[i_el].my((time[i_time]+time[i_time-1])/2)[i_ke,:]) 


						Px = (isnothing(ElementDataIn[i_el].Px(0.0)) ? [0.0;0.0] : ElementDataIn[i_el].Px((time[i_time]+time[i_time-1])/2)[[2*i_ke-1,2*i_ke]])
						Pz = (isnothing(ElementDataIn[i_el].Pz(0.0)) ? [0.0;0.0] : ElementDataIn[i_el].Pz((time[i_time]+time[i_time-1])/2)[[2*i_ke-1,2*i_ke]])
						My = (isnothing(ElementDataIn[i_el].My(0.0)) ? [0.0;0.0] : ElementDataIn[i_el].My((time[i_time]+time[i_time-1])/2)[[2*i_ke-1,2*i_ke]])


						dlF,F,M.gamma1[E[i_el].indx_int[i_ke],i_time],M.gamma2[E[i_el].indx_int[i_ke],i_time],M.gamma3[E[i_el].indx_int[i_ke],i_time] = Tan_Res(E[i_el].xInt[i_ke],E[i_el].wInt[i_ke],M.ux[E[i_el].indx[i_ke],[i_time-1,i_time]],M.uz[E[i_el].indx[i_ke],[i_time-1,i_time]],M.phi[E[i_el].indx[i_ke],[i_time-1,i_time]],M.vx[E[i_el].indx[i_ke],[i_time-1,i_time]],M.vz[E[i_el].indx[i_ke],[i_time-1,i_time]],M.Omg[E[i_el].indx[i_ke],[i_time-1,i_time]], M.gamma1[E[i_el].indx_int[i_ke],i_time-1], M.gamma2[E[i_el].indx_int[i_ke],i_time-1], M.gamma3[E[i_el].indx_int[i_ke],i_time-1], Pvalues[i_el][i_ke], dPvalues[i_el][i_ke], E[i_el].P[i_ke], E[i_el].p0[i_ke], E[i_el].k0[i_ke], ElementDataIn[i_el].C, ElementDataIn[i_el].M, px, pz, my, Px, Pz, My, dt, E[i_el].pb[i_ke], E[i_el].kb[i_ke], E[i_el].L[i_ke], g)
						
						
						Ja[indx_dof[i_el][i_ke],indx_dof[i_el][i_ke]] += hvcat(length(F),dlF...)
						Re[indx_dof[i_el][i_ke]] += vcat(F...)


				end # i_ke
			end # i_el

			Ja = Ja+mass_m

			Dv = -(Ja[indx_solve,indx_solve]\Re[indx_solve])	
			Dv = reshape(Dv,(Int(length(Dv)/(nt-1)),nt-1))
			#display(Ja[indx_solve,indx_solve])
			#display(Re[indx_solve])
			#display(Dv)
			#display(norm(Dv))
			#println()

			#println("\tDv->",norm(Dv))
			#println("\tDv->",norm(Dv),"\tRe->",norm(Re[indx_solve]))
			#print(".")
			#Dv = round.(Dv/dt,digits=12)


			# Popavek hitrost

			if metoda_t_integracije == "timeelement"

				M.vx[indxX,i_time:i_time+nt-2] += Dv[1:3:end,:]
				M.vz[indxZ,i_time:i_time+nt-2] += Dv[2:3:end,:]
				M.Omg[indxP,i_time:i_time+nt-2] += Dv[3:3:end,:]


				for i in 1:nt-1
					M.ux[indxX,i_time+i-1] += Integrate(t->PolyValue(t,Ibtime*hcat(M.vx[indxX,i_time-1],Dv[1:3:end,:] )'),0.0,i*dt)
					M.uz[indxZ,i_time+i-1] += Integrate(t->PolyValue(t,Ibtime*hcat(M.vz[indxZ,i_time-1],Dv[2:3:end,:] )'),0.0,i*dt)
					M.phi[indxP,i_time+i-1] += Integrate(t->PolyValue(t,Ibtime*hcat(M.Omg[indxP,i_time-1],Dv[3:3:end,:] )'),0.0,i*dt)

				end

			elseif metoda_t_integracije == "midpoint"

				M.vx[indxX,i_time] += Dv[indxX_solve]*2
				M.vz[indxZ,i_time] += Dv[indxZ_solve]*2
				M.Omg[indxP,i_time]+= Dv[indxP_solve]*2


				# Popravek pomikov
				M.ux[indxX,i_time] =  M.ux[indxX,i_time-1] + #=Dv[indxX]*dt=# sum(M.vx[indxX,i_time-1:i_time],dims=2)*dt/2
				M.uz[indxZ,i_time] =  M.uz[indxZ,i_time-1] + #=Dv[indxZ]*dt=# sum(M.vz[indxZ,i_time-1:i_time],dims=2)*dt/2
				M.phi[indxP,i_time]=  M.phi[indxP,i_time-1]+ #=Dv[indxP]*dt=# sum(M.Omg[indxP,i_time-1:i_time],dims=2)*dt/2
			end


		end # while norm(Dv) > x
		print("    \titr.cnt.=",count)
	end # i_time
end
println("[  Ok  ]  Račun")

