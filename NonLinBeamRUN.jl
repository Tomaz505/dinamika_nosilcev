#   B R A N J E   P O D A T K O V   I Z   D A T O T E K E

println("Pot do datoteke z podatki")
include("in/"*readline()*".jl")
println()
@info "Vnos podatkov\n\t\t[  Ok  ]"









#   P R E D P O C E S I R A N J E
E = Array{BeamDataProcess}(undef,n_elem)
n_nodes = n_voz
for i = 1:n_elem
	global n_nodes 
	E[i], n_nodes = dataprocess(ElementDataIn[i],VozDataIn[ElementDataIn[i].v],n_nodes;mtd = Integracija)
end
@info "Procesiranje podatkov\n\t\t[  Ok  ]"









#   P L O T
#
#=
konstr_img = plotbeams(E,ElementDataIn,VozDataIn)
display(konstr_img)

@info "Risnaje konstrukcije\n\t\t[  Ok  ]"


println("\n Kako nadaljujem?")
println("0\t\t-> preklici postopek")
println("Enter\t\t-> nadaljuj račun")
if (readline() == "0")
	error("Preklic")
end
=#






Ibtime = re_gramschmid([range(0.0,(nt-1)*dt,length=nt)|>collect])

xt,wt = QuadInt(it)
xt = (xt.+1)*dt/2*(nt-1)
wt = wt*dt/2*(nt-1)

time = collect(ti:dt:tf)
n_time = length(time)
M =  BeamMotion(zeros(n_nodes,n_time),zeros(n_nodes,n_time),zeros(n_nodes,n_time),zeros(n_nodes,n_time),zeros(n_nodes,n_time),zeros(n_nodes,n_time))


begin


	#Komponente za razvoj interpolacije (vključuje odvod če je interpoliran)
	#n_nodes = n_voz + sum(map(i -> .+(size.(E[i].P)...)[1]-2-(length(ElementDataIn[i].div2)-1)*(1+Int(ElementDataIn[i].Ci)),1:n_elem))




	indx = setdiff(vcat(map(j->(VozDataIn[j].i .+[0;n_nodes;2*n_nodes]).*(-VozDataIn[j].Supp'.+1), 1:n_voz )...),[0])
	indx = setdiff(1:3*n_nodes,indx)

	indxX = indx[findall(indx.<= n_nodes)]
	indxZ = indx[findall(n_nodes+1 .<= indx .<= 2*n_nodes)].-n_nodes
	indxP = indx[findall(2*n_nodes+1 .<=indx.<= 3*n_nodes)].-2*n_nodes

	indx_solve = vcat(map(i-> sort(vcat(3*indxX.-2,3*indxZ.-1,3*indxP)) .+n_nodes*(i-1)*3, 2:nt)...)
	println("[  Ok  ]  Priprava na račun")
	#=
	i_time = 2
	i_el=1
	i_ke=1

	dlF,F = Tan_Res(E[i_el].xInt[i_ke],E[i_el].wInt[i_ke],M.ux[E[i_el].indx[i_ke],[i_time-1,i_time]],M.uz[E[i_el].indx[i_ke],[i_time-1,i_time]],M.phi[E[i_el].indx[i_ke],[i_time-1,i_time]],M.vx[E[i_el].indx[i_ke],[i_time-1,i_time]],M.vz[E[i_el].indx[i_ke],[i_time-1,i_time]],M.Omg[E[i_el].indx[i_ke],[i_time-1,i_time]],E[i_el].P[i_ke],E[i_el].p0[i_ke],E[i_el].k0[i_ke],ElementDataIn[i_el].C,reshape(ElementDataIn[i_el].M[i_ke,:],2),ElementDataIn[i_el].px(t[i_time])[i_ke,:],ElementDataIn[i_el].pz(t[i_time])[i_ke,:],ElementDataIn[i_el].my(t[i_time])[i_ke,:],dt,E[i_el].pb[i_ke],E[i_el].kb[i_ke],E[i_el].L[i_ke],g)
	=#
	for i_time in 2:nt-1:n_time
		Dv = [1.;1.]

		println("i_time = ",i_time,"\t\tt = ",time[i_time])
		println()


		# Prediktor
		M.vx[indxX,i_time:i_time+nt-2] = repeat(M.vx[indxX,i_time-1],inner=(1,nt-1))
		M.vz[indxZ,i_time:i_time+nt-2] =  repeat(M.vz[indxZ,i_time-1],inner=(1,nt-1))
		M.Omg[indxP,i_time:i_time+nt-2]=  repeat(M.Omg[indxP,i_time-1],inner=(1,nt-1))


		M.ux[indxX,i_time:i_time+nt-2] = repeat(M.ux[indxX,i_time-1],inner=(1,nt-1))
		M.uz[indxZ,i_time:i_time+nt-2] =  repeat(M.uz[indxZ,i_time-1],inner=(1,nt-1))
		M.phi[indxP,i_time:i_time+nt-2]=  repeat(M.phi[indxP,i_time-1],inner=(1,nt-1))


		count::Int64 = 0
		
		while norm(Dv) > 10.0^dv_norm_tol_exp && count < nwt_iter_max_count
			count += 1
			
			# Tangentna in rezidual
			Ja = spzeros(Float64,3*n_nodes*nt,3*n_nodes*nt)
			Re = zeros(Float64,3*n_nodes*nt,1)

			
			for i_el in eachindex(E)
				for i_ke in eachindex(E[i_el].P)
					indx_dof = vcat(map(i-> n_nodes*(i-1)*3 .+ reshape(hcat(E[i_el].indx[i_ke]*3 .-2,E[i_el].indx[i_ke]*3 .-1, E[i_el].indx[i_ke]*3)',(3*length(E[i_el].indx[i_ke]))),1:nt)...)
					

					JR = Integrate( (t,x)-> TanResAtTX(x,
									    t,
									    M.vx[E[i_el].indx[i_ke],i_time-1:i_time-1+nt-1],
									    M.vz[E[i_el].indx[i_ke],i_time-1:i_time-1+nt-1],
									    M.Omg[E[i_el].indx[i_ke],i_time-1:i_time-1+nt-1],
									    M.ux[E[i_el].indx[i_ke],i_time-1],
									    M.uz[E[i_el].indx[i_ke],i_time-1],
									    M.phi[E[i_el].indx[i_ke],i_time-1],
									    E[i_el].P[i_ke],
									    Ibtime,
									    E[i_el].p0[i_ke][findfirst(E[i_el].xInt[i_ke] .== x)],
									    E[i_el].k0[i_ke][findfirst(E[i_el].xInt[i_ke] .== x)],
									    ElementDataIn[i_el].C,
									    ElementDataIn[i_el].M,
									    (isnothing(ElementDataIn[i_el].px(0.0)) ? 0.0 : [1 x]*[0.5;-0.5;;0.5;0.5]*ElementDataIn[i_el].px(time[i_time]+(t+1)/2*dt)[i_ke,:]), 
									    (isnothing(ElementDataIn[i_el].pz(0.0)) ? 0.0 : [1 x]*[0.5;-0.5;;0.5;0.5]*ElementDataIn[i_el].pz(time[i_time]+(t+1)/2*dt)[i_ke,:]),
									    (isnothing(ElementDataIn[i_el].my(0.0)) ? 0.0 : [1 x]*[0.5;-0.5;;0.5;0.5]*ElementDataIn[i_el].my(time[i_time]+(t+1)/2*dt)[i_ke,:]), 
									    E[i_el].L[i_ke],
									    dt),(xt,wt),(E[i_el].xInt[i_ke],E[i_el].wInt[i_ke]))
					

					pR = Integrate( t-> parResAtT(t,
								       M.vx[E[i_el].indx[i_ke],i_time-1:i_time+nt-2],
								       M.vz[E[i_el].indx[i_ke],i_time-1:i_time+nt-2],
								       M.Omg[E[i_el].indx[i_ke],i_time-1:i_time+nt-2],
								       M.ux[E[i_el].indx[i_ke],i_time-1],
								       M.uz[E[i_el].indx[i_ke],i_time-1],
								       M.phi[E[i_el].indx[i_ke],i_time-1],
								       E[i_el].P[i_ke],
								       Ibtime,
								       E[i_el].pb[i_ke],
								       E[i_el].kb[i_ke],
								       ElementDataIn[i_el].C,
								       (isnothing(ElementDataIn[i_el].Px(0.0)) ? [0.0;0.0] : ElementDataIn[i_el].Px(time[i_time-1]+(t+1)/2*dt)[[i_ke,i_ke+1]]), 
								       (isnothing(ElementDataIn[i_el].Pz(0.0)) ? [0.0;0.0] : ElementDataIn[i_el].Pz(time[i_time-1]+(t+1)/2*dt)[[i_ke,i_ke+1]]),
								       (isnothing(ElementDataIn[i_el].My(0.0)) ? [0.0;0.0] : ElementDataIn[i_el].My(time[i_time-1]+(t+1)/2*dt)[[i_ke,i_ke+1]]), 
								       E[i_el].L[i_ke],
								       dt),(xt,wt))


					Ja[indx_dof,indx_dof] += JR[:,1:end-1]
					Re[indx_dof] += JR[:,end]+ pR
					

					#=
					# Obtežba ov vemsnem času
					px = (isnothing(ElementDataIn[i_el].px(0.0)) ? [0.0;0.0] : ElementDataIn[i_el].px((t[i_time]+t[i_time-1])/2)[i_ke,:]) 	
					pz = (isnothing(ElementDataIn[i_el].pz(0.0)) ? [0.0;0.0] : ElementDataIn[i_el].pz((t[i_time]+t[i_time-1])/2)[i_ke,:]) 
					my = (isnothing(ElementDataIn[i_el].my(0.0)) ? [0.0;0.0] : ElementDataIn[i_el].my((t[i_time]+t[i_time-1])/2)[i_ke,:]) 
					if length(E[i_el].P) == 1
						Px = (ElementDataIn[i_el].Px((t[i_time]+t[i_time-1])/2)[1],ElementDataIn[i_el].Px((t[i_time]+t[i_time-1])/2)[2])
						Pz = (ElementDataIn[i_el].Pz((t[i_time]+t[i_time-1])/2)[1],ElementDataIn[i_el].Pz((t[i_time]+t[i_time-1])/2)[2])
						My = (ElementDataIn[i_el].My((t[i_time]+t[i_time-1])/2)[1],ElementDataIn[i_el].My((t[i_time]+t[i_time-1])/2)[2])
					elseif i_ke == 1
						Px = (ElementDataIn[i_el].Px((t[i_time]+t[i_time-1])/2)[1],0.)
						Pz = (ElementDataIn[i_el].Pz((t[i_time]+t[i_time-1])/2)[1],0.)
						My = (ElementDataIn[i_el].My((t[i_time]+t[i_time-1])/2)[1],0.)
					elseif i_ke == length(E[i_el].P)	
						Px = (0.,ElementDataIn[i_el].Px((t[i_time]+t[i_time-1])/2)[2])
						Pz = (0.,ElementDataIn[i_el].Pz((t[i_time]+t[i_time-1])/2)[2])
						My = (0.,ElementDataIn[i_el].My((t[i_time]+t[i_time-1])/2)[2])
					else
						Px = (0.,0.)
						Pz = (0.,0.)
						My = (0.,0.)
					end


					dlF,F = Tan_Res(E[i_el].xInt[i_ke],E[i_el].wInt[i_ke],M.ux[E[i_el].indx[i_ke],[i_time-1,i_time]],M.uz[E[i_el].indx[i_ke],[i_time-1,i_time]],M.phi[E[i_el].indx[i_ke],[i_time-1,i_time]],M.vx[E[i_el].indx[i_ke],[i_time-1,i_time]],M.vz[E[i_el].indx[i_ke],[i_time-1,i_time]],M.Omg[E[i_el].indx[i_ke],[i_time-1,i_time]],E[i_el].P[i_ke],E[i_el].p0[i_ke],E[i_el].k0[i_ke],ElementDataIn[i_el].C,ElementDataIn[i_el].M, px, pz, my, Px, Pz, My,dt,E[i_el].pb[i_ke],E[i_el].kb[i_ke],E[i_el].L[i_ke],g)
					
					
				
					
					Ja[indx_dof,indx_dof] += hvcat(length(F),dlF...)
					Re[indx_dof] += vcat(F...)
					=#

				end # i_ke
			end # i_el
			#Ja = round.(Ja,digits = 13)
			#Re = round.(Re,digits = 13)
			Dv = -(Ja[indx_solve,indx_solve]\Re[indx_solve])	
			Dv = reshape(Dv,(Int(length(Dv)/(nt-1)),nt-1))

			#Dv = round.(Dv/dt,digits=12) 
			


			println("\t\tDv = ")
			display(Dv)
			println("\t|Dv| = ",norm(Dv),"\n")
			println("\t\tJa = ")
			display(Ja[indx_solve,indx_solve])
			println("\n")
			println("\t\tRe = ")
			display(Re[indx_solve])
			println("\t|Re| = ",norm(Re[indx_solve]),"\n")


			# Popavek hitrosti
			M.vx[indxX,i_time:i_time+nt-2] += Dv[1:3:end,:]
			M.vz[indxZ,i_time:i_time+nt-2] += Dv[2:3:end,:]
			M.Omg[indxP,i_time:i_time+nt-2] += Dv[3:3:end,:]


			for i in 1:nt-1
				M.ux[indxX,i_time+i-1] += Integrate(t->PolyValue(t,Ibtime*hcat(M.vx[indxX,i_time-1],Dv[1:3:end,:] )'),0.0,i*dt)
				M.uz[indxZ,i_time+i-1] += Integrate(t->PolyValue(t,Ibtime*hcat(M.vz[indxZ,i_time-1],Dv[2:3:end,:] )'),0.0,i*dt)
				M.phi[indxP,i_time+i-1] += Integrate(t->PolyValue(t,Ibtime*hcat(M.Omg[indxP,i_time-1],Dv[3:3:end,:] )'),0.0,i*dt)
			
			end	

			#= Popravek pomikov
			M.ux[indxX,i_time:i_time+length(time_interpolation[1])-1] += Dv[1:3:end,:]*dt
			M.uz[indxZ,i_time:i_time+length(time_interpolation[1])-1] += Dv[2:3:end,:]*dt
			M.phi[indxP,i_time:i_time+length(time_interpolation[1])-1]+= Dv[3:3:end,:]*dt
			=#

			#sleep(1)
		end # while norm(Dv) > x
	end # i_time
end
println("[  Ok  ]  Račun")

