#iniciacija funkcij
#   Moduli

include("NonLinBeamMOD.jl");
using LinearAlgebra, .NonLinBeam, SparseArray

# TILE NISO NUJNI ||
	#using BenchmarkTools, Plots, GraphRecipes,Graphs



# SPRAVI TA MACRO V MODUL !
macro assignto(DatVar::Any,ei::Any,Prop::Expr,StrucField::Any)
    ei = eval(ei)
    for i in ei
        eval(
            Meta.parse(string(eval(DatVar))*"["*string(i)*"]."*string(eval(StrucField))*"="*string(eval(Prop)))
        )
    end
end;    




# Precompile
datainit([1 2],Float64.([0 0;0 -1]));
InterpolKoeff([0,1]);


