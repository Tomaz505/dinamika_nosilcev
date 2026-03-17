#iniciacija funkcij
#   Moduli


include("NonLinBeamMOD.jl");
#=
import REPL
using REPL.TerminalMenus
options = ["Modal","Dynamic","Plot",...]
analize = request(MultiSelectMenu(options))

Uporaba Set
=#

using .NonLinBeam, SparseArrays,LinearAlgebra, Latexify
using Plots
theme(:dao)

#unicodeplots()
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
