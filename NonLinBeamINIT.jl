#iniciacija funkcij

#   Moduli
begin
    include("NonLinBeamMOD.jl");
    using LinearAlgebra, .NonLinBeam, SparseArrays, Graphs
    using BenchmarkTools, Plots, GraphRecipes
    macro assignto(DatVar::Any,ei::Any,Prop::Expr,StrucField::Any)
        ei = eval(ei)
        for i in ei
            eval(
                Meta.parse(string(eval(DatVar))*"["*string(i)*"]."*string(eval(StrucField))*"="*string(eval(Prop)))
            )
        end
    end    
    #pgfplotsx()
end;
