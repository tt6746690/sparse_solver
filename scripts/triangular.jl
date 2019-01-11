using SparseArrays
using MatrixMarket
using DelimitedFiles
using LinearAlgebra

matrs = ["s_small", "s_medium", "s_torso", "s_tsopf"]


for matr in matrs
    L = MatrixMarket.mmread("data/$(matr)L.mtx")
    @show cond(Matrix(L))

    b = MatrixMarket.mmread("data/$(matr)b.mtx")
    x = L \ b
    
    open("data/sol_$(matr)", "w") do io
        DelimitedFiles.writedlm(io, x)
    end
end