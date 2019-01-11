using SparseArrays
using MatrixMarket
using DelimitedFiles

# verify correctness of triangular solve ...
# by runing it in another programming language

matrs = ["s_small", "s_medium", "s_torso", "s_tsopf"]


for matr in ["s_torso"]
    L = MatrixMarket.mmread("data/$(matr)L.mtx")
    @show cond(Matrix(L))

    # b = MatrixMarket.mmread("data/$(matr)b.mtx")
    # x = L \ b
    
    # y = sparse(L*x);

    # @show y
    # @show b

    # open("data/sol_$(matr)", "w") do io
    #     DelimitedFiles.writedlm(io, x)
    # end

end