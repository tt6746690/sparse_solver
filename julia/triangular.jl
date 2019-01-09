using MatrixMarket
using DelimitedFiles

for matr in ["s_small", "s_torso", "s_tsopf"]
    L = MatrixMarket.mmread("data/$(matr)L.mtx")
    b = MatrixMarket.mmread("data/$(matr)b.mtx")
    x = L \ b
    open("data/sol_$(matr)", "w") do io
        DelimitedFiles.writedlm(io, x)
    end
end