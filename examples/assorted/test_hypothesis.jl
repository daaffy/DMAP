using LinearAlgebra

n = 2000
B = rand(10*n,10*n)
A = rand(n,n)
@time eigen(A)
