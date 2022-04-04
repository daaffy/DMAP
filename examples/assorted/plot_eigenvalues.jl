using Plots

inds = 1:20
display(plot(inds,real(λ)[end:-1:end-inds[end]+1]))