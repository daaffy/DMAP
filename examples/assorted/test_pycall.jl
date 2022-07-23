
using PyCall

# using Pkg
# ENV["PYTHON"] = "Users/jackh/opt/anaconda3"
# Pkg.build("PyCall")

function times_two(x::Number)
    py"""
    
    def times_two(x):
        return 2*x
    
    """

    return py"times_two"(x)
end

function make_sparse()
    py"""
    import numpy as np
    import scipy
    # from scipy.sparse import csr_matrix

    def test():

        x = np.array([[2,3,4]])
        y = x*2
        return y
    """

    py"test"()

end

make_sparse()