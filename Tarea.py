import sympy as sp

"""
las funciones por como estan escritas solo valen para cuando el sistema tiene un unico vector en nullspace

"""
def solve_M(A,B):
    M = sp.Matrix(sp.BlockMatrix([sp.eye(A.shape[0]) - A,-B]))
    return sp.Matrix(M.nullspace())

def eqpair(A,B):
    x = solve_M(A,B)[:-B.shape[1]]
    u = solve_M(A,B)[A.shape[0]]
    return x,u


if __name__ == '__main__':
    A = sp.Matrix([[1,-1],[1/2,-1]])
    B = sp.Matrix([-1,1])
    # sp.pprint(solve_M(A,B))
    print('El par de equilibro es', eqpair(A,B))

