import sympy as sp

def solve_M(A: sp.Matrix,B: sp.Matrix):
    M = sp.Matrix(sp.BlockMatrix([sp.eye(A.shape[0]) - A,-B]))  
    return M.nullspace()


def eqpair(A,B):
    espacioNulo = solve_M(A,B)
    listaPares=[]
    for v in espacioNulo:
        listaPares.append((v[:-B.shape[1]],v[B.shape[0]:]))
    return listaPares

if __name__ == '__main__':
    A = sp.Matrix([[1,-1],[1/2,-1]])
    # sp.pprint(A)
    B = sp.Matrix([-1,1])
    # sp.pprint(B)
    for par in eqpair(A,B):
        sp.pprint(par)
