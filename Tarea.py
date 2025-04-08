import sympy as sp


def solve_M(A,B): #conviene usar Np.Arrayj
    M = sp.Matrix(sp.BlockMatrix([sp.eye(A.rows) - A,-B])) #eye_like(A)
    #usar concatenate 
    return M.nullspace()

def eqPair(A,B):
    espacioNulo = solve_M(A,B)
    listaPares=[]
    for v in espacioNulo:
        listaPares.append((v[:-B.cols],v[A.rows:]))
    return listaPares

if __name__ == '__main__':
    A = sp.Matrix([[1,-1],[1/2,-1]])
    # sp.pprint(A)
    B = sp.Matrix([-1,1])
    # sp.pprint(B)
    for par in eqPair(A,B):
        sp.pprint(par)
