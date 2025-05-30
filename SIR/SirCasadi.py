from casadi import *
from libreriaSIR import *


desde = 3
cant_dias = 2
N=100
dt = cant_dias/N

sim = SIRsim([999,1,0],(0,desde), params=(1,0.003))
sim.simulacion()

def EE(sols=False):
    nu = sim.params[0]
    beta = sim.params[1]
    sir = Opti()
    # Definimos nuestras variables
    S = sir.variable(N+1)
    I = sir.variable(N+1)
    R = sir.variable(N+1)
    u = sir.variable(N)

    # Valores iniciales
    sir.subject_to(S[0] == sim.S[-1])
    sir.subject_to(I[0] == sim.I[-1])
    sir.subject_to(R[0] == sim.R[-1])

    # Ecuaciones del modelo (Euler explícito)
    for k in range(N):
        sir.subject_to(u[k]<=1)
        sir.subject_to(u[k]>=0)
        sir.subject_to(S[k+1] == S[k] - dt * (1 - u[k]) * beta * S[k] * I[k])
        sir.subject_to(I[k+1] == I[k] + dt * ((1 - u[k]) * beta * S[k] * I[k] - nu * I[k]))
        sir.subject_to(R[k+1] == R[k] + dt * nu * I[k])

    for k in range(1,N-1):
        sir.subject_to(u[k]==u[k+1])

    # Restricciones sobre las variables

    sir.subject_to(S >= 0)
    sir.subject_to(I >= 0)
    sir.subject_to(R >= 0)

    # Función objetivo
    S_star = sim.R_star  # valor deseado de S al final
    sir.minimize((S[-1] - S_star)**2 + I[-1]**2)

    # Resolver
    p_opts = {}
    s_opts = {'print_level': 0, 'sb': 'yes'}
    sir.solver('ipopt',p_opts,s_opts)
    sol = sir.solve()

    

    sirS = SIRsim([999,1,0],(0,16),(desde,desde+cant_dias),control=sol.value(u[-1]))
    sirS.simulacion()
    if sols:
        print(f'{'control':<7}={sirS.control}')
        print(f"{'S_c':<7}={sirS.dc[0]}")
        print(f'{'I_c':<7}={sirS.dc[1]}')
        print(f'{'R_c':<7}={sirS.dc[2]}')
        print(f'{'S_star':<7}={S_star}')
    
    sirS.grafica_temporal()
    sirS.show()

def RK2(sols=False):
    nu = sim.params[0]
    beta = sim.params[1]
    sirRK2 = Opti()

    S = sirRK2.variable(N+1)
    I = sirRK2.variable(N+1)
    R = sirRK2.variable(N+1)
    u = sirRK2.variable(N)

    # Valores iniciales
    sirRK2.subject_to(S[0] == sim.S[-1])
    sirRK2.subject_to(I[0] == sim.I[-1])
    sirRK2.subject_to(R[0] == sim.R[-1])

    # Ecuaciones del modelo (RK2)
    for k in range(N):
        sirRK2.subject_to(u[k]<=1)
        sirRK2.subject_to(u[k]>=0)
        # k1
        S_k = S[k]
        I_k = I[k]
        R_k = R[k]
        u_k = u[k]

        dS1 = -(1 - u_k) * beta * S_k * I_k
        dI1 = (1 - u_k) * beta * S_k * I_k - nu * I_k
        dR1 = nu * I_k

        # estado intermedio
        S_mid = S_k + dt/2 * dS1
        I_mid = I_k + dt/2 * dI1
        R_mid = R_k + dt/2 * dR1

        # k2 (en el punto medio)
        dS2 = -(1 - u_k) * beta * S_mid * I_mid
        dI2 = (1 - u_k) * beta * S_mid * I_mid - nu * I_mid
        dR2 = nu * I_mid

        # actualización final
        sirRK2.subject_to(S[k+1] == S_k + dt * dS2)
        sirRK2.subject_to(I[k+1] == I_k + dt * dI2)
        sirRK2.subject_to(R[k+1] == R_k + dt * dR2)


    for k in range(1,N-1):
        sirRK2.subject_to(u[k]==u[k+1])

    # Restricciones sobre las variables
    sirRK2.subject_to(S >= 0)
    sirRK2.subject_to(I >= 0)
    sirRK2.subject_to(R >= 0)

    # Función objetivo
    S_star = sim.R_star  # valor deseado de S al final
    sirRK2.minimize((S[-1] - S_star)**2 + I[-1]**2)

    # Resolver
    p_opts = {}
    s_opts = {'print_level': 0, 'sb': 'yes'}
    sirRK2.solver('ipopt',p_opts,s_opts)
    solRK2 = sirRK2.solve()

    sirS = SIRsim([999,1,0],(0,16),(desde,desde+cant_dias),control=solRK2.value(u[-1]))
    sirS.simulacion()
    if sols:
        print(f'{'control':<7}={sirS.control}')
        print(f"{'S_c':<7}={sirS.dc[0]}")
        print(f'{'I_c':<7}={sirS.dc[1]}')
        print(f'{'R_c':<7}={sirS.dc[2]}')
        print(f'{'S_star':<7}={S_star}')
    sirS.grafica_temporal()
    sirS.show()


EE(True)
RK2(True)
