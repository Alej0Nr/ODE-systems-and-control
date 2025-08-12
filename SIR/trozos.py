from casadi import *
from parametros import *
from libreriaSIR import *
import numpy as np
import matplotlib.pyplot as plt
import os
os.makedirs("graficos", exist_ok=True)


sim = SIRsim([999,1,0],(0,desde), params=(1,0.003))
sim.simulacion()


def RK2():
    nu = sim.params[0]
    beta = sim.params[1]
    sirRK2 = Opti()

    S = sirRK2.variable(N+1)
    I = sirRK2.variable(N+1)
    R = sirRK2.variable(N+1)
    u = sirRK2.variable(N+1)

    # Valores iniciales
    sirRK2.subject_to(S[0] == sim.S[-1])
    sirRK2.subject_to(I[0] == sim.I[-1])
    sirRK2.subject_to(R[0] == sim.R[-1])
    #sirRK2.subject_to(u[0] == 0)

    # Ecuaciones del modelo (RK2)
    for k in range(N):
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


    # Restricciones sobre las variables
    sirRK2.subject_to(sirRK2.bounded(0,u,1))
    sirRK2.subject_to(S >= 0)
    sirRK2.subject_to(I >= 0)
    sirRK2.subject_to(R >= 0)
    flag=desde
    for k in range(N):
        if int(flag)==int(flag+dt):
            #print(k, flag, flag+dt, 'TRUE')
            sirRK2.subject_to(u[k]==u[k+1])
        flag+=dt
    # Función objetivo
    S_star = sim.R_star  # valor deseado de S al final
    sirRK2.minimize((S[-1] - S_star)**2) if fcosto==1 else sirRK2.minimize((S[-1] - S_star)**2 + sumsqr(u))

    # Resolver
    p_opts = {}
    s_opts = {'print_level': 1, 'sb': 'yes'}
    sirRK2.solver('ipopt',p_opts,s_opts)
    solRK2 = sirRK2.solve()
    # print(solRK2.value(I[-1]))

    tiempos = np.linspace(desde, desde+cant_dias, N+1)
    plt.subplot(211)
    plt.plot([0,desde],[0,0], color='navy')
    plt.plot(tiempos,solRK2.value(u),label=r'$u(t)$',color='navy', drawstyle='steps')
    plt.plot([desde+cant_dias,desde+cant_dias+2],[0,0],color='navy')
    plt.ylim(-0.05,1.05)
    plt.yticks(np.arange(0,1.05,.1))
    plt.xticks(np.arange(0,desde+cant_dias+2.1,1))
    plt.grid()
    plt.legend()

    plt.subplot(212)
    sim.grafica_temporal()
    plt.plot(tiempos,solRK2.value(S),color='mediumblue')
    plt.plot(tiempos,solRK2.value(I),color='red')
    plt.plot(tiempos,solRK2.value(R),color='green')
    sim.cambio_Ci([solRK2.value(S)[-1],solRK2.value(I)[-1],solRK2.value(R)[-1]])
    sim.cambio_tiempo_sim((desde+cant_dias,desde+cant_dias+2))
    sim.simulacion()
    sim.grafica_temporal(label=False)
    plt.xticks(np.arange(0,desde+cant_dias+2.1,1))
    plt.axhline(y = S_star,color= 'black', linestyle = ':', label = r'$S^\ast$')
    # plt.plot(desde,0,marker= 6, color='dodgerblue')
    # plt.plot(desde+cant_dias,0,marker=6,color='dodgerblue')
    plt.grid()
    plt.legend()
    plt.savefig(f'graficos/trozos{fcosto}.svg') if guardar else None
    plt.show() if graficas else None
    plt.clf()

    if planoSI:
        sim2 = SIRsim([999,1,0],(0,desde), params=(1,0.003))
        sim2.simulacion()
        sim2.grafica_SI()
        plt.axvspan(sim2.S[-1],solRK2.value(S)[-1], color='lavender')
        plt.plot(solRK2.value(S),solRK2.value(I),color="purple")
        sim2.cambio_Ci([solRK2.value(S)[-1],solRK2.value(I)[-1],solRK2.value(R)[-1]])
        sim2.cambio_tiempo_sim((desde+cant_dias,desde+cant_dias+10))
        sim2.simulacion()
        sim2.grafica_SI()
        plt.xlabel('Suceptibles')
        plt.ylabel('Infectados')
        plt.savefig(f'graficos/SItrozos{fcosto}.svg') if guardar else None
        plt.show() if graficas else None
        plt.clf()


RK2()