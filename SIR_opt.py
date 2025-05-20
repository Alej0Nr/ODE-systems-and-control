import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fmin
import os


nu, beta = 1, 0.003
intervalo_tiempo  = [0,20]
intervalo_control = [3,5]

## defino t_eval en funcion de los que estan arriba
t_eval1=np.round(np.arange(intervalo_tiempo[0],intervalo_control[0]+0.1,0.1), 1)
t_eval2=np.round(np.arange(intervalo_control[0],intervalo_control[1]+0.1,0.1), 1)
t_eval3=np.round(np.arange(intervalo_control[1],intervalo_tiempo[1]+0.1,0.1),1)



def SIR(t, z, u=0):
    s, i, r = z #valores de S, I y R
    u = u[0] if isinstance(u, (list, np.ndarray)) else u  # Asegurarse de que u es escalar
    #definicion del sistema
    dsdt = -(1-u)*beta*s*i 
    didt = (1-u)*beta*s*i - nu*i
    drdt = nu*i

    return np.array([dsdt, didt, drdt])



def f(u):
    u = np.array([u]) if not isinstance(u, (list, np.ndarray)) else u
    solucion = solve_ivp(SIR, [0,3], [999,1,0], t_eval=t_eval1, method='RK45')
    s, i, r = solucion.y[0][-1],solucion.y[1][-1],solucion.y[2][-1]
    solucion = solve_ivp(SIR, [3,5], [s,i,r], args=(u), t_eval=t_eval2, method='RK45')
    s = solucion.y[0][-1]

    return abs(s-(nu/beta))

u_opt = fmin(f, x0=0.5, xtol=1e-4, disp=0)
u_opt = u_opt[0] if isinstance(u_opt, (list, np.ndarray)) else u_opt
print(f"Valor óptimo de u: {u_opt}")


def solucionario(u=0):
    solucion = solve_ivp(SIR, [intervalo_tiempo[0],intervalo_control[0]], [999,1,0], t_eval=t_eval1, method='RK45')
    s1, i1, r1 = solucion.y[0],solucion.y[1],solucion.y[2]
    l=len(s1)
    solucion = solve_ivp(SIR, [intervalo_control[0],intervalo_control[1]], [s1[-1],i1[-1],r1[-1]], args=[u], t_eval=t_eval2, method='RK45')
    s2, i2, r2 = solucion.y[0],solucion.y[1],solucion.y[2]
    solucion = solve_ivp(SIR, [intervalo_control[1],intervalo_tiempo[1]], [s2[-1],i2[-1],r2[-1]], t_eval=t_eval3, method='RK45')
    s3, i3, r3 = solucion.y[0],solucion.y[1],solucion.y[2]
    s = np.append(np.append(s1, s2), s3)
    i = np.append(np.append(i1, i2), i3)
    r = np.append(np.append(r1, r2), r3)
    return s,i,r,l



def grafica_SI(us=[0]):
    for j in range(len(us)):
        s,i,r,l = solucionario(us[j])
        if j==len(us)-1:
            plt.plot(s, i, label=f'{us[j]}')
        else:
            plt.plot(s[l:],i[l:], label=f'{us[j]}')
    plt.xlabel('Susceptible')
    plt.ylabel('Infected')
    plt.ylim(0,1000)
    plt.xlim(0,1000)
    plt.grid()
    plt.axvline(x = nu/beta,color= 'black', linestyle = ':', label = r'$R_0$')
    plt.legend()
    os.makedirs("graficos", exist_ok=True)
    plt.savefig(f'graficos/SIR_opt.svg')
    plt.clf()

grafica_SI([u_opt,.5,.3,.6,0])

# grafica_SI(np.arange(0,1.001,.001)) #testing




