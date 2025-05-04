import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp



def SIR(t, z, nu, beta, u=0):
    s, i, r = z
    dsdt = -(1-u)*beta*s*i
    didt = (1-u)*beta*s*i - nu*i
    drdt = nu*i

    return [dsdt, didt, drdt]

def plano_temporal(parameters, desde = None, hasta = None):
    if desde == None or hasta == None:
        solucion = solve_ivp(SIR, t_span, z0, t_eval= np.linspace(t_span[0], t_span[1], 1000), args= parameters, method='RK45')
        plt.plot(solucion.t, solucion.y[0], label='S(t)', color = 'mediumblue')
        plt.plot(solucion.t, solucion.y[1], label='I(t)', color = 'green')
        plt.plot(solucion.t, solucion.y[2], label='R(t)', color = 'red')
        print([solucion.y[0][-1],solucion.y[1][-1],solucion.y[2][-1]])
    elif len(parameters)==3:
        solucion1 = solve_ivp(SIR, (t_span[0],desde), z0, t_eval= np.linspace(t_span[0], desde, 1000), args= parameters[:-1], method='RK45')
        z0n=[solucion1.y[0][-1],solucion1.y[1][-1],solucion1.y[2][-1]]
        solucion2 = solve_ivp(SIR, (desde, hasta), z0n, t_eval= np.linspace(desde, hasta, 1000), args= parameters, method='RK45')
        z0n=[solucion2.y[0][-1],solucion2.y[1][-1],solucion2.y[2][-1]]
        solucion3 = solve_ivp(SIR, (hasta,t_span[1]), z0n, t_eval= np.linspace(hasta, t_span[1], 1000), args= parameters[:-1], method='RK45')
        plt.plot(solucion1.t, solucion1.y[0], label='S(t)', color = 'mediumblue')
        plt.plot(solucion1.t, solucion1.y[1], label='I(t)', color = 'green')
        plt.plot(solucion1.t, solucion1.y[2], label='R(t)', color = 'red')
        plt.plot(solucion2.t, solucion2.y[0], color = 'mediumblue')
        plt.plot(solucion2.t, solucion2.y[1], color = 'green')
        plt.plot(solucion2.t, solucion2.y[2], color = 'red')
        plt.plot(solucion3.t, solucion3.y[0], color = 'mediumblue')
        plt.plot(solucion3.t, solucion3.y[1], color = 'green')
        plt.plot(solucion3.t, solucion3.y[2], color = 'red')
        plt.axvspan(desde,hasta, color='lavender')
    else:
        return print("Recuerde ingresar un control")
    
    plt.legend()
    plt.grid()
    plt.show()

def plano_fase(parameters, desde = None, hasta = None):

    if desde == None:
        solucion = solve_ivp(SIR, (0,100), [999,1,0], t_eval= np.linspace(0, 100, 1000), args= parameters, method='RK45')
        plt.plot(solucion.y[0], solucion.y[1])


    elif hasta == None:
        solucion = solve_ivp(SIR, (0,100), [999,1,0], t_eval= np.linspace(0, desde, 1000), args= parameters[:-1], method='RK45')
        plt.plot(solucion.y[0], solucion.y[1], color= 'mediumblue')
        z0c=[solucion.y[0][-1],solucion.y[1][-1],solucion.y[2][-1]]
        solucion = solve_ivp(SIR, (desde,100), z0c, t_eval= np.linspace(desde, 100, 1000), args= parameters, method='RK45')
        plt.plot(solucion.y[0], solucion.y[1], color = 'mediumblue')
        plt.axvspan(z0c[0], solucion.y[0][-1], color='lavender')
    else:
        solucion = solve_ivp(SIR, (0,100), [999,1,0], t_eval= np.linspace(0, desde, 1000), args= parameters[:-1], method='RK45')
        plt.plot(solucion.y[0], solucion.y[1], color = 'mediumblue')
        z0c=[solucion.y[0][-1],solucion.y[1][-1],solucion.y[2][-1]]
        solucion = solve_ivp(SIR, (desde, hasta), z0c, t_eval= np.linspace(desde, hasta, 1000), args= parameters, method='RK45')
        plt.plot(solucion.y[0], solucion.y[1], color = 'mediumblue')
        z0n=[solucion.y[0][-1],solucion.y[1][-1],solucion.y[2][-1]]
        solucion = solve_ivp(SIR, (hasta, 100), z0n, t_eval= np.linspace(hasta, 100, 1000), args= parameters[:-1], method='RK45')
        plt.plot(solucion.y[0], solucion.y[1], color = 'mediumblue')
        plt.axvspan(z0c[0], z0n[0], color='lavender')
    plt.xlabel('S')
    plt.ylabel('I')
    plt.ylim(0,1000)
    plt.xlim(0,1000)
    plt.grid()
    plt.show()

beta = 0.003
nu = 1
z0 = [999, 1, 0] 
t_span = (0, 12)
t_eval = np.linspace(t_span[0], t_span[1], 1000)

"""si se desea solo graficar el sistema sin control use"""

# plano_temporal([nu, beta])

"""si desea ver como afecta un control use"""

control = 0.30
# plano_temporal([nu, beta, control])

"""si desea ver como el control afecta en un intervalo de tiempo use"""

# plano_temporal([nu, beta, control], desde=2, hasta=4)




""" si se desea graficar con o sin control use """ 
# plano_fase([nu,beta])
# plano_fase([nu,beta,control])

""" si se desea graficar controlado en un intervalo use """
# plano_fase([nu,beta,control], desde= 2)
# plano_fase([nu,beta,control], desde= 2, hasta= 5)


