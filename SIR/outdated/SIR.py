import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import os
os.makedirs("alcanzabilidad", exist_ok=True)

n = 1

def SIR(t, z, nu, beta, u=0):
    s, i, r = z
    dsdt = -(1-u)*beta*s*i
    didt = (1-u)*beta*s*i - nu*i
    drdt = nu*i

    return [dsdt, didt, drdt]

def plano_temporal(parameters, desde = None, hasta = None, t_span = (0, 12), z0 = [999,1,0]):
    def ploter(intervalo , z0, control = False, add_label= False):
        sol = solve_ivp(SIR, intervalo, z0, t_eval= np.linspace(intervalo[0], intervalo[1], np.ceil(intervalo[1]-intervalo[0])*10), args= parameters[:2] if control==False else parameters, method='RK45')
        plt.plot(sol.t, sol.y[0], color = 'mediumblue', label = r'$S(t)$' if add_label else None)
        plt.plot(sol.t, sol.y[1], color = 'green', label = r'$I(t)$' if add_label else None)
        plt.plot(sol.t, sol.y[2], color = 'red', label = r'$R(t)$' if add_label else None)
        z0=[sol.y[0][-1],sol.y[1][-1],sol.y[2][-1]]
        return z0
    if desde == None or hasta == None:
        ploter(t_span, z0, add_label= True, control=True)

    elif len(parameters)==3:
        z0 = ploter((t_span[0],desde), z0, add_label= True)
        z0 = ploter((desde,hasta), z0, control= True)
        ploter((hasta, t_span[1]), z0)

        plt.axvspan(desde,hasta, color='lavender')
    else:
        return print("Recuerde ingresar un control")
    plt.xticks(np.arange(t_span[0],t_span[1]+0.1,1))
    plt.axhline(y = nu/beta,color= 'black', linestyle = ':', label = r'$S^\ast$')
    plt.legend()
    plt.grid()
    plt.savefig(f'graficos/{n}.svg')
    plt.show()



def plano_SI(parameters, ci = [999,1,0] , intervalo_tiempo = (0,12), intervalo_control = None, direcciones= False,color = 'blue'):
    """
    Grafica el plano SvsI con o sin intervalo de control, por default usa como condición inicial S0=999, I0=1 y R0=0.
    
    Por default, si se le pasa 3 parametros (nu, beta, control) se grafica el sistema controlado en todo el intervalo,
    si se quiere controlar entre 2 y 4 use intervalo_control= (2,4)

    Ejemplos:
            plano_SI([nu,beta], direcciones=True)
            plano_SI([nu,beta], intervalo_tiempo=(0,5))
            plano_SI([nu,beta,control],direcciones= True)
            plano_SI([nu,beta,control], intervalo_control= (2,4), intervalo_tiempo=(0,8))
            plano_SI([nu,beta,control], intervalo_control= (2,5))
    """
    def ploter(z0, intervalo, control = False):
        solucion = solve_ivp(SIR, intervalo, z0, t_eval= np.linspace(*intervalo, np.ceil((intervalo[1]-intervalo[0]))*10), args= parameters if control else parameters[:2], method='RK45')
        plt.plot(solucion.y[0], solucion.y[1], color= color)
        z0 = [solucion.y[0][-1],solucion.y[1][-1],solucion.y[2][-1]]
        return z0
    
    if intervalo_control == None:
        ploter(ci, intervalo_tiempo, True)

    else:
        z0 = ploter(ci, (intervalo_tiempo[0],intervalo_control[0]))
        z0c = ploter(z0, (intervalo_control[0],intervalo_control[1]), control= True)
        plt.axvspan(z0[0], z0c[0], color='lavender')
        z0 = ploter(z0c, (intervalo_control[1], intervalo_tiempo[1]))
    
    if direcciones:
        u = 0 if len(parameters)==2 else parameters[2]
        nu, beta = parameters[:2]
        S, I = np.meshgrid(np.arange(0,1000, 20),np.arange(0,1000, 20))
        dsdt = -(1-u)*beta*S*I
        didt = (1-u)*beta*S*I - nu*I
        norm= np.hypot(dsdt,didt)
        dSdt = dsdt/norm
        dIdt = didt/norm
        plt.quiver(S,I,dSdt,dIdt, color='lightgrey')
    # plt.plot(parameters[0]/parameters[1],0,'o')
    plt.xlabel('S')
    plt.ylabel('I')
    plt.ylim(0,1000)
    plt.xlim(0,1000)
    plt.grid()







beta = 0.003
nu = 1
z0 = [999, 1, 0] 
"""si se desea solo graficar el sistema sin control use"""

# plano_temporal([nu, beta],t_span=(0,12))

"""si desea ver como afecta un control use"""

# control = 0.34
# plano_temporal([nu, beta, control], t_span=(0,15))

"""si desea ver como el control afecta en un intervalo de tiempo use"""

# plano_temporal([nu, beta, control], desde=2, hasta=5)


""" Para el plano SI """
plano_SI([nu,beta], direcciones=True)
# plano_SI([nu,beta], intervalo_tiempo=(0,5))
# plano_SI([nu,beta,control],direcciones= True)
# control = 0.4103481292724609 #deepseek
# plano_SI([nu,beta,control], intervalo_control= (3,5), intervalo_tiempo=(0,12))
# control = 0.410347747802734
# plano_SI([nu,beta,control], intervalo_control= (3,5), intervalo_tiempo=(0,12),color='red')
# plt.axvline(x = nu/beta, color = 'black', linestyle = ':', label = r'$S^\ast$')
# plt.legend()
# plt.savefig(f'graficos/{n}.svg')
# plt.show()
# plano_SI([nu,beta,control], intervalo_control= (2,5))

for u in np.arange(0,1.01,0.01):
    plano_SI([nu,beta,u], intervalo_control=(2,3), intervalo_tiempo=(0,4))
plt.show()