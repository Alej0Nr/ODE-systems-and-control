import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fmin

"""
software implementation of
Sontag - Systems_biology_notes_8.0.11
Chapter 2.5
"""

a, b, c, d = 0.4807, 0.02482, 0.02756, 0.9272
# Condiciones iniciales: (N_0, P_0)
z0 = [30, 4] 
# Intervalo de tiempo (ej: de t=0 a t=20)
t_span = (0, 20)
t_eval = np.linspace(t_span[0], t_span[1], 500)



def sistema(t, z, a, b, c, d):
    n, p = z
    dpdt = p * (c * n - d)
    dndt = n * (a - b * p)
    return [dndt, dpdt]



def lv_campo_direcciones(parameters, points = True):
    a, b, c, d = parameters
    def campo(p_arange,n_arange, color):
        P, N = np.meshgrid(p_arange,n_arange)
        dN= N*(a-b*P)           
        dP= P*(c*N-d)
        norm = np.sqrt(dP**2 + dN**2)
        dNu = dN/norm
        dPu = dP/norm
        plt.quiver(P,N,dPu,dNu,color = color)

    plt.subplot()
    campo(np.arange(-5,40,1),np.arange(-5,40,1),'gray')
    campo(np.arange(0,40,1),np.arange(0,40,1),'blue')
    if points:
        plt.plot(0,0,'r.')
        plt.plot(a/b,d/c,'r.')
    plt.xlim(-5,40)
    plt.ylim(-5,40)
    plt.grid(True)

    plt.show()

#### GPT version 
def lv_campo_direcciones(parameters, points=True):
    a, b, c, d = parameters

    # create full grid
    x = np.arange(-5, 40, 1)
    y = np.arange(-5, 40, 1)
    P, N = np.meshgrid(x, y)

    # compute vector field
    dN = N * (a - b * P)
    dP = P * (c * N - d)
    norm = np.hypot(dP, dN)
    dNu = dN / norm
    dPu = dP / norm

    # mask arrays
    positive_mask = (P >= 0) & (N >= 0)

    # gray where NOT both ≥0
    dPu_gray = np.ma.masked_where(positive_mask, dPu)
    dNu_gray = np.ma.masked_where(positive_mask, dNu)

    # blue where both ≥0
    dPu_blue = np.ma.masked_where(~positive_mask, dPu)
    dNu_blue = np.ma.masked_where(~positive_mask, dNu)

    # plot
    plt.figure()
    plt.quiver(P, N, dPu_gray, dNu_gray, color='gray')
    plt.quiver(P, N, dPu_blue, dNu_blue, color='blue')

    if points:
        plt.plot(0, 0, 'r.')
        plt.plot(a/b, d/c, 'r.')

    plt.xlim(-5, 40)
    plt.ylim(-5, 40)
    plt.grid(True)
    plt.show()

def lv_plano_temporal(parameters):
    solucion = solve_ivp(sistema, t_span, z0, t_eval= t_eval, args= parameters, method='BDF')
    plt.plot(solucion.t, solucion.y[0], label='Presa N(t)')
    plt.plot(solucion.t, solucion.y[1], label='Predador P(t)')
    plt.xlabel('Año')
    plt.ylabel('Población')
    plt.legend()
    plt.xticks(np.arange(0,t_span[1]+1,2))
    plt.grid()
    plt.show()


def lv_plano_fase(parameters, z0, punto = False):
    """
    Toma los parametros del sistema y un vector de condiciones iniciales y plotea
    la solución para cada una de estas
    """
    for ci in z0:
        solucion = solve_ivp(sistema, t_span, ci, args= parameters, t_eval=t_eval, method='BDF', rtol= 1e-11)
        plt.plot(solucion.y[0], solucion.y[1])
    if punto:
        a, b, c, d = parameters
        plt.plot(d/c, a/b, 'x', color= 'black')
    plt.xlabel('Liebres')
    plt.ylabel('Linces')
    plt.title('Diagrama de fase del sistema L-V')
    plt.grid()
    plt.show()


t_eval2 = np.linspace(t_span[0], t_span[1], 21)

def lverr(parameters, Predator, Prey):

    a, b, c, d = parameters
    sis = solve_ivp(sistema, t_span, [30,4], args=(a,b,c,d), t_eval=t_eval2, method='BDF')
    value = (sis.y[1]-Predator)**2 + (sis.y[0]-Prey)**2
    error = sum(value)
    return error

def lv_fit(guess, Predator, Prey):
    parameters = fmin(lverr, guess, args=(Predator, Prey), disp = False)
    sol = solve_ivp(sistema, t_span, z0, args=parameters, t_eval=t_eval2, method='BDF')
    plt.plot(sol.t, sol.y[1],'-,', label='Modelo Lince', color = 'seagreen')
    plt.plot(sol.t, sol.y[0], '-,', label='Modelo Liebre', color = 'red')
    plt.plot(H,'.', color = 'darkorange', label = 'Linces')
    plt.plot(L,'.', color = 'teal', label = 'Liebres')
    plt.xlabel('Año')
    plt.ylabel('Población')
    plt.xticks(np.arange(0,21,2))
    plt.legend()
    plt.grid()
    plt.show()

L = np.array([4, 6.1, 9.8, 35.2, 59.4, 41.7, 19, 13, 8.3, 9.1, 7.4, 8, 12.3, 19.5, 45.7, 51.1, 29.7, 15.8, 9.7, 10.1, 8.6]) # Linces obs
H = np.array([30, 47.2, 70.2, 77.4, 36.3, 20.6, 18.1, 21.4, 22, 25.4, 27.1, 40.3, 57, 76.6, 52.3, 19.5, 11.2, 7.6, 14.6, 16.2, 24.7]) # Liebres obs
guess = (1,.02,.02,1) # punto de partida para empezar a buscar el minimo
#fitted_parameters = fmin(lverr, guess, disp = False, ftol= 1e-6) # a, b, c, d que minimizan el error

def lv_fit_tf(guess, Predator, Prey):
    ## P(t) y N(t)
    plt.subplot(211)
    parameters = fmin(lverr, guess, args= (Predator, Prey), disp = False)
    sol = solve_ivp(sistema, t_span, z0, args=parameters, t_eval=t_eval2, method='BDF')
    plt.plot(sol.t, sol.y[0],'-', label='Modelo Liebre', color = 'red')
    plt.plot(sol.t, sol.y[1], '-', label='Modelo Lince', color = 'seagreen')
    plt.plot(L,'.', color = 'teal', label = 'Linces')
    plt.plot(H,'.', color = 'darkorange', label = 'Liebres')
    plt.xlabel('Año')
    plt.ylabel('Población')
    plt.xticks(np.arange(0,21,2))
    plt.legend()
    plt.grid()
    
    # plano PvsN
    plt.subplot(212)
    solucion = solve_ivp(sistema, t_span, z0, args= parameters, t_eval=t_eval, method='BDF', rtol= 1e-11)
    plt.plot(solucion.y[0], solucion.y[1],color= 'indigo')
    plt.xlabel('Liebres')
    plt.ylabel('Linces')
    plt.grid()
    plt.plot(H, L, '.', color= 'maroon')
    plt.xticks(np.arange(0,81,20))
    plt.yticks(np.arange(0,71,10))
    
    plt.show()

def cis(n):
    return [(d/c , a/b - i) for i in np.linspace(0, a/b, n + 2) if (a/b -i != 0)]


lv_campo_direcciones([a,b,c,d])
# lv_plano_temporal([a,b,c,d])
# lv_plano_fase([a,b,c,d], cis(8), True)
# lv_fit(guess, L, H)
# lv_fit_tf(guess, L, H)
