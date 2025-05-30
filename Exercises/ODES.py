import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


class ODE:
    def __init__(self, dxdt, dydt):
        self.dxdt = dxdt
        self.dydt = dydt
        self.t_span = None
        self.t_eval = None

        
    def sistema(self, t, z):
        x, y = z
        dxdt = self.dxdt(x, y)
        dydt = self.dydt(x, y)
        return [dxdt, dydt]
    
    def t_data(self, t_span, n):
        self.t_span = t_span
        self.t_eval = np.linspace(t_span[0], t_span[1], n)

    def campo_direcciones(self, x_arange, y_arange, color = 'blue'):
        X, Y = np.meshgrid(x_arange,y_arange)
        dX = self.dxdt(X, Y)
        dY = self.dydt(X, Y)
        norm = np.sqrt(dX**2 + dY**2)
        norm = np.where(norm == 0, 1e-8, norm)
        dXu = dX/norm
        dYu = dY/norm
        plt.quiver(X,Y,dXu,dYu,color = color)
        plt.xlim(min(x_arange),max(x_arange))
        plt.ylim(min(y_arange),max(y_arange))
        plt.grid(True)
    
    def plano_temporal(self, z0):
        solucion = solve_ivp(self.sistema , self.t_span, z0, t_eval= self.t_eval, method='BDF')
        plt.plot(solucion.t, solucion.y[0])
        plt.plot(solucion.t, solucion.y[1])
        plt.xticks(np.arange(0, self.t_span[1]+1,2))
        plt.grid(True)

    def plano_fase(self, z0s, color = 'mediumblue'):
        for ci in z0s:
            solucion = solve_ivp(self.sistema, self.t_span, ci, method='RK45', rtol= 1e-11)
            plt.plot(solucion.y[0], solucion.y[1], color= color)
        plt.grid(True)


        

def cis(n):
    return [(i,i) for i in range(0,n,1)]


def predador_presa():
    #### Definimos el sistema Predador Presa
    a, b, c, d = 0.4807, 0.02482, 0.02756, 0.9272
    z0 = [30, 4]

    lv = ODE(lambda x, y: x * (c * y - d), lambda x , y: y * (a - b * x))
    lv.t_data((0, 20), 500)

    #### ploteamos distintos graficos
    lv.campo_direcciones(np.arange(0,40,1),np.arange(0,40,1))
    plt.show()

    lv.plano_temporal(z0)
    plt.show()

    lv.plano_fase(cis(20))
    plt.show()


# predador_presa()


def ejemplo1():
    #### definimos dy/dx = x**2 - y**2
    ex = ODE(lambda x,y: 1, lambda x, y: x**2-y**2)
    ex.campo_direcciones(np.arange(-5,6,.5),np.arange(-5,6,.5))
    plt.show()

    ex.t_data((0,2), 500)
    z0 = [.0001,1]

    ex.plano_fase([z0])
    plt.show()
    
# ejemplo1()

