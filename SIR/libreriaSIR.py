import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def SIR(t, z, nu=1, beta=0.003, u=0):
    s, i, r = z
    dsdt = -(1-u)*beta*s*i
    didt = (1-u)*beta*s*i - nu*i
    drdt = nu*i

    return [dsdt, didt, drdt]

class SIRsim:
    def __init__(self, Ci, tiempo_sim, tiempo_con=None, params=[1,0.003], control=0):
        self.ci = Ci                        # condiciones iniciales
        self.tiempo_sim = tiempo_sim        # tiempo que se desea simular
        self.params = params                # parametros del sistema
        self.R_star = params[0]/params[1]   # calculamos y guardamos R estrella
        self.tiempo_con = tiempo_con        # tiempo en el que se va a controlar
        self.control = control              # control que se va a aplicar
        self.f = None                       # guardara el estado de S, I y R al final de la simulacion
        self.dc = None                      # guardara el estado de S, I y R despues de controlar
        self.t = []                         # guardara vector de tiempo     
        self.S = []                         # guardara S(t)
        self.I = []                         # guardara I(t)
        self.R = []                         # guardara R(t)
        


    def simulacion(self, rtol=1e-7, method='RK45'):
        self.f, self.dc = [None]*2
        self.t, self.S, self.I, self.R = [[]]*4

        if self.tiempo_con==None:
            solucion = solve_ivp(SIR, self.tiempo_sim, self.ci, t_eval= None, rtol= rtol,args= self.params , method= method)
            self.t = solucion.t
            self.S, self.I, self.R = solucion.y[0],solucion.y[1],solucion.y[2]
            self.f = self.S[-1], self.I[-1], self.R[-1]
        else:
            solucion = solve_ivp(SIR, (self.tiempo_sim[0],self.tiempo_con[0]), self.ci, t_eval= None, rtol= rtol,args= self.params , method= method)
            self.t = np.append(self.t, solucion.t)
            self.S = np.append(self.S, solucion.y[0])
            self.I = np.append(self.I, solucion.y[1])
            self.R = np.append(self.R, solucion.y[2])
            solucion = solve_ivp(SIR, (self.tiempo_con[0], self.tiempo_con[1]), [self.S[-1],self.I[-1],self.R[-1]], t_eval= None, rtol= rtol,args= (*self.params,self.control) , method= method)
            self.t = np.append(self.t, solucion.t)
            self.S = np.append(self.S, solucion.y[0])
            self.I = np.append(self.I, solucion.y[1])
            self.R = np.append(self.R, solucion.y[2])
            self.dc = self.S[-1],self.I[-1],self.R[-1]
            solucion = solve_ivp(SIR, (self.tiempo_con[1], self.tiempo_sim[1]), [self.S[-1],self.I[-1],self.R[-1]], t_eval= None, rtol= rtol,args= self.params , method= method)
            self.t = np.append(self.t, solucion.t)
            self.S = np.append(self.S, solucion.y[0])
            self.I = np.append(self.I, solucion.y[1])
            self.R = np.append(self.R, solucion.y[2])
            self.f = self.S[-1], self.I[-1], self.R[-1]

    def cambio_control(self, control):
        self.control = control

    def grafica_temporal(self, colores=['mediumblue','red','green']):
        plt.plot(self.t,self.S,label=r'$S(t)$', color = colores[0])
        plt.plot(self.t,self.I,label=r'$I(t)$', color = colores[1])
        plt.plot(self.t,self.R,label=r'$R(t)$', color = colores[2])
        plt.xticks(np.arange(*self.tiempo_sim,1))
        
    
    def grafica_SI(self, color = 'purple', R_star=True):
        plt.axvline(x = self.R_star,color= 'black', linestyle = ':', label = r'$R^\ast$') if R_star == True else None
        plt.plot(self.S,self.I,color = color, label=r'$SvsI$')
        plt.ylim(0,sum(self.ci))
        plt.xlim(0,sum(self.ci))

    def show(self):
        plt.grid()
        plt.legend()
        plt.show()
        plt.clf()


if __name__=='__main__':
    sir = SIRsim(Ci=(999,1,0), tiempo_sim=(0,16), tiempo_con=(3,5), params=[1,0.003], control=.60)
    sir.simulacion()

    print(sir.dc)   # estado despues de controlar
    print(sir.f)    # estado final
    sir.grafica_temporal()
    sir.show()

    sir.cambio_control(control=0.4938655741319217)
    sir.simulacion()

    print(sir.dc)   # estado despues de controlar
    print(sir.f)    # estado final
    sir.grafica_temporal()
    sir.show()

    sir.grafica_SI()
    sir.show()