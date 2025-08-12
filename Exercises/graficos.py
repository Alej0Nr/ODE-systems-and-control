import numpy as np
import matplotlib.pyplot as plt
import ODES as ODES
import os

plot    = False
guardar = True
n=1

def guardar(nombre=None):
    os.makedirs("graficos", exist_ok=True)
    global n
    if nombre==None:
        plt.savefig(f'graficos/grafico_{n}.svg')
        n+=1
    else:
        plt.savefig(f'graficos/{nombre}.svg')

def silla(gamma = -1, mu = 1, nombre=None):
    
    silla = ODES.ODE(lambda x, y: gamma*x, lambda x, y: mu*y)
    silla.t_data((0, 20), 500)
    silla.campo_direcciones(np.arange(-10,11,1),np.arange(-10,11,1), color = 'grey')
    silla.plano_fase([[10,1],[10,.5],[-10,1],[-10,.5],[10,-1],[10,-.5],[-10,-1],[-10,-.5]])
    guardar(nombre) if guardar else None
    plt.show() if plot else None
    plt.clf()

def nodo_estable_1(gamma=-1, nombre=None):
    nodo = ODES.ODE(lambda x, y: gamma*x, lambda x, y: gamma*y)
    nodo.t_data((0, 20), 500)
    nodo.campo_direcciones(np.arange(-10,11,1),np.arange(-10,11,1), color = 'grey')
    nodo.plano_fase([[10,10],[-10,-10],[10,-10],[-10,10],[-10,30],[-10,-30],[10,30],[10,-30]])
    guardar(nombre) if guardar else None
    plt.show() if plot else None
    plt.clf()

def nodo_estable_2(gamma = -2, mu = -1, nombre=None):
    nodo = ODES.ODE(lambda x, y: gamma*x, lambda x, y: mu*y)
    nodo.t_data((0, 20), 500)
    nodo.campo_direcciones(np.arange(-10,11,1),np.arange(-10,11,1), color = 'grey')
    nodo.plano_fase([[10,10],[-10,-10],[10,-10],[-10,10],[-10,30],[-10,-30],[10,30],
                     [10,-30],[10,2.5],[10,-2.5],[-10,2.5],[-10,-2.5]])
    guardar(nombre) if guardar else None
    plt.show() if plot else None
    plt.clf()

def nodo_estable_3(gamma = -2, nombre=None):
    nodo = ODES.ODE(lambda x, y: gamma*x+y, lambda x, y: gamma*y)
    nodo.t_data((0, 20), 500)
    nodo.campo_direcciones(np.arange(-10,11,1),np.arange(-10,11,1), color = 'grey')
    nodo.plano_fase([[0,10],[0,-10],[-5,10],[-5,-10],[5,10],[5,-10],[-7.5,10],[-7.5,-10],[7.5,10],[7.5,-10]])
    guardar(nombre) if guardar else None
    plt.show() if plot else None
    plt.clf()

def foco_estable(a = -2, b = 1, nombre=None):
    foco = ODES.ODE(lambda x, y: a*x-b*y, lambda x, y: b*x+a*y)
    foco.t_data((0, 20), 500)
    foco.campo_direcciones(np.arange(-10,11,1),np.arange(-10,11,1), color = 'grey')
    foco.plano_fase([[10,10],[-10,10],[-10,-10],[10,-10],
                     [-4,10],[-4,-10],[4,10],[4,-10],
                     [10,5],[-10,-5]])
    guardar(nombre) if guardar else None
    plt.show() if plot else None
    plt.clf()



def centro(b = 1, nombre=None):
    centro = ODES.ODE(lambda x, y: -b*y, lambda x, y: b*x)
    centro.t_data((0, 20), 500)
    centro.campo_direcciones(np.arange(-10,11,1),np.arange(-10,11,1), color = 'grey')
    centro.plano_fase([[0,1],[0,2],[0,3],[0,4],[0,5]])
    guardar(nombre) if guardar else None
    plt.show() if plot else None
    plt.clf()

def matriz(a, b, c, d, Ci = None, nombre=None):
    """
    el sistema x'=Ax esta dado por
    A=  a b
        c d
    x = x
        y
    """
    matriz = ODES.ODE(lambda x, y: a*x+b*y, lambda x, y: c*x+d*y)
    matriz.t_data((0, 20), 500)
    matriz.campo_direcciones(np.arange(-10,11,1),np.arange(-10,11,1), color = 'grey')
    if Ci == None:
        guardar(nombre) if guardar else None
        plt.show() if plot else None
        plt.clf()
        return
    matriz.plano_fase(Ci)
    guardar(nombre) if guardar else None
    plt.show() if plot else None
    plt.clf()


# silla()
# nodo_estable_1()
# nodo_estable_2()
# nodo_estable_3()
# foco_estable()
# centro()
# matriz(3,-2,1,1)
# matriz(3,1,-1,1)
