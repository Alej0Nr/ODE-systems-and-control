import numpy as np
import matplotlib.pyplot as plt

#logistica discreta
def logistica():
    N=10
    x_0=0.1
    r=2.8
    puntos=[x_0]

    for i in range(N):
        x=puntos[-1]
        puntos.append(r*x*(1-x))

    plt.plot(range(N+1),puntos,'.',color='blue')
    plt.show()

# otra
def ejemplo1():
    N=8

    for ci in [1,3,5,-1,-3,-5,0]:
        puntos=[ci]    
        for i in range(N):
            x=puntos[-1]
            puntos.append(0.5*x)

        plt.plot(range(N+1),puntos,'.',label=r'$x_0=$ '+f'{ci}')
    plt.legend()
    plt.ylim((5.3,-5.3))
    plt.yticks(np.arange(-5,5.1,1))
    plt.grid()
    plt.savefig(f'graficos/PFej1.svg')


logistica()