## Ejercicio 1.1.3 Perko
import numpy as np
import matplotlib.pyplot as plt


t = np.arange(-3, 3, .01)

x1, y1 = np.exp(t), np.exp(-t)

x2, y2 = np.exp(t), np.exp(0*t)

x3, y3 = np.exp(t), np.exp(t)

Condiciones_Iniciales = [-5,-4,-3,-2,-1,0,1,2,3,4,5]


plt.subplot(221)

for c1 in Condiciones_Iniciales:
    for c2 in Condiciones_Iniciales:
        line, = plt.plot(c1*x1, c2*y1, lw=2, color='b')
plt.xlim(-5,5)
plt.ylim(-5, 5)
plt.grid(True)
plt.xlabel("a=-1")

plt.subplot(222)
for c1 in Condiciones_Iniciales:
    for c2 in Condiciones_Iniciales:
        line, = plt.plot(c1*x2, c2*y2, lw=2, color='r')
plt.xlim(-5,5)
plt.ylim(-5, 5)
plt.grid(True)
plt.xlabel("a=0")

plt.subplot(223)
for c1 in Condiciones_Iniciales:
    for c2 in Condiciones_Iniciales:
        line, = plt.plot(c1*x3, c2*y3, lw=2, color='green')
plt.xlabel("a=1")

plt.xlim(-5,5)
plt.ylim(-5, 5)
plt.grid(True)
plt.show()
