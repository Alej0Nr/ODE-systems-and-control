## Ejercicio 1.1.3 Perko
import numpy as np
import matplotlib.pyplot as plt


t = np.arange(-3, 3, .01)
a = -1
x1, y1 = np.exp(t), np.exp(a*t)
c_i = np.arange(-5,5,1) 

## plot diagrama de fase
plt.subplot()

for c1 in c_i:
    for c2 in c_i:
        line, = plt.plot(c1*x1, c2*y1, lw=2, color='purple')
plt.xlim(-5,5)
plt.ylim(-5, 5)
plt.grid(True)
plt.xlabel(f"a={a}")

### plot campo de direcciones
plt.subplot()
x = np.arange(-5,5,1)
y = np.arange(-5,5,1)
X, Y = np.meshgrid(x,y)
dy= a*Y
dx= X
norm = np.sqrt(dx**2 + dy**2)
dyu = dy/norm
dxu = dx/norm
plt.quiver(X,Y,dxu,dyu,color='gray')

plt.xlim(-5,5)
plt.ylim(-5, 5)
plt.grid(True)

plt.show()
