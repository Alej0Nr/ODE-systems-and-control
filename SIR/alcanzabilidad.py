import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fmin
from PIL import Image
import os

n=1


for recheability in np.arange(0.1,5.1,0.1):
    recheability= round(recheability,1)
    nu, beta = 1, 0.003
    intervalo_tiempo  = [0,20]
    intervalo_control = [3,3+recheability]

    ## defino t_eval en funcion de los que estan arriba
    t_eval1=np.arange(intervalo_tiempo[0],intervalo_control[0]+0.01,0.01)
    t_eval2=np.arange(intervalo_control[0],intervalo_control[1]+0.01,0.01)

    l = len(t_eval1)



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
        solucion = solve_ivp(SIR, [intervalo_tiempo[0],intervalo_control[0]], [999,1,0], t_eval=t_eval1, method='RK23')
        s, i, r = solucion.y[0][-1],solucion.y[1][-1],solucion.y[2][-1]
        solucion = solve_ivp(SIR, [intervalo_control[0],intervalo_control[1]], [s,i,r], args=(u), t_eval=t_eval2, method='RK23')
        s = solucion.y[0][-1]

        return abs(s-(nu/beta))

    u_opt = fmin(f, x0=0.5, xtol=1e-4, disp=0)
    u_opt = u_opt[0] if isinstance(u_opt, (list, np.ndarray)) else u_opt
    # print(f"Valor óptimo de u: {u_opt}")


    def solucionario(u=0):
        solucion = solve_ivp(SIR, [intervalo_tiempo[0],intervalo_control[0]], [999,1,0], t_eval=t_eval1, method='RK23')
        s1, i1, r1 = solucion.y[0],solucion.y[1],solucion.y[2]
        ini=(s1[-1],i1[-1])
        solucion = solve_ivp(SIR, [intervalo_control[0],intervalo_control[1]], [s1[-1],i1[-1],r1[-1]], args=[u], t_eval=t_eval2, method='RK23')
        s2, i2, r2 = solucion.y[0],solucion.y[1],solucion.y[2]
        s = np.append(s1, s2)
        i = np.append(i1, i2)
        r = np.append(r1, r2)
        return s,i,r,ini



    def grafica_SI(us=[0]):
        for j in range(len(us)):
            s,i,r,ini = solucionario(us[j])
            plt.plot(s[-1], i[-1],'.',color='blue')
            plt.plot(ini[0],ini[1],'.',color="black")
        plt.xlabel('Susceptible')
        plt.ylabel('Infected')
        plt.ylim(0,1000)
        plt.xlim(0,1000)
        plt.grid()
        plt.plot(0,0,color='white',label=f'T = {recheability}')
        plt.axvline(x = nu/beta,color= 'black', linestyle = ':', label = r'$R_0$')
        plt.legend()
        os.makedirs("alcanzabilidad", exist_ok=True)
        global n
        plt.savefig(f'alcanzabilidad/{n}.png')
        n+=1
        plt.clf()
    grafica_SI(np.arange(0,1.001,.001)) #testing

def makegif():
    png_files = []
    for i in range(1,51):
        png_files.append(f"{i}.png")
    output_gif = "alcanzabilidad/alcanzabilidad.gif"
    frame_duration = 100  # milliseconds per frame
    loop = 0  # 0 = infinite loop

    # Load all PNG images
    images = []
    for png_file in png_files:
        try:
            img = Image.open("alcanzabilidad/"+png_file)
            images.append(img)
        except Exception as e:
            print(f"Error loading {png_file}: {e}")

    # Create GIF if we have images
    if images:
        # Save first image and append the rest
        images[0].save(
            output_gif,
            save_all=True,
            append_images=images[1:],
            duration=frame_duration,
            loop=loop,
            optimize=True
        )
        print(f"GIF created successfully: {output_gif}")
        print(f"Size: {images[0].size} (width x height)")
        print(f"Frames: {len(images)}")
        print(f"Duration: {frame_duration}ms per frame")
    else:
        print("No PNG images found. GIF not created.")

# grafica_SI([u_opt,.5,.3,.6,0])

# grafica_SI(np.arange(0,1.001,.001)) #testing
# makegif()