import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as  optimize

'''
Nombre Variables y proceso sacado de:
+http://www.phys.lsu.edu/~tohline/PHYS7412/sod.html
+Gas Dynamics: The Riemann Problem and Discontinuous Solutions:
Application to the Shock Tube Problem   http://www.astro.uu.se/~hoefner/astro/teach/ch10.pdf


'''

#contact discontinuity
x_0 = 2.0
#Tiempo de la solucion
t = 1
#puntos
n=5000
#const(gamma)
g=1.4
g2=(g-1)/(2*g)
g1=(g - 1)/(g + 1)
#Condiciones iniciales
denL=1.0
Pl = 1.0
denR = 0.125
Pr = 0.1
U = 0.0

#(10.4)
al= np.sqrt(g*Pl/denL)
#(10.29)

def compability_eq(Ms):
    return (al/g1)*(1-(((Pr/Pl)*((Ms*Ms/g2)- g1))**g2)) - (Ms-(1/Ms))
Ms = optimize.fsolve(compability_eq,0.5)#Sacar raices
#print (Ms)


#(10.28)

U2 = (Ms-(1/Ms))/(g2*g)
a2 = al- U2*g2*g
v_2 = U2 - a2


#(10.30-10.32-10.33)
x1 = x_0 - al*t
x2 = x_0 + v_2*t
x3 = x_0 + U2*t
x4 = x_0 + Ms*t 



x = np.linspace(0, 4, n)
den = np.zeros(n)
P = np.zeros(n)
u = np.zeros(n)
e = np.zeros(n)

#energia de la ecuacion (10.19)

for i in range(0, n):
    #R
    if(x[i]>=x4):
        den[i] = denR
        P[i] = Pr
        u[i] = U
        e[i] = P[i]/(g - 1) + 0.5*den[i]*u[i]**2
    #1
    elif(x[i] >= x3 and x[i] < x4):
        
        den[i] = denR/((2/(g+1))*(1/Ms**2) + g1)
        P[i] = Pl*(a2/al)**(1/g2)
        u[i] = U2
        e[i] = P[i]/(g - 1) + 0.5*den[i]*u[i]**2
    #2
    elif(x[i] >= x2 and x[i] < x3):
        P[i] = Pl*(a2/al)**(1/g2)
        den[i] = denL*(P[i]/Pl)**(1/g)
        u[i] = U2
        e[i] = P[i]/(g - 1) + 0.5*den[i]*u[i]**2
    #E
    elif(x[i] >= x1 and x[i] < x2 ):
        u[i] = (2/(g + 1))*(al + (x[i] - x_0)/t)
        a = al - (g - 1)*u[i]*0.5
        P[i] = Pl*(a/al)**(1/g2)
        den[i] = g*P[i]/a**2
        e[i] = P[i]/(g - 1) + 0.5*den[i]*u[i]**2
    #L
    else:
        den[i] = denL
        P[i] = Pl
        u[i] = U
        e[i] = P[i]/(g - 1) + 0.5*den[i]*u[i]**2
        
#Importacion de datos en archivo
datos = np.loadtxt('UpwindGodunov_step_5.txt')
den_num = datos[:,1]
u_num = datos[:,2]
e_num = datos[:,3]
P_num = datos[:,4]


#Creacion de graficas
plt.figure(1)
plt.plot(x, P,label="Analitica")
plt.plot(x, P_num, label=" Godunov")
plt.title('$Presion P$')
plt.legend()
plt.xlabel(r'$x$')
plt.ylim(0,1.1)
plt.ylabel(r'$P$')
plt.savefig("Presion.png")


plt.figure(2)
plt.plot(x, u,label="Analitica")
plt.plot(x, u_num, label=" Godunov")
plt.title('$Velocidad \vec{u}$')
plt.legend()
plt.xlabel(r'$x$')
plt.ylim(-0.1,1)
plt.ylabel(r'$Velocidad$')
plt.savefig("Velocidad.png")

plt.figure(3)
plt.plot(x, den,label="Analitica")
plt.plot(x, den_num, label=" Godunov")
plt.legend()
plt.title(r'$Densidad \rho$')
plt.xlabel(r'$x$')
plt.ylabel(r'$\rho$')
plt.ylim(0,1.1)
plt.savefig("Densidad.png")



plt.figure(4)
plt.plot(x, e,label="Analitica")
plt.plot(x, e_num, label=" Godunov")
plt.title('$Energia E$')
plt.legend()
plt.xlabel(r'$x$')
plt.ylabel(r'$Energia$')
plt.savefig("Energia.png")





        
