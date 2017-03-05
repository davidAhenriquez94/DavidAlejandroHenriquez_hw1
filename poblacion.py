import genetica
import matplotlib.pyplot as plt
import numpy as np 
import scipy.special as ssp
exp = genetica.Expresion()
fileout_1 = open('r_poblacion.dat','w')
fileout_2 = open('p_poblacion.dat','w')
r = []
p = []
for i in range(0,200):
    exp.resuelve()
    x = exp.num_ARNm[-1]
    y = exp.num_proteinas[-1]
    fileout_1.write(str(x))
    fileout_2.write(str(y))
    r.append(x)
    p.append(y)

fig = plt.figure()
ax = plt.axes()
ax.set_xlabel('numero de ARNm')
ax.set_ylabel('frecuencia')
ax.set_title('histograma ARNm')
plt.hist(r)

z = []
for i in np.linspace(0,10,num=40):
    y = exp.K_r/exp.Y_r
    x = 100*((y**i)*np.exp(-y))/ssp.gamma(i+1)
    z.append(x)

plt.plot(np.linspace(0,12,num=40),z)

filename = "r_histograma"
plt.savefig(filename + '.pdf', format = 'pdf')
plt.close()

fig = plt.figure()
ax = plt.axes()
ax.set_xlabel('numero de proteinas')
ax.set_ylabel('frecuencia')
ax.set_title('histograma proteinas')
plt.hist(p)

filename = "p_histograma"
plt.savefig(filename + '.pdf', format = 'pdf')
plt.close()
