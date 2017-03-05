import numpy as np
import matplotlib.pyplot as plt
class Expresion:
    def __init__(self, K_r = 1.0, K_p = 60.0, Y_r = 1.0/5.0 , Y_p = 1.0/30.0, r_0 = 0.0, p_0 = 0.0, dt = 0.0003, T_f = 150):
        self.K_r = K_r
        self.K_p = K_p
        self.Y_r = Y_r
        self.Y_p = Y_p
        self.r_0 = r_0
        self.p_0 = p_0
        self.dt = dt
        self.T_f = T_f
        self.tiempo = np.arange(0,T_f,dt)
        self.num_ARNm = []
        self.num_proteinas = []

    def resuelve(self):
        n_1 = 0
        n_2 = 0
        for i in range(0,self.tiempo.size):
            x_1 = np.random.random()
            x_2 = np.random.random()
            x_3 = np.random.random()
            x_4 = np.random.random()
            
            if x_1 < self.K_r*self.dt and x_2 > self.Y_r*n_1*self.dt:
                n_1 = n_1 + 1
                self.num_ARNm.append(n_1)
            elif x_1 > self.K_r*self.dt and x_2 < self.Y_r*n_1*self.dt:
                if n_1 == 0:
                    self.num_ARNm.append(n_1)
                else :
                    n_1 = n_1 - 1
                    self.num_ARNm.append(n_1)
            else :
                self.num_ARNm.append(n_1)
           
            if x_3 < self.K_p*n_1*self.dt and x_4 > self.Y_p*n_2*self.dt:
                n_2 = n_2 + 1
                self.num_proteinas.append(n_2)
            elif x_3 > self.K_p*n_1*self.dt and x_4 < self.Y_p*n_2*self.dt:
                if n_2 == 0:
                    self.num_proteinas.append(n_2)
                else :
                    n_2 = n_2 - 1
                    self.num_proteinas.append(n_2)
            else :
                self.num_proteinas.append(n_2)


    def grafica(self,analitica):
        
        fig = plt.figure()
        ax = plt.axes()
        ax.set_xlabel("$t$")
        ax.set_ylabel("$r(t)$")
        ax.set_title("$r(t)\ vs\ t$ ")
        plt.plot(self.tiempo,self.num_ARNm, label = "$Stochastic\ Solution$" )
        ax.legend(loc='best')

        if analitica == True :
            r = []            
            for i in self.tiempo:
                r.append((self.K_r/self.Y_r)*(1-np.exp(-self.Y_r*i)))
            plt.plot(self.tiempo, r , label = "$Analytic\ Solution$")
            ax.legend(loc='best')
                       
        filename = "r_t"
        plt.savefig( filename + '.pdf', format = 'pdf')
        plt.close()
            
        fig = plt.figure()
        ax = plt.axes()
        ax.set_xlabel("$t$")
        ax.set_ylabel("$p(t)$")
        ax.set_title("$p(t)\ vs\ t$")
        plt.plot(self.tiempo,self.num_proteinas, label = "$Stochastic\ Solution$")
        ax.legend(loc='best')

        if analitica == True :
            p = []
            for i in self.tiempo:
                p.append(((self.K_r*self.K_p)/(self.Y_r*self.Y_p))*(1 - np.exp(-self.Y_p*i) + ((self.Y_p)/(self.Y_p + self.Y_r))*(np.exp(-self.Y_p*i) + np.exp(-self.Y_r*i))))
            plt.plot(self.tiempo, p, label = "$Analytic\ Solution$")
            ax.legend(loc='best')
      
        filename = "p_t"
        plt.savefig(filename + '.pdf', format = 'pdf')
        plt.close()

