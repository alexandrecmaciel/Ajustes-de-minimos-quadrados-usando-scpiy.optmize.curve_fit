# ROTINA PARA RESOLVER NUMERICAMENTE A EQ. DE DIODO SCHOTTKY NA PRESENCA DE RESISTENCIAS PARASITAS
# ALEXANDRE DE CASTRO MACIEL, DEP. DE FISICA DA UFPI
# OUTUBRO DE 2015

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as optm
import __main__
from scipy import special
from matplotlib import font_manager

def JVdiode(vv, Rs, Rp, m, i0):
    kT    = __main__.kT
    log_ii_teo  = []
    for v in vv:
        r      = (v-i0*Rp)/(Rp+Rs)
        a0     = (Rp+Rs)/(i0*Rp)*np.exp(-v/(m*kT))
        c      = Rs/(m*kT)
        z      = c/a0*np.exp(-c*r)
        i      = r+(1/c)*special.lambertw(z)
        if(abs(i.real) > abs(i0)):
            log_ii_teo.append(np.log10(abs(i.real)))
        else:
            log_ii_teo.append(np.log10(i0))
    nlog_ii_teo=np.array(log_ii_teo)
    return nlog_ii_teo
    
# PARAMETROS FIXOS
T        = 330 #in K
k        = 8.617e-5 #in eV/K
kT       = k*T
# PARAMETROS VARIÁVEIS
Rs       = 4.0e2 # in Ohms
Rp       = 8.5e7 # in Ohms
m        = 1.35
i0       = 1.0e-10 # in Amperes
# INICIALIZANDO VETORES

vv  = []
ii_exp  = []
ii_exp_abs  = []
ii_teo  = []
log_ii_exp  = []

# 
directory = "diretorio/para/seu/dado/experimental/" # Directory of files. End with /
path = directory+"nome_do_arquivo.csv" #arquivo em duas colunas separadas por virgula.
look = 0
with open(path, 'r') as data:
    for line in data.readlines():
        if(look):
            vv.append(float(line.split(",")[0]))
            ii_exp.append((float(line.split(",")[1])))
            ii_exp_abs.append(abs((float(line.split(",")[1]))))
            log_ii_exp.append(np.log10(abs(float(line.split(",")[1]))))
        if(line[0:7] == "Voltage"):
            look = 1

nvv = np.array(vv)
nii_exp = np.array(ii_exp)
nii_exp_abs = np.array(ii_exp_abs)
nlog_ii_exp=np.array(log_ii_exp)

# ESCOLHER INTERVALO DE INTERESSE EM TORNO DO CENTRO DO GRAFICO
ran = 0
nvv         = nvv[ran:-1-ran]
nii_exp     = nii_exp[ran:-1-ran]
nii_exp_abs = nii_exp_abs[ran:-1-ran]
nlog_ii_exp = nlog_ii_exp[ran:-1-ran]

x0  = np.array([Rs, Rp, m, i0])
sol, cov = optm.curve_fit(JVdiode, nvv, nlog_ii_exp, x0)

Rs       = sol[0]
Rp       = sol[1]
m        = sol[2]
i0       = sol[3]

nlog_ii_teo = JVdiode(nvv, Rs, Rp, m, i0)

for i in nlog_ii_teo:
    ii_teo.append(10.0**i)

nii_teo = np.array(ii_teo)

font = font_manager.FontProperties(family ='Purisa')
fontx = {'fontname':'Purisa'}

plt.clf()
plt.axvline(0, color='black')
plt.xlabel('Voltagem (V)', **fontx)
plt.ylabel('Corrente (A)', **fontx)
plt.text(0.9*nvv[0], 2e-3, r'$\alpha$ = $%.2f$ ' % m,fontsize=18)
plt.text(0.9*nvv[0], 2e-4, r'$R_P$ = $%.1E \Omega$' % Rp,fontsize=18)
plt.text(0.9*nvv[0], 2e-5, r'$R_S$ = $%.1E \Omega$' % Rs,fontsize=18)
plt.text(0.9*nvv[0], 2e-6, r'$i_0$ = $%.1E A$' % i0,fontsize=18)
plt.semilogy(nvv , nii_exp_abs , 'g-o', label = 'Dados Experimentais')
plt.semilogy(nvv , nii_teo , 'r-', label = 'Ajuste teorico: Diodo', alpha=.5)
plt.legend(loc='lower left', shadow=True, prop = font)
plt.show()