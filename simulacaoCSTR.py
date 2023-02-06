# Importando bibliotecas necessárias
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import fsolve

# Valores Parâmetros do CSTR
F0 = 40                 # ft^3/hr
F = 40                  # ft^3/hr
V = 48                  # ft^3
Fj = 49.9               # ft^3/hr
Vj = 3.85               # ft^3
alpha = 7.08*10**10     # hr^-1
E = 30000               # BTU/mol
R = 1.99                # BTU/mol ºR
U = 150                 # BTU/hr-ft^2-ºR
A = 250                 # ft^2
lambd = -30000          # BTU/mol
Cp = 0.75               # BTU/lbm-ºR
Cj = 1                  # BTU/lbm-ºR
rho = 50                # lmb/ft^3
rhoj = 62.3             # lmb/ft^3


# Condições iniciais
Ca0 = 0.50              # mol/ft^3
Tj0 = 530               # ºR
T0 = 530                # ºR


# Coeficiênte da taxa da reação
def k(T):
    return alpha*np.exp(-E/R/T)
  

# Estado estacionário
def cstr_ss(v):
    Ca, Tj, T = v    # vetor [Ca0, Tj0, T0] com as condições iniciais
    f = np.zeros(3)
    f[0] = F0*(Ca0-Ca)-V*k(T)*Ca
    f[1] = rho*Cp*F0*(T0-T)-lambd*V*k(T)*Ca-U*A*(T-Tj)
    f[2] = rhoj*Cj*Fj*(Tj0-Tj)+U*A*(T-Tj)
    return f    # retorna o sistema de equações do estado estacionário

roots = fsolve(cstr_ss, [Ca0, Tj0, T0])     # função fsolve() encontra as raízes do sistema dadas condições iniciais


print("\nSolução:\n\nCa = " + str(roots[0]) + "\nTj = " + str(roots[1]) + "\nT = " + str(roots[2]))


# Simulção da Dinâmica do Sistema
def cstr_din(v, t):
    Ca, T, Tj = v   # vetor [Ca0, Tj0, T0] com as condições iniciais
    dCa_dt = (F0*(Ca0-Ca)-V*k(T)*Ca)/V
    dT_dt = (rho*Cp*F0*(T0-T)-lambd*V*k(T)*Ca-U*A*(T-Tj))/(rho*Cp*V)
    dTj_dt = (rhoj*Cj*Fj*(Tj0-Tj)+U*A*(T-Tj))/(rhoj*Cj*Vj)
    return [dCa_dt, dT_dt, dTj_dt]

t = np.linspace(0, 10)
v = odeint(cstr_din, [Ca0, Tj0, T0], t)     # função odeint() resolve o sistema de ODEs


fig, axs = plt.subplots(3, 1)
fig.suptitle('Simulação do Modelo do Reator CSTR', fontsize=16)

axs[0].plot(t, v[:, 0], label='Ca', color='orange')
axs[0].set_ylabel('Concentração (Ca)')
axs[0].legend()
axs[0].grid()

axs[1].plot(t, v[:, 1], label='T', color='blue')
axs[1].set_ylabel('Temperatura Interna da Reação (T)')
axs[1].legend()
axs[1].grid()

axs[2].plot(t, v[:, 2], label='Tj', color='green')
axs[2].set_xlabel('tempo (t)')
axs[2].set_ylabel('Temperatura da Jaqueta (Tj)')
axs[2].legend()
axs[2].grid()

plt.show()
