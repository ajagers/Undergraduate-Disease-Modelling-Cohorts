import math
import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import numpy as np

# Constants
# u_ is natural growth/decay rate 
u_c = 0.1
u_a = 1
u_p = 0.03

# machine learn ofr alpha_c so experiment 
alpha_c = 0.5
beta = 159.7

# Natural Decay 
S_a = 0.015

c_0= 0
a_0 = 0
p_0= 17.5 * (10 ** 6)

# percentage inhibition from Eliza 
p=0.23

t = np.linspace(0,100,1000)

def dCdt(C, t):
    return (1-p)*u_c * C - (S_a + alpha_c) * a_0 + (beta * p_0)

y0 = [c_0]

solve_C = odeint(dCdt, y0, t=t)
C_sol = solve_C.T[0]


plt.plot(t, C_sol, color = 'blue')
plt.title("$C(t)$ Graph with Constant A and P")
plt.xlabel('$t$', fontsize=12)
plt.ylabel('$C(t)$', fontsize=12)
plt.yscale('log')
plt.grid()
plt.show()