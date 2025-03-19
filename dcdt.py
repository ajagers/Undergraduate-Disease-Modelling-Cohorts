
import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import odeint
import numpy as np

def odes (x, t): 
    # constants (Everything is changed to per min now)
    u_c = 0.5  #1.8 per min    # 0.03 per hr 
    u_p = 0.9      # 0.015 per hr 
    u_i = 1.5         # 0.025

    # Deltas
    S_a = 5.004    # 0.0834  
    S_c = 0.0009     # 0.015 
    S_p = 9e-6   #0.00399995 
    S_i = 3.162    # 0.0527

    alpha_c = 0   #keep curcumin at a range of 0.1 to 0.2

    Beta = 0.0000012 #0.0001597

    r = 1.8    # 0.03

    lambda_p = 18  # 0.3
    lambda_i = 0.3    #0.005
   
    A = x[0]
    P = x[1]
    C = x[2]
    I = x[3]

    # our ODEs 
    dPdt = lambda_p - S_p * P - Beta * C * P   #good 
    dCdt = u_c * (1 - alpha_c) * C  - S_c * (1+ alpha_c) * I * C + Beta * C * P #good
    dAdt = S_p * P + S_c * (1 + alpha_c) * I * C - r * I * A #good
    dIdt = lambda_i - S_i * I + u_c * C #good 

    return [dAdt, dPdt, dCdt, dIdt]

#  initial conditions 
# x_0 = [a_0, p_0, c_0, i_0] 

x_0 = [66.666667, 1.75e6, 40000, 1] 


print(odes(x_0,0))
tf = 4320   # 72 hr to 4320 min
t = np.linspace(0,tf,1000)
x = odeint(odes, x_0, t)

A = x[:,0]
P = x[:,1]
C = x[:,2]
I = x[:,3]

plt.plot(t, C, color = 'blue')
plt.title("$C(t)$ Graph")
plt.xlabel('$t$', fontsize=12)
plt.ylabel('$C(t)$', fontsize=12)
plt.yscale('log')
plt.grid()
plt.show()