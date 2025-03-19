
import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import odeint
import numpy as np


def odes (x, t): 
    # constants (Everything is changed to per min now)
    # mu 
    u_c = 1.8      # 1.8    # 0.03 per hr 
    u_p = 0.9      #0.9     # 0.015 per hr 
    u_i = 1.5      #1.5     # 0.025

    # Deltas
    S_a = 5.004    # 5.004  # 0.0834  
    S_c = 0.9      # 0.9    # 0.015 
    S_p = 9e-6     # 9e-6   # 0.00399995 
    S_i = 3.162    #3.162   # 0.0527

    alpha_c = 0    #keep curcumin at a range of 0.1 to 0.2

    Beta = 0.0000012     # 0.00012   # 0.0001597

    r = 1.8             # 1.8     # 0.03

    lambda_p = 18       #18       # 0.3
    lambda_i = 0.3      #0.3      # 0.005
   
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
# apoptotic cells changed from 4000 cells/hr to 66.66667 cells/min

# Problem with changeing apoptotic to cells/min because immune cells go negative 

x_0 = [4000, 1.75e6, 40000, 1] 


print(odes(x_0,0))
tf = 72   # 72 hr to 4320 min
t = np.linspace(0,tf,1000)
x = odeint(odes, x_0, t)

A = x[:,0]
P = x[:,1]
C = x[:,2]
I = x[:,3]

# # # Plotting each equation in a separate subplot with a log y-axis
fig, axis = plt.subplots(2, 2, figsize=(12, 10))

# Plot P(t)
axis[0, 0].plot(t, P, color='blue')
axis[0, 0].set_title("$P(t)$ Graph")
axis[0, 0].set_yscale('log')
axis[0, 0].set_xlabel('$t$')
axis[0, 0].set_ylabel('$P(t)$')
axis[0, 0].grid()

# Plot A(t)
axis[0, 1].plot(t, A, color='green')
axis[0, 1].set_title("$A(t)$ Graph")
axis[0, 1].set_yscale('log')
axis[0, 1].set_xlabel('$t$')
axis[0, 1].set_ylabel('$A(t)$')
axis[0, 1].grid()

# # Plot C(t)

# plt.plot(t, C, color = 'red')
# plt.yscale('log')
# plt.title('$C(t)$ Graph')
# plt.axhline(y=40000, color='black', linestyle='--', linewidth=2, label="40,000 Cancer Cells")

axis[1, 0].plot(t, C, color='red')
axis[1, 0].set_title("$C(t)$ Graph")
axis[1, 0].set_yscale('log')
axis[1, 0].set_xlabel('$t$')
axis[1, 0].set_ylabel('$C(t)$')
axis[1, 0].grid()



axis[1,1].plot(t, I, color='purple')
axis[1,1].set_title('$I(t)$ Graph')
axis[1, 1].set_xlabel('$t$')
# axis[1,1].set_yscale('log')
axis[1,1].set_ylabel('$I(t)$')
axis[1,1].grid()

# Adjust layout for better spacing
plt.tight_layout()
plt.show()
