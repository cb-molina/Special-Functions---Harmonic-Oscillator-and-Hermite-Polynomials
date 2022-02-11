'''
    Christian B Molina
    PHY 104B - WQ2021
    Project 2 - Special Functions Harmonic Oscillator
'''

import numpy as np
import math
from scipy import constants
import matplotlib.pyplot as plt

# Initialization ==================================================
k_order = 4
ndmn = k_order
nary = np.zeros((ndmn+1,ndmn+1))
x = np.linspace(-10,10,1000)
y = np.linspace(-5,5,1000)
zeros = np.zeros((1000))

# nary[n][k]
nary[0][0] = 1
nary[1][0] = 0
nary[1][1] = 2


# Functions =======================================================

## Calculating coefficients
def hermite(korder, a, n_val):
    '''
        Application of equation (4-16)
    '''
    n_val-=1
    coeff = 2*a[n_val][korder-1] - 2*n_val*a[n_val-1][korder]
    return coeff

## Generating the polynomial function
def hermitepoly(korder,n_val,a,x_val):
    ki = 0
    poly = 0
    '''
        Equation (4-11)
        Whatever n_val is upon entering this function, it stays the same 
        korder influences how many times the while loop will iterate. 
    '''
    while ki <= korder:
        poly += a[n_val][ki]*np.power(x_val,ki)
        ki+=1
    
    '''
        When n_val = 0, the return statement results in a zero division.
        To avoid undefined behavior, the return statement will be just have n = 1
    '''
    if n_val == 0:
        return poly / np.power(1,3)
    else:
        return poly / np.power(n_val,3)

## Phi
def harmonicOsci(n_val, x_val):
    '''
        Equation (4-9)
        we're not using their defined rho since it is unitless
    '''
    sArray = np.zeros(len(x_val))
    j = 0
    while j < len(x_val):
        sArray[j] = (1 / (math.sqrt(math.pow(2,n_val)*math.factorial(n_val)*math.sqrt(math.pi)))) * math.exp(-math.pow(x_val[j],2) / 2)
        j += 1
    return sArray


# Coefficient Calculations ========================================
'''
    n will be iterated by one after completely iterating through k-orders
'''
n = 2
while n <= ndmn:
    k = 0
    while k <= k_order:
        nary[n][k] = hermite(k,nary,n)
        k +=1
    n+=1
print(nary)

# Plotting ========================================================
plt.rcParams["font.family"] = "serif"
'''
    Replicating Figure 4-2
    One-dimensional harmonic oscillator wave function phi_n(p) for 
    n = 0, 1, 2, 3 and 4.
'''
fig1 = plt.figure(1)
plt.grid(True)
plt.ylim(-0.7, 0.8)
plt.xlim(-4, 4)
plt.xlabel(r'$\rho$', fontsize = 13, font = "serif")
plt.ylabel(r'$\Psi_n(\rho)$'"\n[Energy]",rotation = 0, fontsize = 13, font = "serif")
plt.title(r'One-dimensional harmonic oscillator wave function $\Psi_n (\rho)$ ', fontsize=12)

i = 0
while i <= k_order:
    plt.plot(x, harmonicOsci(i,x)*hermitepoly(i,i,nary,x))
    plt.legend(('n=0','n=1','n=2','n=3','n=4','n=5'))
    i += 1
'''
    ^ Loop through the plots
'''

plt.plot(x, zeros, "k-.")
plt.plot(zeros,y,"k-.")

'''
    Replicating Figure 4-3
'''
plt.figure(2)
plt.grid(True)
plt.ylim(-1, 3)
plt.xlim(-0.5, 3)


plt.xlabel(r'$\rho$', fontsize = 13, font = "serif")
plt.ylabel(r'$\frac{H_n(\rho)}{n^3}$', rotation = 0, fontsize = 13, font = "serif")
plt.title(r'Hermite polynomials H$_n (\rho)$ for n = 1 to 4 ', fontsize=12)

i = 1
while i <= k_order:
    plt.plot(x,hermitepoly(i,i,nary,x))
    # plt.legend(('n=1','n=2','n=3','n=4','n=5'))
    i += 1
plt.plot(x, zeros, "k-.")
plt.plot(zeros,y,"k-.")



# fig3 = plt.figure(3)
# plt.grid(True)
# # plt.ylim(-0.7, 0.8)
# # plt.xlim(-4, 4)
# plt.xlabel(r'$\rho$', fontsize = 13, font = "serif")
# plt.ylabel(r'$\Psi_n(\rho)$'"\n[Energy]",rotation = 0, fontsize = 13, font = "serif")
# plt.title(r'One-dimensional harmonic oscillator wave function $\Psi_n (\rho)$ ', fontsize=12)


# i_n = 30
# f_n = 40

# plt.plot(x, harmonicOsci(i_n,x)*hermitepoly(i_n,i_n,nary,x), label = f"n={i_n}")
# plt.plot(x, harmonicOsci(f_n,x)*hermitepoly(f_n,f_n,nary,x), label = f"n={f_n}")

# plt.legend()

# plt.plot(x, zeros, "k-.")
# plt.plot(zeros,y,"k-.")


plt.show()


''' Try to increase the excitation level of the harmonic oscillator as far as you can, and explain why there is a limit.
    http://hyperphysics.phy-astr.gsu.edu/hbase/quantum/hosc6.html

    The limit exists at the excitation level of n = 10. At this point, the overall trend of finding the particle is equivalent to
    the finding its classical probability. This phenomenon is an example of the correspondence principle - in which the probability
    of finding the oscillator at some point converges for its classical and quantum probability. 
'''