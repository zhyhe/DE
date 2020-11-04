import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import math

def func(f):
    f0=200.
    S0=1.449e-52
    p1=-4.05
    p2=-0.69
    a1=185.62;a2=232.56
    b1=31.18;b2=-64.72;b3=52.24;b4=-42.16;b5=10.17;b6=11.53
    c1=13.58;c2=-36.46;c3=18.56;c4=27.43
    x=f/f0
    return S0*(x**p1+a1*x**p2+a2*((1+b1*x+b2*x**2+b3*x**3+b4*x**4+b5*x**5+b6*x**6)/(1+c1*x+c2*x**2+c3*x**3+c4*x**4)))

x=np.linspace(1,1500,1500)
plt.axes(yscale='log')
plot2=plt.semilogx(x,np.sqrt(func(x)),'b',label='curve_fit')
print(np.sqrt(func(100.)))
plt.title('Shd')
plt.show()
