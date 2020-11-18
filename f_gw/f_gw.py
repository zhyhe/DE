import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

def loadData(filename):
    inFile=open(filename,'r')
    x=[]
    y=[]
    for line in inFile:
        deltad_L = line.split(',')
        x.append(float(deltad_L[0]))
        y.append(float(deltad_L[1]))
    return (x,y)

def func(x,a,b,c):
    return a*x+b*x**2+c*x**3

(x,y)=loadData('../A.txt')
length=len(x)
popt,pcov=curve_fit(func,x,y)
a,b,c=popt[0],popt[1],popt[2]

print(a,b,c)

def A_12(z):
    return a*z+b*z**2+c*z**3
def A(z):
    return 0.1449*z-0.0118*z**2+0.0012*z**3
print(A(2.))
plt.semilogy(x,y,'.',color='red',label='Monte Carlo sampling results')

#plt.plot(x,A_12(x),'g',label='fit results')

#plt.plot(x,A(x),'b',label='in article')

plt.xlabel("z")
plt.xlim(0.,0.25)
plt.xticks([0,0.4,0.8,1.2,1.6,2.0])
plt.title('A^(-1/2)')
plt.legend()
plt.show()
