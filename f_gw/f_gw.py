import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

def loadData(filename):
    inFile=open(filename,'r')
    y=[]
    for line in inFile:
        y.append(float(line))
    return y

def func(x,a,b,c):
    return a*x+b*x**2+c*x**3

y=loadData('../A.txt')
length=len(y)
x=np.linspace(1,20,length)
x=x/10
popt,pcov=curve_fit(func,x,y)
a,b,c=popt[0],popt[1],popt[2]

print(a,b,c)

def A_12(z):
    return a*z+b*z**2+c*z**3
def A(z):
    return 0.1449*z-0.0118*z**2+0.0012*z**3
print(A(2.))
plt.plot(x,y,'o',color='red',label='Monte Carlo sampling results')

plt.plot(x,A_12(x),'g',label='fit results')

plt.plot(x,A(x),'b',label='in article')

plt.xlabel("z")
plt.xlim(0.,0.25)
plt.xticks([0,0.4,0.8,1.2,1.6,2.0])
plt.title('A^(-1/2)')
plt.legend()
plt.show()
