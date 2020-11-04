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
plot1=plt.plot(x,y,'*',label='original values')
plot2=plt.plot(x,A_12(x),'r',label='curve_fit')

plt.title('A^(-1/2)')
plt.show()
