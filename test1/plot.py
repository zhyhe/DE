from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt

import numpy as np

 
def loadData(flieName):
  inFile = open(flieName, 'r')#以只读方式打开某fileName文件
  #定义两个空list，用来存放文件中的数据
  theta = []
  phi = []
  for line in inFile:
    trainingSet = line.split(' ') #对于每一行，按','把数据分开，这里是分成两部分
    theta.append(float(trainingSet[1])) #第一部分，即文件中的第一列数据逐一添加到list theta 中
    phi.append(float(trainingSet[2])) #第二部分，即文件中的第二列数据逐一添加到list phi 中
  return (theta, phi)  # X,y组成一个元组，这样可以通过函数一次性返回

(theta,phi)=loadData('test.txt')
'''
    def randrange(n, vmin, vmax):

        Helper function to make an array of random numbers having shape (n, )

        with each number distributed Uniform(vmin, vmax).

        return (vmax - vmin)*np.random.rand(n) + vmin

    '''
fig = plt.figure()

ax = fig.add_subplot(111, projection='3d')
# draw sphere
u, v = np.mgrid[0:2*np.pi:40j, 0:2*np.pi:40j]
x1 = np.cos(u)*np.sin(v)
y1 = np.sin(u)*np.sin(v)
z1 = np.cos(v)
ax.plot_wireframe(x1, y1, z1, color="0.5",linewidth=0.1)

 

length = len(theta)
counter=1
while counter <= length:
    counter+=1
    xs=np.sin(theta)*np.cos(phi)
    ys=np.sin(theta)*np.sin(phi)
    zs=np.cos(theta)
    ax.scatter(xs, ys, zs, s=2, c='r', marker='.')
 
# For each set of style and range settings, plot n random points in the box

# defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
    '''
for c, m, zlow, zhigh in [('r', 'o', -1, -1), ('b', '^', 1, 1)]:

    xs = randrange(n, 23, 32)

    ys = randrange(n, 0, 100)

    zs = randrange(n, zlow, zhigh)

    ax.scatter(xs, ys, zs, c=c, marker=m)
    ax.scatter(xs, ys, zs, c=c, marker=m)
    '''
 

ax.set_xlabel('X Label')

ax.set_ylabel('Y Label')

ax.set_zlabel('Z Label')

#plt.savefig('sphere.png',dpi=520) 

plt.show()
