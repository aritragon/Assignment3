import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
def f (x,y):
	return (np.exp(-x*x-y*y))
xmax = 10.0
xmin = -10.0

ymax = 10.0
ymin = -10.0

npts = 64
deltax = (xmax-xmin)/(npts-1)
deltay = (ymax-ymin)/(npts-1)

xarray = np.zeros(npts)
yarray = np.zeros(npts)
sdata = np.zeros((npts,npts))
#xarray[i]=xmin+i*deltax
#		yarray[]=ymin+i*deltay

for i in range(0,npts):
	for j in range(0,npts):
		sdata[i][j]=f(xmin+i*deltax,ymin+j*deltay)
		
	
		
nft = np.fft.fft2(sdata,norm='ortho')

kxarr = np.fft.fftfreq(npts, d=deltax)
kxarr = 2*np.pi*kxarr
kyarr = np.fft.fftfreq(npts, d=deltay)
kyarr = 2*np.pi*kyarr
factorx = np.exp(-1j*kxarr*xmin)
factory = np.exp(-1j*kyarr*ymin)
#aft = 100*(np.sqrt(2))*deltax*deltay*math.pow(np.sqrt(npts/(2.0*np.pi)),2)*nft
aft = np.zeros((npts,npts))
for p in range (npts):
	for k in range (npts):
		aft[p][k] =deltax*deltay*math.pow(np.sqrt(npts/(2.0*np.pi)),2)*factorx[p]*factory[k]*nft[p][k]
F=np.zeros(npts*npts)
q=0
kx=[]
ky=[]
print(kxarr.shape)

for p in range (npts):
	for k in range (npts):
		F[q]= aft[p][k]
		kx.append(kxarr[p])
		ky.append(kyarr[k])
		q=q+1
		#print(aft[p][k],'\n')
#F1=np.asarray(F)
#print(F)
q=0
x=[]
y=[]
z=np.zeros(100*100)
x1=np.array(np.linspace(-10,10,100))
y1=np.array(np.linspace(-10,10,100))
for p in range (100):
	for k in range (100):
		z[q]=0.5*np.exp(-0.25*(x1[p]*x1[p]+y1[k]*y1[k]))
		x.append(x1[p])
		y.append(y1[k])
		q=q+1
fig=plt.figure()
ax=plt.axes(projection='3d')
ax.scatter(kx,ky,F,color='red',label="Numerical FT")
ax.plot3D(x,y,z,'gray',label="Analytic FT")
ax.set_xlabel("kx")
ax.set_ylabel("ky")
ax.set_zlabel("F(k)")
plt.legend()
plt.show()


