import numpy as np
import matplotlib.pyplot as plt
import math
import time
def f(x):
	return(np.exp(-x*x))
point1=4
point2=100	
T=np.zeros((point2-point1+1,2))
t=np.zeros((point2-point1+1,1))
e=0
for i in range (point1,point2+1):
	npts =int(i)
	t[e][0]=npts
	sdata=np.zeros(npts)
	ftman= np.zeros(npts)
	xmax = 10.0
	xmin = -10.0
	delta = (xmax-xmin)/(npts-1)
	xarray=np.asarray(np.linspace(xmin,xmax,npts))

	for q in range(npts):
		sdata[q]=f(xarray[q])

	karray = np.fft.fftfreq(npts, d=delta)
	karray = 2*np.pi*karray
	factor = np.exp(-1j*karray*xmin)
	t1=time.time()
	ftpy = delta*np.sqrt(npts/(2.0*np.pi))*factor*np.fft.fft(sdata, norm='ortho')
	T[e][0]=(time.time()-t1)
	t2=time.time()

	for p in range (0,npts):
		for k in range (0,npts):
			ftman[p]=ftman[p]+sdata[k]*np.exp(-1j*karray[p]*xarray[k])		
		ftman[p]=(1/np.sqrt(npts))*ftman[p]
	T[e][1]=(time.time()-t2)
	e=e+1
plt.plot(t,T[:,0],'*',color='r',label="Time for np.fft.fft")
plt.plot(t,T[:,1],'*',color='b',label="manually computing FT")
plt.legend()
plt.xlabel('no. of points')
plt.ylabel('Time')
plt.show()	
