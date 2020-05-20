import numpy as np
import matplotlib.pyplot as plt

def s (x):
	if(x>=-1 and x<=1):
		return (1)		#box fuction
	else:
		return(0)
def f (x):
	if(x>=-1 and x<=1):
		return (1)
	else:
		return(0)

xmax = 10.0
xmin = -10.0

npts = 256
delta = (xmax-xmin)/(npts-1)  #setting sampling frequency

xarray = np.zeros(npts)
sdata = np.zeros(npts)
fdata = np.zeros(npts)

for i in range(npts):
	sdata[i]=s(xmin+i*delta)
	fdata[i]=f(xmin+i*delta)
	xarray[i]=xmin+i*delta
	
Inft=delta*np.convolve(sdata,fdata, 'same')   #convolution command in numpy	
		

plt.plot(xarray,Inft,'red',label="convolution f(x)*g(x)")
plt.plot(xarray,sdata,'blue',label="box fuction")		#plotting
plt.xlabel('x')
plt.legend()
plt.show()
