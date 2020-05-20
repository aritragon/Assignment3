import numpy as np
import matplotlib.pyplot as plt

def s (x):
	if(x>=-1 and x<=1):
		return (1)
	else:
		return(0)
def f (x):
	if(x>=-1 and x<=1):
		return (1)
	else:
		return(0)
def f1(x):
	return(np.exp(-x*x))
def s1(x):
	return(np.exp(-4*x*x))

xmax = 10.0
xmin = -10.0

npts = 256
delta = (xmax-xmin)/(npts-1)

xarray = np.zeros(npts)
sdata = np.zeros(npts)
fdata = np.zeros(npts)
nft = np.zeros(npts)

for i in range(npts):
	sdata[i]=s(xmin+i*delta)
	fdata[i]=f(xmin+i*delta)
	xarray[i]=xmin+i*delta
	
	
		
nft1 = (np.fft.fft(sdata,norm='ortho'))
nft2 = (np.fft.fft(fdata,norm='ortho'))

for i in range(npts):
	nft[i] = nft1[i]*nft2[i]

Inft=delta*np.sqrt(npts)*(2/(np.pi))*(np.fft.ifftshift(nft))
plt.plot(xarray,Inft,'red')
plt.plot(xarray,sdata,'blue')
plt.plot(xarray,fdata,'green')
#print(karray.shape)
plt.xlabel('k')
plt.ylabel('f(k)')
plt.show()
