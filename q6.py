import numpy as np
import matplotlib.pyplot as plt

def f (x):
	return (1)

xmax = 500.0
xmin = -500.0

npts = 256
delta = (xmax-xmin)/(npts-1)	 #define sample spacing

xarray = np.zeros(npts)
sdata = np.zeros(npts)


for i in range(npts):
	sdata[i]=f(xmin+i*delta)
	xarray[i]=xmin+i*delta
	
		

nft = (np.fft.fft(sdata,norm='ortho'))     #performing DFT
karray = np.fft.fftfreq(npts, d=delta)		#sampling the frequency
karray = 2*np.pi*karray
factor = np.exp(-1j*karray*xmin)
aft = delta*np.sqrt(npts/(2.0*np.pi))*factor*nft	#Numerical fourier transform

k=np.linspace(-8,8,npts)
rkarr=np.asarray(k)
rft= np.zeros(npts)
for i in range(0,npts):
	if(rkarr[i]>=-1 and rkarr[i]<=1):		
		rft[i]=np.sqrt(np.pi/2.0)

#plt.plot(karray,aft,'-*',color='r')
#plt.plot(rkarr,rft)
#print(karray.shape)
#plt.xlabel('k')
#plt.ylabel('f(k)')
#plt.show()
