import numpy as np
import matplotlib.pyplot as plt

def f (x):
	if(x!=0):
		return (np.sin(x)/x)  #define f(x)
	else:
		return(0.0)
def ftilde (x):
	if(abs(x)<=1):
		return (np.sqrt(np.pi/2.0)) #define analytic fourier transform
	else:
		return(0.0)
####################################################################
file1=open("fourier2.dat","r")
file2=open("fourier3.dat","r")
lines1=file1.readlines()
lines2=file2.readlines()
file1.close()
file2.close()
freq=[]
cft=[]
freq2=[]
cft2=[]
for line in lines1[1:]:			#To match with c code fourier transform using fftw
	p=line.split()
	freq.append(float(p[0]))
	cft.append(float(p[1]))
freq=np.asarray(freq)
cft=np.asarray(cft)

for line in lines2[1:]:			#To match with c code fourier transform using gsl
	p=line.split()
	freq2.append(float(p[0]))
	cft2.append(float(p[1]))
freq2=np.asarray(freq2)
cft2=np.asarray(cft2)
#####################################################################

xmax = 50.0
xmin = -50.0

npts = 256
delta = (xmax-xmin)/(npts-1)  #define sample spacing

xarray = np.zeros(npts)
sdata = np.zeros(npts)
rft= np.zeros(npts)

for i in range(npts):
	sdata[i]=f(xmin+i*delta)
	xarray[i]=xmin+i*delta
	
		
nft = (np.fft.fft(sdata,norm='ortho'))    #performing DFT
karray = np.fft.fftfreq(npts, d=delta)		#sampling the frequency
karray = 2*np.pi*karray
factor = np.exp(-1j*karray*xmin)
aft = delta*np.sqrt(npts/(2.0*np.pi))*factor*nft	#Numerical fourier transform

k=np.linspace(-8,8,npts)
rkarr=np.asarray(k)


for i in range(0,npts):
	rft[i]=ftilde(karray[i])

#plt.plot(karray,aft,'-b^',markersize='3.5',label='python_numerical')
#plt.plot(freq,cft,'-m*',markersize='2.5',label='c_numerical_fftw')
plt.plot(freq2,cft2,'-go',markersize='1.5',label='c_numerical_gsl')
plt.plot(karray,rft,'-r',markersize='0.5',label='analytic')	#Plotting

plt.xlabel('k')
plt.ylabel('f(k)')
plt.legend()
plt.show()
