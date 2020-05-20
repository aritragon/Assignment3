import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
myfile=open("noise.dat","r")		# reading the datafile
sdata=myfile.readlines()


npts=int(len(sdata))
for i in range (npts):
	sdata[i]=np.float(sdata[i])  #to make the values of float datatype

P_S=np.zeros(npts)

#f, Pxx_den = signal.periodogram(sdata,return_onesided=False,scaling='spectrum')

nft = (np.fft.fft(sdata,norm='ortho'))
karray = np.fft.fftfreq(npts, d=1)		#sampling the frequency
karray = karray
for i in range (0,npts-1):
	P_S[i]=(1/npts)*np.absolute(nft[i])*np.absolute(nft[i]) #power spectrum using periodogram


fig1=plt.figure()
plt.plot(sdata,'mo') 
plt.title('noise data')

num_bins = 10
fig2=plt.figure()
plt.hist(P_S,num_bins, facecolor='green') #plotting histogram of power spectrum
plt.title('Histogram of Power spectrum')

fig3=plt.figure()
plt.plot(karray,'-r')
#plt.plot(karray,nft,'ro')
plt.xlabel('k')		#plotting the DFT
plt.ylabel('f(k)')
plt.show()
	
fig4=plt.figure()
plt.plot(karray,P_S,'-y')
plt.xlabel('k')			#plotting power spectrum
plt.ylabel('Power spectrum')
plt.show()
