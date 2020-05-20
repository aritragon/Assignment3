#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define xmax 50.0
#define xmin -50.0
#define npts 256
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

double generate_data (int i,double delta)
{
  double x;

    x=xmin+i*delta;		/* Generating sampling for f(x)*/
    if(x!=0)
    {
        return(sin(x)/x);
    }
    else
        return (1.0);
}

int main ()
{
  FILE *mptr;
  mptr=fopen("fourier3.dat","w");	/*opening file to write mode*/

  int i; 
  double sdata[2*npts];
  double freq;
  double delta=(xmax-xmin)/(npts-1);
  double aft_real,aft_imag;

  for (i = 0; i < npts; i++)
    {
       REAL(sdata,i) = generate_data(i,delta); 
       IMAG(sdata,i) = 0.0;
    }

  gsl_fft_complex_radix2_forward (sdata, 1, npts);	 /*algorithm for gsl*/

 for (i=0;i<npts;i++)
  	{
	   if (i<=npts/2-1)
	   {
		freq=(2*M_PI/(npts*delta))*(i);
	   }
	   else
	   {
		freq=(2*M_PI/(npts*delta))*(i-npts);
	   }
	   aft_real=(1/sqrt(npts))*REAL(sdata,i);
   	   aft_imag=(1/sqrt(npts))*IMAG(sdata,i);
	   aft_real= delta*sqrt(npts/(2*M_PI))*(cos(freq*xmin)*aft_real+sin(freq*xmin)*aft_imag);
	   aft_imag= delta*sqrt(npts/(2*M_PI))*(cos(freq*xmin)*aft_imag-sin(freq*xmin)*aft_real);
	   fprintf(mptr,"%e\t%e\t%e\n",freq,aft_real,aft_imag); 	/*writing data in file*/
	  
 	}

  return 0;
}
