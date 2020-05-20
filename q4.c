#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fftw3.h>
#define xmax 50.0
#define xmin -50.0
#define npts 256
#define REAL 0
#define IMAG 1


void generate_data (fftw_complex* sdata,double delta)
{
  int i;
  double x;
  for (i=0;i<npts;i++)			/* Generating sampling for f(x)*/
  {
    x=xmin+i*delta;
    sdata[i][REAL]=exp(-x*x);
    sdata[i][IMAG]=0;
  }
}
double FT( double x)
{

	return((1/sqrt(2))*exp(-0.25*x*x));	
}
int main()
{
	FILE *mptr;
    	mptr=fopen("fourier4.dat","w");		/*opening file to write mode*/
	
	fftw_complex *sdata;
	fftw_complex *result;		/*initiallzing arrays to store initial and final data*/
	sdata  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npts);	 /*allocating memory*/
	result = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npts);
	double k,freq;
	double delta=(xmax-xmin)/(npts-1);
	double aft_real,aft_imag,rft;
  	fftw_plan p;  	/*algorithm for fftw*/
	
	p = fftw_plan_dft_1d(npts, sdata, result, FFTW_FORWARD, FFTW_ESTIMATE);

	generate_data(sdata,delta)		/*algorithm for fftw*/
  	fftw_execute(p); 	/*creates DFT*/

	int i;
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
	   aft_real=(1/sqrt(npts))*result[i][REAL];
   	   aft_imag=(1/sqrt(npts))*result[i][IMAG];    /*creating FT using proper normalisation*/
	   aft_real= delta*sqrt(npts/(2*M_PI))*(cos(freq*xmin)*aft_real+sin(freq*xmin)*aft_imag);
	   aft_imag= delta*sqrt(npts/(2*M_PI))*(cos(freq*xmin)*aft_imag-sin(freq*xmin)*aft_real);
	
	   rft=FT(freq);
	   fprintf(mptr,"%e\t%e\t%e\t%e\n",freq,aft_real,aft_imag,rft);  /*writing data in a file*/
	  
 	}
	
	fftw_destroy_plan(p);
	fftw_free(sdata); 
	fftw_free(result);
	return (0);
}


