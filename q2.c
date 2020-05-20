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
  for (i=0;i<npts;i++)
  {
    x=xmin+i*delta;			/* Generating sampling for f(x)*/
    if(x!=0)
    {
        sdata[i][REAL]=sin(x)/x;
    }
    else
        sdata[i][REAL]=1; 
    sdata[i][IMAG]=0;
  }
}

int main()
{
	FILE *mptr;
    	mptr=fopen("fourier2.dat","w");		/*opening file to write mode*/
	fftw_complex *sdata;
	fftw_complex *result;		/*initiallzing arrays to store initial and final data*/
	double k,freq;
	double delta=(xmax-xmin)/(npts-1);
	double aft_real,aft_imag;

	sdata  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npts);    /*allocating memory*/
	result = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npts);
	
  	fftw_plan p;    /*algorithm for fftw*/
	
	p = fftw_plan_dft_1d(npts, sdata, result, FFTW_FORWARD, FFTW_ESTIMATE);

	generate_data(sdata,delta);

  	fftw_execute(p); 		/*creates the DFT*/

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
	   aft_real=(1/sqrt(npts))*result[i][REAL];   /*making FT by proper normalisation factor*/
   	   aft_imag=(1/sqrt(npts))*result[i][IMAG];
	   aft_real= delta*sqrt(npts/(2*M_PI))*(cos(freq*xmin)*aft_real+sin(freq*xmin)*aft_imag);
	   aft_imag= delta*sqrt(npts/(2*M_PI))*(cos(freq*xmin)*aft_imag-sin(freq*xmin)*aft_real);
	   fprintf(mptr,"%e\t%e\t%e\n",freq,aft_real,aft_imag);
	  
 	}
	
	fftw_destroy_plan(p);
	fftw_free(sdata); 
	fftw_free(result);
	return (0);
}


