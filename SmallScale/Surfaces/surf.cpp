/*
This program generates and analyses 2D random rough surfaces.
The grid size is nxm arbitrarily; the program also allows the user to enter in specific rms slopes which the program can then normalise to.

Compile using ./compile_surf.sh
Run using ./surf.exe < surf.in
*/


#include <math.h>
#include <time.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <ctype.h>


// RANDOM NUMBER GENERATOR
#include "../../Utilities/RandomC/randomc.h"		//enables random number generators (1)
#include "../../Utilities/RandomC/mersenne.cpp"		//Mersenne Twister random number generator (1)
#define RANDOM_GENERATOR CRandomMersenne // defines the generator to use; 'Mersenne Twister'
#include "../../Utilities/RandomC/stocc.h"                     // define random library classes
#include "../../Utilities/RandomC/stoc1.cpp"                   // random library source code
#include "../../Utilities/RandomC/userintf.cpp"                // define system specific user interface

using namespace std;
#include "/usr/local/include/fftw3.h"
#include "../../Utilities/constants.cpp"
#include "../../Utilities/useful_functions.cpp"
#include "../../Utilities/fftw_arr.cpp"
#include "../../Utilities/arr.cpp"

//############# INITIALISES THE RANDOM NUMBER GENERATORS ###############	
	unsigned long int seed=time(0);			// A seed
	StochasticLib1 nrand(seed);		// The generator objects
	CRandomMersenne rg(seed+342);

//void gen_surface(int n1, int n2, double H, double *surf, char *spec_file);
//void spec_check(int n1, int n2, int which, double *surface, double *mod_phase);
void spec_check(int n1, int n2, int which, double *surface, double *mod_phase, double *line, fftw_complex *spec);
void gen_surface2(int n1, int n2, double H, double *surf, char *spec_file);
void arr_write2D_m(int d1, int d2, fftw_complex *arr, char *outfile, int app);
void write2D(int d1, int d2, double H, double d, double Srms, double *arr, char *outfile);
void spec_master(int n1, int n2, double *surface, char *check_file);
void reorder(int n1, int n2, double *surf);
void reorder2(int n1, int n2, double *surf);
double slope_master(int n1, int n2, double *surf, char *slope_file, int s1, int s2, int ds, int smethod, int writeyn);
void norm_surf(int size, double *surface, double measured, double want);
void write_slice(int s1, int s2, int ds, char *file, int n1, int n2, double *surf, double d);
double dD;
int mfl=128;
void write2D_text(int d1, int d2, double H, double d, double Srms, double *arr, char *outfile);

int main()
	{
	int n1, n2,i,j,dim, norm_yn, s1,s2,ds,smethod, check_yn, spec_yn;
	double *surface, H;
	int surf_size;
	ifstream in;
	char surf_file[mfl], check_file[mfl], spec_file[mfl], slope_file[mfl];
	double tanSrms, d, d1slope, sscale;
	
	int slice1, slice2, dslice, sliceyn,opmode;
	char slice_file[mfl];
	
	in.open("surf.in");
	in>>n1>>n2>>H; // n1 and n2 are dimension size; H is the Hurst parameter (0.5=Brownian noise???)
	in>>norm_yn>>d>>sscale>>tanSrms; // d is smallest size of surface, sscale is scale at which the slopes are tanSrms
	in>>opmode>>surf_file; // file to output surface to. Will be over-written.
	cout<<"surf file is: "<<surf_file<<"\n";
	in>>check_yn>>check_file; // file to write check to. Will be over-written.
//	in>>spec_yn>>spec_file; // file to write spectrum to. Will be over-written.
	in>>s1>>s2>>ds>>smethod>>slope_file;
	in>>sliceyn>>slice1>>slice2>>dslice>>slice_file;
	in.close();
	tanSrms=tanSrms*pow(sscale/d,1-H); // alters the expected tanSrms to the appropriate size scale
	
	ofstream out;
	dD=d;
	surf_size=sizeof(double)*n1*n2;
	cout<<"Generating surface, size "<<n1<<"x"<<n2<<" = "<<n1*n2<<" (or "<<surf_size<<" bytes)"<<endl;
	surface=(double*) fftw_malloc(surf_size); // allocates memory to surface
	if (surface ==0) {cout<<"Not enough memory for surface!!!!"<<endl;}
	
	gen_surface2(n1,n2,H,surface,spec_file);
	cout<<"Normalising slopes..."<<endl;
	if (norm_yn==1)
		{
		d1slope=slope_master(n1,n2, surface, slope_file,1,2,1,1,0);
		cout<<d1slope<<"\n";
		norm_surf(n1*n2, surface,d1slope, tanSrms);
		}
	cout<<"Now analysing the slopes.\n";
	if (smethod < 5 && smethod > 0) {d1slope=slope_master(n1,n2, surface, slope_file, s1,s2,ds,smethod,1);}
	if (check_yn == 1) {spec_master(n1,n2,surface,check_file);} // this calculates spectra in individual dimensions
	if (sliceyn > 0) {write_slice(slice1,slice2, dslice, slice_file, n1, n2, surface, d);}
	if (opmode == 1)
		{
		write2D(n1,n2,H,d,tanSrms,surface,surf_file);
		}
	else
		{
		write2D_text(n1,n2,H,d,tanSrms,surface,surf_file);
		}
	}


void write_slice(int s1, int s2, int ds, char *file, int n1, int n2, double *surf, double d)
	{
	int i,j;
	ofstream out;
	
	
	out.open(file);
	out<<"# slices through a 2d "<<n1<<"x"<<n2<<" surface\n";
	for (i=s1; i<s2; i+= ds)
		{
		out<<"#this slice through n1/n2= "<<i<<" respectively.\n";
		for (j=0; j<n2; j++)
			{
			out<<j*d<<" "<<surf[i*n2+j]<<" "<<surf[i+n2*j]<<"\n";
			}
		out<<"\n\n\n\n";
		
		}
	out.close();
	}

void norm_surf(int size, double *surface, double measured, double want)
	{
	int i;
	double factor=want/measured;
	
	for (i=0; i<size; i++)
		{
		surface[i] *= factor;
		}
	}


/// note that the slopes should follow RMS^2 = 2*H-2, typically  1.56-2 = -0.44 for RMS^2, or -0.22 otherwise
double slope_master(int n1, int n2, double *surf, char *slope_file, int s1, int s2, int ds, int smethod, int writeyn)
	{
	int i,j,k;
	double temp;
	ofstream out;
	double slope1, slope2, partn, dmax, d1slope;
	
	if (smethod==4)
		{
		smethod=1;
		ds=1;
		s1=1;
		s2=min(n1,n2);
		}
	
	if (smethod==3)
		{
		s1=1;
		s2=min(n1,n2);
		dmax=(double) s2;
		ds=(int) pow(dmax,0.5);
		smethod=1;
		}
	if (smethod == 2) {ds=max(ds,2); s1=max(s1,1);} // this creates a log-spacing with interval multiple of 2 or more
	if (writeyn == 1) {out.open(slope_file); out<<"#Format: [length scale (m)] [MS slope tangent (x-dir)] [MS slope tangent (y-dir)] \n";}
	
	s2=min(s2,min(n1,n2));
	
	for (k=s1; k<s2; k += 0)
		{
		cout<<"Analysing slopes for k="<<k<<endl;
		slope1=0.;
		slope2=0.;
		for (i=0; i<n1-k; i++)
			{
			for (j=0; j<n2-k; j++)
				{
				slope1 += pow((surf[i*n2+j]-surf[i*n2+j+k])/(k*dD),2);
				slope2 += pow((surf[i*n2+j]-surf[(i+k)*n2+j])/(k*dD),2);
				}
			}
		temp = n2-k;
		temp *= n1-k;
		slope1 /= temp; // just accounts for the number being summed
		slope2 /= temp;
		if (writeyn == 1) {out<<k*dD<<" "<<slope1<<" "<<slope2<<"\n";} // sill squared when written
		slope1 = pow(slope1,0.5);
		slope2 = pow(slope2,0.5);
		if (k == 1) {d1slope=(slope1+slope2)/2.;}
		if (smethod == 2) {k *= ds;} else {k += ds;}
		}
	if (writeyn == 1) {out.close();}
	return d1slope;
	}

// runs a check of the spectrum for 1-D slices along the surface, and writes these to a file
// currently, it only outputs a mean power spectrum in each direction, so things aren't too crazy.
void spec_master(int n1, int n2, double *surface, char *check_file)
	{
	ofstream out;
	int i,j, dim1=(int) (n1/2)+1, dim2=(int) (n2/2)+1;
	double *mod_phase, *line, sum1[dim1], sum2[dim2];
	fftw_plan plan;
	fftw_complex *spec=0;
	
	out.open(check_file);
	
	mod_phase=new double[dim1*2];
	line=(double*) fftw_malloc(sizeof(double)*n1);
	spec = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*dim1);
	plan=fftw_plan_dft_r2c_1d(n1,line,spec,FFTW_ESTIMATE);
	
	out<<"# Dim1 (x): "<<n1<<", Dim2 (y): "<<n2<<"\n";
	out<<"# First outputs spectra in x-direction, then y-direction. Outputs are power, NOT fft!\n";
	out<<"#Format is [frequency (m^-1)] [power - arb units]\n\n\n";
	arr_zero(sum1,dim1);
	set_zero(spec,dim1);
	
	for (i=0; i<n2;  i++)
		{
		for (j=0; j<n1; j++)
			{
			line[j]=surface[i + j*n2];
			}
		fftw_execute(plan);
		arr_ri_to_mp(spec,mod_phase,dim1);
		
		for (j=0; j<dim1; j++)
			{
			sum1[j] += mod_phase[2*j]*mod_phase[2*j]; // adds the magnitude
			}
		}
	// correct normalisation: sum x^2 dx = sum f^2 df, here div by n2 and NX gives sum x^2 = sum f^2. So need
	// sum f^2 = dx/df sum x^2; df = 1/DX = 1/(NX*dx), so 1/df = DX = NX dx, mult is dx^2 NX, but should div by NX do mult by dx^2.
	arr_div(sum1, n2/dD/dD, dim1);
	for (j=0; j<dim1; j++) {out<<j/(dim1*dD)<<" "<<sum1[j]<<"\n";}
	
	delete [] mod_phase;
	fftw_free(line);
	fftw_free(spec);
	fftw_destroy_plan(plan);
	
// ########### now doing it for the other direction ##############

	mod_phase=new double[dim2*2];
	line=(double*) fftw_malloc(sizeof(double)*n2);
	spec = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*dim2);
	plan=fftw_plan_dft_r2c_1d(n2,line,spec,FFTW_ESTIMATE);
	out<<"\n\n\n\n";
	arr_zero(sum2,dim2);
	
	for (i=0; i<n1;  i ++)
		{
		for (j=0; j<n2; j++)
			{
			line[j]=surface[i*n1 + j];
			}
		fftw_execute(plan);
		arr_ri_to_mp(spec,mod_phase,dim1);
		
		for (j=0; j<dim2; j++)
			{
			sum2[j] += mod_phase[2*j]*mod_phase[2*j];
			}
		}
	arr_div(sum2, n1/dD/dD, dim2);
	out<<"\n\n\n";
	for (j=0; j<dim2; j++) {out<<j/(dim2*dD)<<" "<<sum2[j]<<"\n";}
	arr_div(sum2, n1,dim2);
	delete [] mod_phase;
	fftw_free(line);
	fftw_free(spec);
	fftw_destroy_plan(plan);
	
	out.close();
	}

// calculates the spectrum for a single 1-D surface slice
void spec_check(int n1, int n2, int which, double *surface, double *mod_phase, double *line, fftw_complex *spec)
	{
	int d1,d2,i;
	
	if (which < 0)
		{
		d1=n2;
		for (i=0; i<n2; i++)
			{
			line[i]=surface[(which)*n1+i];
			}
		}
	else
		{
		d1=n1;
		for (i=0; i<n1; i++)
			{
			line[i]=surface[which + i*n2];
			}
		}
	d2=(int) (d1/2) +1;
	}

// generates a random rough surface by directly generating the spectrum and transforming.
void gen_surface2(int n1, int n2, double H, double *surf, char *spec_file)//, double delta, double slope_at_delta)
	{
	fftw_plan Plan;
	int d1=n1, d2=(int) (n2/2)+1,i,j,spec_size,i2, hd1=d1/2;
	fftw_complex *spectrum;
	double D, Pwr_slope, exp, expon2, height, amp, phase;
	ofstream out;
	
	D=3-H; // D is the fractal dimension
	Pwr_slope = 2*D-8.; // it should be 2D - 5: I have no idea where 2D-8 comes from...
	
	exp=Pwr_slope/2.; // we do this in the fft domain
	expon2=exp/2.; // we square the positions first
	
	
	spec_size=sizeof(fftw_complex)*d1*d2;
	cout<<"Initialising spectrum, size "<<d1<<"x"<<d2<<" = "<<d1*d2<<" (or "<<spec_size<<" bytes) ..."<<endl;
	spectrum=(fftw_complex*) fftw_malloc(spec_size);
	
	Plan=fftw_plan_dft_c2r_2d(n1,n2,spectrum,surf,FFTW_ESTIMATE);
	
/* generates a random spectrum.
The sum runs from both +ve and -ve frequencies in dim1, but over +ve frequencies only in dim2.
Note that the array set-up gives a strange arrangement of frequencies.
*/
	
	for (i=0; i<d1; i++)
		{
		if (i >= hd1) {i2=d1-i;} else {i2=i;}
		for (j=0; j<d2; j++)
			{
			if (i2==0 && j==0) {height=0;} else {
			height = pow(i2*i2+j*j,expon2);} // generates +ve frequencies only
			
			amp = nrand.Normal(0,height);
			phase = twopi*rg.Random();
			
			spectrum[i*d2+j][0] = amp*cos(phase);
			spectrum[i*d2+j][1] = amp*sin(phase);
			}
		}
	cout<<"    ... finished analysing spectrum.\nExecuting fourier transform..."<<endl;
	fftw_execute(Plan);
	cout<<" ... finished execution. Reordering the array to something sensible and writing surface..."<<endl;
	reorder(n1,n2,surf); // the surface outputs aren't in the correct order - this is just the way fftw runs
	reorder2(n1,n2,surf); // flips things in the other direction now, to *finally* get the correct order.
	arr_write2D_m(d1,d2,spectrum,spec_file,0); // writes the surface
	cout<<"Done. The surface is generated!"<<endl;
	fftw_free(spectrum);
	}

void reorder(int n1, int n2, double *surf)
	{
	int i, h=n1*n2/2;
	double *temp;
	temp = new double[h];
	for (i=0; i<h; i++) {temp[i]=surf[i+h];}
	for (i=0; i<h; i++) {surf[i+h]=surf[i];}
	for (i=0; i<h; i++) {surf[i]=temp[i];}
	delete [] temp;
	}

void reorder2(int n1, int n2, double *surf)
	{
	long int i,j, hn2;
	double *temp;
	temp = new double[(long int) (n1*n2)];
	hn2=n2/2;
	for (i=0; i<n1*n2; i++) {temp[i]=surf[i];}
	for (i=0; i<n1; i++)
		{
		for (j=0; j<hn2; j++)
			{
			surf[i*n1+j]=temp[i*n1+j+hn2];
			}
		for (j=0; j<hn2; j++)
			{
			surf[i*n1+j+hn2]=temp[i*n1+j];
			}
		}
	delete [] temp; temp=0;
	}

void arr_write2D_m(int d1, int d2, fftw_complex *arr, char *outfile, int app)
	{
	ofstream out;
	int i,j;
	
	if (app == 0)
		{
		out.open(outfile);
		}
	else
		{
		out.open(outfile, ios::app);
		out<<"\n\n\n\n";
		}
	out<<"# "<<d1<<" "<<d2<<"\n";
	
	
	for (i=0; i<d1; i++)
		{
		for (j=0; j<d2; j++)
			{
			out<<pow(arr[i*d2+j][0]*arr[i*d2+j][0]+arr[i*d2+j][1]*arr[i*d2+j][1],0.5)<<" ";
			}
		out<<"\n";
		}
	out<<"\n\n\n\n";
	out.close();
	}


void write2D(int d1, int d2, double H, double d, double Srms, double *arr, char *outfile)
	{
	ofstream out;
	float temp;
	int i,j;
	
	out.open(outfile, ios::binary);
	out.write((char *) &d1, sizeof(int));
	out.write((char *) &d2, sizeof(int));
	temp = (float) H;
	out.write((char *) &temp, sizeof(float));
	temp = (float) d;
	out.write((char *) &temp, sizeof(float));
	temp = (float) Srms;
	out.write((char *) &temp, sizeof(float));
	
//	out<<"# "<<d1<<" "<<d2<<" "<<H<<" "<<d<<" "<<Srms<<"\n";
	
	for (i=0; i<d1; i++)
		{
		for (j=0; j<d2; j++)
			{
			temp=(float) arr[i*d2+j];
			out.write((char *) &temp, sizeof(float));
			}
		}
	out.close();
	
	}


void write2D_text(int d1, int d2, double H, double d, double Srms, double *arr, char *outfile)
	{
	ofstream out;
	float temp;
	int i,j;
	
	out.open(outfile);
	out<<"# "<<d1<<" "<<d2<<" "<<H<<" "<<d<<" "<<Srms<<"\n";
	
	for (i=0; i<d1; i++)
		{
		for (j=0; j<d2; j++)
			{
			out<<arr[i*d2+j]<<" ";
			}
		out<<"\n";
		}
	out.close();
	
	}
