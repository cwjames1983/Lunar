void set_zero(fftw_complex *arr, int n);

void arr_add_c2c(fftw_complex *from, fftw_complex *to, int n);

void arr_add_m2d(fftw_complex *from, double *to, int n);
void arr_add_sqr2d(fftw_complex *from, double *to, int n);
void arr_ri2mod2(fftw_complex *from, double *to, int n);
void arr_ri2mod(fftw_complex *from, double *to, int n);
void arr_zero_i(fftw_complex from[], int n);

void arr_add_r2d(fftw_complex *from, double *to, int n);
void arr_sub_d2r(double *sub, fftw_complex *from, int n);
void arr_print_r(fftw_complex *arr, int n);
double arr_rms_r(fftw_complex *arr, int n, double mean);
double arr_mean_r(fftw_complex *arr, int n);
void arr_czero(fftw_complex *arr, int n);
void arr_div_c(fftw_complex *arr, int n, double div);
void arr_mul_c(fftw_complex *arr, int n, double mul);
void arr_copy_d2r(double *d, fftw_complex *c, int n);
void arr_copy_c2c(fftw_complex *from, fftw_complex *to, int n);
void arr_pow_r(fftw_complex *c, double pwr, int n);
void arr_add_mag2d(fftw_complex *from, double *to, int n);
void arr_ri_to_mp(fftw_complex *from, double *to, int n);
void multiply(fftw_complex fromone, fftw_complex fromtwo, fftw_complex to, int n);

void arr_write_r(fftw_complex *arr, int n);
void arr_writel_rc(fftw_complex *arr, int n, char *filename);
void arr_writel2d_m(fftw_complex *arr, int n1, int n2,  char *filename, int old);
void arr_writel2d_m_grid(fftw_complex *arr, int n1, int n2,  char *filename, int old);
void arr_writel_m(fftw_complex *arr, int n, char *filename);
void arr_writel_rc(fftw_complex *arr, int n, char *filename, int old);
void arr_writel_m(fftw_complex *arr, int n, char *filename, int old);
void arr_cmultiply(fftw_complex *one, fftw_complex *two, fftw_complex *three, int n);
void arr_writel2d_m(fftw_complex *arr, int n1, int n2,  char *filename);

void set_zero(fftw_complex *arr, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		arr[i][0]=0.;
		arr[i][1]=0.;
		}
	}

void arr_ri_to_mp(fftw_complex *from, double *to, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		to[2*i]=pow(from[i][0]*from[i][0]+from[i][1]*from[i][1],0.5);
		to[2*i+1]=my_atan(from[i][0],from[i][1]);
		}
	}
	
void arr_add_mag2d(fftw_complex *from, double *to, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		to[i] += from[i][0]*from[i][0]+from[i][1]*from[i][1];
		}
	}

void arr_pow_r(fftw_complex *c, double pwr, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		c[i][0]=pow(c[i][0],pwr);
		}
	}

void arr_copy_d2r(double *d, fftw_complex *c, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		c[i][0]=d[i];
		}
	}

void arr_copy_c2c(fftw_complex *from, fftw_complex *to, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		to[i][0]=from[i][0];
		to[i][1]=from[i][1];
		}
	}

void arr_div_c(fftw_complex *arr, int n, double div)
	{
	int i;
	for (i=0; i<n; i++)
		{
		arr[i][1]/= div;
		arr[i][0]/=div;
		}
	}
	
void arr_mul_c(fftw_complex *arr, int n, double mul)
	{
	int i;
	for (i=0; i<n; i++)
		{
		arr[i][1]*= mul;
		arr[i][0]*=mul;
		}
	}
	
void arr_czero(fftw_complex *arr, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		arr[i][1]=0.;
		arr[i][0]=0.;
		}
	}

void arr_zero_i(fftw_complex from[], int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		from[i][1]=0;
		}
	}

void arr_add_c2c(fftw_complex *from, fftw_complex *to, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		to[i][0] += from[i][0];
		to[i][1] += from[i][1];
		}
	}

void arr_add_sqr2d(fftw_complex *from, double *to, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		to[i] += from[i][0]*from[i][0]+from[i][1]*from[i][1];
		}
	}

void arr_ri2mod(fftw_complex *from, double *to, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		to[i] = from[i][0]*from[i][0]+from[i][1]*from[i][1];
		to[i]=pow(to[i],0.5);
		}
	}
void arr_ri2mod2(fftw_complex *from, double *to, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		to[i] = from[i][0]*from[i][0]+from[i][1]*from[i][1];
		}
	}
	
void arr_add_m2d(fftw_complex from[], double *to, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		to[i] += pow(from[i][0]*from[i][0]+from[i][1]*from[i][1],0.5);
		}
	}

void arr_add_r2d(fftw_complex *from, double *to, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		to[i]+=from[i][0];
		}
	}



double arr_mean_r(fftw_complex *arr, int n)
	{
	double mean=0.;
	int i;
	for (i=0; i<n; i++)
		{
		mean += arr[i][0];
		}
	mean /= n;
	return mean;
	}

double arr_rms_r(fftw_complex *arr, int n, double mean)
	{
	double ssqr=0.;
	int i;
	for (i=0; i<n; i++)
		{
		ssqr += arr[i][0]*arr[i][0];
		}
	ssqr -= n*mean*mean;
	ssqr/= n-1;
	return pow(ssqr,0.5);
	}

void arr_print_r(fftw_complex *arr, int n)
	{
	int i;
	for(i=0; i<n; i++)
		{
		cout<<arr[i][0]<<" ";
		}
	}

void arr_sub_d2r(double *sub, fftw_complex *from, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		from[i][0] -= sub[i];
		}
	}


void arr_write_r(fftw_complex *arr, int n)
	{
	int i;
	for(i=0; i<n ;i++)
		{
		cout<<arr[i][0]<<" ";
		}
	cout<<"\n";
	}

void arr_writel_rc(fftw_complex *arr, int n, char *filename)
	{
	int i;
	ofstream out(filename);
	for(i=0; i<n ;i++)
		{
		out<<i<<" "<<arr[i][0]<<" "<<arr[i][1]<<"\n";
		}
	out.close();
	}

void arr_writel_m(fftw_complex *arr, int n, char *filename, int old)
	{
	int i;
	ofstream out;
	if (old == 0)
		{
		out.open(filename);
		}
	else
		{
		out.open(filename, ios_base::app);
		out<<"\n\n\n";
		}
	for(i=0; i<n;i++)
		{
		out<<i<<" "<<pow(arr[i][0]*arr[i][0]+arr[i][1]*arr[i][1],0.5)<<"\n";
		}
	out.close();
	}

void arr_writel2d_m(fftw_complex *arr, int n1, int n2,  char *filename, int old)
	{
	int i,j;
	ofstream out;
	if (old == 0)
		{
		out.open(filename);
		}
	else
		{
		out.open(filename, ios_base::app);
		out<<"\n\n\n";
		}
	for(i=0; i<n1;i++)
		{
		for (j=0; j<n2; j++)
			{
			out<<i<<" "<<j<<" "<<pow(arr[i*n2+j][0]*arr[i*n2+j][0]+arr[i*n2+j][1]*arr[i*n2+j][1],0.5)<<"\n";
			}
		}
	out.close();
	}

void arr_writel2d_m_grid(fftw_complex *arr, int n1, int n2,  char *filename, int old)
	{
	int i,j;
	ofstream out;
	if (old == 0)
		{
		out.open(filename);
		}
	else
		{
		out.open(filename, ios_base::app);
		out<<"\n\n\n";
		}
	for(i=0; i<n1;i++)
		{
		for (j=0; j<n2; j++)
			{
			out<<" "<<pow(arr[i*n2+j][0]*arr[i*n2+j][0]+arr[i*n2+j][1]*arr[i*n2+j][1],0.5);
			}
		out<<"\n";
		}
	out.close();
	}

void arr_writel_rc(fftw_complex *arr, int n, char *filename, int old)
	{
	int i;
	ofstream out;
	if (old == 0)
		{
		out.open(filename);
		}
	else
		{
		out.open(filename, ios_base::app);
		out<<"\n\n\n";
		}
	for(i=0; i<n ;i++)
		{
		out<<i<<" "<<arr[i][0]<<" "<<arr[i][1]<<"\n";
		}
	out.close();
	}

// here we multiply by taking arr two starred
void arr_cmultiply(fftw_complex *one, fftw_complex *two, fftw_complex *three, int n)
	{
	int i;
	for (i=0; i<n ; i++)
		{
		three[i][0]=one[i][0]*two[i][0]+one[i][1]*two[i][1];
		three[i][1]=-one[i][0]*two[i][1]+one[i][1]*two[i][0];
		}
	}
