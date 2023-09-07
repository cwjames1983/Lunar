void arr_write2l(double *arr1, double *arr2, int n, char *filename, int old);
void arr_add(double *to, int toadd, int n);
void arr_cout(int *arr, int n);
void arr_cout(int *arr, int n1, int n2);
void arr_zero(double *arr, int n);
void arr_zero(int *arr, int n);
void arr_add(double *from, double *to, int n);
void arr_add(int *from, int *to, int n);
void arr_add(int *from, double *to, int n);
void arr_add(double *to, double toadd, int n);
void arr_equal(double *from, double *to, int n);
void arr_safe_div(int *div, double *arr, int n);
void arr_safe_sqrt_div(int *div, double *arr, int n);

void arr_div(double *arr, double d, int n);
void arr_div(int *div, double *arr, int n);

void arr_sub(double *sub, double *from, int n);
double arr_rms(double *arr, int n, double mean);
double arr_mean(double *arr, int n);
void arr_sub(double *sub,double *from, int n);
void arr_add_r2d(double *from, double *to, int n);
void arr_write(double *arr, ofstream &out, int n, double x1, double x2);
void arr_write(double *arr, int n);
void arr_writel(double *arr, int n, char *filename);
void arr_pow(double *arr, double pwr, int n);
void arr_mpow(double *arr, double pwr, int n);
void arr_div(double *div, double *arr, int n);
double arr_sum(double *arr, int n1, int n2);
double arr_copy(double *from, double *to, int n);
void arr_mult(double mult, double *arr, int n);
void arr_writel(double *arr, int n, char *filename, int old);
void arr_write2d(double *arr, int n1, int n2, char *filename, int old);
void arr_write2d_grid(double *arr1, int n1, int n2, char *filename, int old);
void arr_writel_rev(double *arr, int n, char *filename, int old);
void arr_mult(double *arr, double *by, int n);

void arr_writel(double *arr1, double *arr2, int n, char *filename, int old);
void arr_writel_rev(double *arr1, double *arr2, int n, char *filename, int old);

void arr_lin_eq(double *arr1, double *arr2, double *to, int n, double x1, double x2);
void arr_lin_eq(double *arr1, double *arr2, double *to, int n, double x1, double x2)
	{
	int i;
	for (i=0; i<n; i++)
		{
		to[i]=x1*arr1[i]+x2*arr2[i];
		}
	}

void arr_write(double *arr, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		cout<<arr[i]<<" ";
		}
	cout<<"\n";
	}

void arr_writel(double *arr, int n, char *filename)
	{
	int i;
	ofstream out;
	out.open(filename);
	for (i=0; i<n; i++)
		{
		out<<arr[i]<<"\n";
		}
	out.close();
	}

void arr_mult(double *arr, double *mult, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		arr[i] *= mult[i];
		}
	}

void arr_mult(double mult, double *arr, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		arr[i] *= mult;
		}
	}

double arr_copy(double *from, double *to, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		to[i]=from[i];
		}
	}

double arr_sum(double *arr, int n1, int n2)
	{
	int i;
	double sum=0.;
	for(i=n1; i<=n2; i++)
		{
		sum += arr[i];
		}
	return sum;
	}


void arr_div(double *div, double *arr, int n)
	{
	int i;
	for(i=0; i<n; i++)
		{
		arr[i] /= div[i];
		}
	}

void arr_div(int *div, double *arr, int n)
	{
	int i;
	for(i=0; i<n; i++)
		{
		arr[i] /= div[i];
		}
	}

void arr_safe_div(int *div, double *arr, int n)
	{
	int i;
	for(i=0; i<n; i++)
		{
		if (div[i] != 0)
			{
			arr[i] /= div[i];
			}
		else
			{
			if (arr[i] != 0)
				{
				arr[i]=-1e30;
				}
			}
		}
	}
	
	
void arr_safe_sqrt_div(int *div, double *arr, int n)
	{
	int i;
	for(i=0; i<n; i++)
		{
		if (div[i] != 0)
			{
			arr[i] /= sqrt((double) div[i]);
			}
		else
			{
		/*	if (arr[i] != 0)
				{
				arr[i]=-1e30;
				} */
			}
		}
	}

void arr_pow(double *arr, double pwr, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		arr[i]=pow(arr[i],pwr);
		}
	}
	
void arr_mpow(double *arr, double pwr, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		arr[i]=pow(fabs(arr[i]),pwr);
		}
	}

void arr_write(double *arr, ofstream &out, int n, double x1, double x2)
	{
	int i;
	double x=x1, dx=(x2-x1)/(n-1);
	for (i=0; i<n; i++)
		{
		out<<x<<" "<<arr[i]<<"\n";
		x += dx;
		}
	}
	


void arr_add_r2d(double *from, double *to, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		to[i]+=from[i];
		}
	}

void arr_sub(double *sub, double *from, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		from[i] -= sub[i];
		}
	}

double arr_mean(double *arr, int n)
	{
	double mean=0.;
	int i;
	for (i=0; i<n; i++)
		{
		mean += arr[i];
		}
	mean /= n;
	return mean;
	}

double arr_rms(double *arr, int n, double mean)
	{
	double ssqr=0.;
	int i;
	for (i=0; i<n; i++)
		{
		ssqr += arr[i]*arr[i];
		}
	ssqr -= n*mean*mean;
	ssqr/= n-1;
	return pow(ssqr,0.5);
	}

void arr_sub_d2d(double *sub, double *from, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		from[i] -= sub[i];
		}
	}

void arr_div(double *arr, double d, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		arr[i]/=d;
		}
	}

void arr_zero(double *A, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		A[i]=0;
		}
	}

void arr_zero(int *A, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		A[i]=0;
		}
	}
	
void arr_add(double *to, double toadd, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		to[i]+=toadd;
		}
	}
	
void arr_add(double *to, int toadd, int n)
	{
	int i;
	double tadd=(double) toadd;
	for (i=0; i<n; i++)
		{
		to[i]+=tadd;
		}
	}

void arr_add(double *from, double *to, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		to[i]+=from[i];
		}
	}

void arr_add(int *from, double *to, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		to[i]+=from[i];
		}
	}

void arr_add(int *from, int *to, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		to[i]+=from[i];
		}
	}
void arr_equal(double *from, double *to, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		to[i]=from[i];
		}
	}

void arr_writel(double *arr, int n, char *filename, int old)
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
		out<<i<<" "<<arr[i]<<"\n";
		}
	out.close();
	}

void arr_writel_rev(double *arr, int n, char *filename, int old)
	{
	int i, half=n/2;
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
	for(i=half; i<n ;i++)
		{
		out<<i-half<<" "<<arr[i]<<"\n";
		}
	for(i=0; i<half ;i++)
		{
		out<<i+half<<" "<<arr[i]<<"\n";
		}
	out.close();
	}

void arr_write2l(double *arr1, double *arr2, int n, char *filename, int old)
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
		out<<arr1[i]<<" "<<arr2[i]<<"\n";
		}
	out.close();
	}

void arr_writel(double *arr1, double *arr2, int n, char *filename, int old)
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
		out<<i<<" "<<arr1[i]+arr2[i]<<"\n";
		}
	out.close();
	}

void arr_write2d(double *arr1, int n1, int n2, char *filename, int old)
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
		out<<"\n\n\n\n";
		}
	for(i=0; i<n1 ;i++)
		{
		for (j=0; j<n2; j++)
			{
			out<<i<<" "<<j<<" "<<arr1[i*n2+j]<<"\n";
			}
		}
	out.close();
	}

void arr_write2d_grid(double *arr1, int n1, int n2, char *filename, int old)
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
		out<<"\n\n\n\n";
		}
	for(i=0; i<n1 ;i++)
		{
		for (j=0; j<n2; j++)
			{
			out<<arr1[i*n2+j]<<" ";
			}
		out<<"\n";
		}
	out.close();
	}

void arr_writel_rev(double *arr1, double *arr2, int n, char *filename, int old)
	{
	int i, half=n/2;
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
	for(i=half; i<n ;i++)
		{
		out<<i-half<<" "<<arr1[i]+arr2[i]<<"\n";
		}
	for(i=0; i<half ;i++)
		{
		out<<i+half<<" "<<arr1[i]+arr2[i]<<"\n";
		}
	out.close();
	}

void arr_cout(int *arr, int n)
	{
	int i;
	for (i=0; i<n; i++)
		{
		cout<<arr[i]<<" ";
		}
	cout<<"\n";
	}

void arr_cout(int *arr, int n1, int n2)
	{
	int i,j;
	for (i=0; i<n1; i++)
		{
		for (j=0; j<n2; j++)
			{
			cout<<arr[i*n2+j]<<" ";
			}
		cout<<"\n";
		}
	cout<<"\n\n";
	}


