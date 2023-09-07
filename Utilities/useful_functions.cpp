const int maxlength=1024, timeout=3;


//double quad (double a, double b, double c, int want);
void myitoa(int n, char*string, int max);
void arraymult (double array[], int dim1, int dim2, double mult);
double dot(double vector1[], double  vector2[], int dim);
double dot(double vector1[3], double vector2[3]);
void printarray(double array[], int n1, int n2);
void norm(double vector[],int n);
double sqr(double n);
double modvector(double vector[], int dim);
void equal(double array1[], double array2[], int dim1, int dim2);
void cross(double vector1[3], double vector2[3], double output[3]);
int intpow(int toraise, int power);
//double abs(double value) {if (value<0) {return -value;} else {return value;}}
double total_angle(double theta, double phi);
long int first_element_lt(long int n_elements,double value, double *table);
int whichpow(long int value);
inline double quadminus (double a, double b, double c);
inline double quadplus (double a, double b, double c);
double quadminus_check (double a, double b, double c);
double quadplus_check (double a, double b, double c);
inline double mod(double v[3]);
void norm(double vector[3]);

double my_atan(double x, double y);

//atan(y/x)
double my_atan(double x, double y)
	{
	double result;
	if (x == 0)
		{
		if (y > 0) {result = piontwo;} else {result=-piontwo;}
		}
	else
		{
		result=atan(y/x);
		if (x < 0)
			{
			if (result < 0) {result += M_PI;}
			else {result -= M_PI;}
			}
		}
	if (isnan(result)) {cout<<"isnan!!! "<<x<<" "<<y<<"\n";}
	return result;
	}

// this set created for telescope.cpp, but could have other uses
void gen_weights(double *frequencies, double *fweights, double band_min, double band_max, int n_freqs);
double gaussian(double mean, double x, double std_dev);
double gaussian(double dist, double std_dev);

inline double cosine_rule(double a, double b, double c);
inline double cosine_angle(double a, double b, double c);
inline double cosine_dist(double a, double b, double theta);
void weightings(double min,double max, double ref1, double ref2,double weights[2]);
double inv_gauss(double result, double std_dev);

// DATA INPUT  
int getnumber(int N0, int Nmax, char abort);
int getyn (char **pmessages);
long int getlong(char **pmessages);
int getint(char **pmessages);
double getdouble(char **pmessages, double max);
void readline(double *vector, int ninputs);
double interpolate(double x, double x1, double y1, double x2, double y2);
double interpolate(double x, double x1, double y1, double x2, double y2)
	{
	double k,y;
	k=(x2-x)/(x2-x1);
	y=k*y1+(1.-k)*y2;
	return y;
	}

// ARRAY MANIPULATION ROUTINES

double modvector(double vector[], int dim) {
	int c1;
	double total=0;
	for (c1=0 ; c1<dim ; c1++) {
		total += sqr(vector[c1]);}
	total=sqrt(total);
	return total;}
	
// code here is r/i/m/p for real/imaginary/mod/phase, which imply a complex type, and d, which implies double
//void arr_add_r2d(fftw_complex *from, double *to, int n)
//	{


long int first_element_lt(long int n_elements,double value, double *table)
	{
	int i;
	long int posn, final=0;
	i=whichpow(n_elements);
	while (n_elements >1)
		{
		i=0;
		
		posn=intpow(2,i)-1;
		if (table[posn+final]< value)
			{
			final += intpow(2,i);
			n_elements -= intpow(2,i);
			i=whichpow(n_elements);
			}
		else
			{
			n_elements=intpow(2,i);
			i-=1;
			}
		}
	return final;
	}

int whichpow(long int value)
	{
	int i;
	while (intpow(2,i+1)<value) {i++;}
	return i;
	}

void myitoa(int n, char*string, int max)
	{
	int count;
	for (count=0; count<max; count++)
		{
		string[max-count-1]=48+(n%(intpow(10L,count+1))-n%(intpow(10L,count)))/intpow(10L,count);
		}
	string[max]='\0';
	}
	/*
double quad (double a, double b, double c, int want)
	{
	float part1,part2;
	part1 = -b/(2.*a);
	part2 = sqr(b)-4.0*a*c;
	if (part2 < 0)
		{
		cout<< "\n sqrt -ve! \n";
		part2=0;
		}
	else
		{
		part2 = sqrt(part2)/(2.0*a);
		}
	switch (want)
		{
		case 0: return part1-part2; break;
		case 1: return part1+part2; break;
		}
	}
*/
inline double quadplus (double a, double b, double c)
	{
	return (-b+sqrt(b*b-4.*a*c))/(2.*a);
	}

double quadplus_check (double a, double b, double c)
	{
	double temp=b*b-4.*a*c;
	if (temp <= 0.) {return 0.;} else {return (-b+sqrt(temp))/(2.*a);}
	}

double quadminus_check (double a, double b, double c)
	{
	double temp=b*b-4.*a*c;
	if (temp <= 0.) {return 0.;} else {return (-b-sqrt(temp))/(2.*a);}
	}

inline double quadminus (double a, double b, double c)
	{
	return (-b-sqrt(b*b-4.*a*c))/(2.*a);
	}

void equal(double array1[], double array2[], int dim1, int dim2)
	{
	int i, j;
	for (i=0 ; i<dim1 ; i++)
		{
		if (dim2 > 0)
			{
			for (j=0 ; j<dim2 ; j++)
				{
				array1[i,j]=array2[i,j];
				}
			}
	
		else
			{
			array1[i]=array2[i];
			}
		}
	}

void equal(double *vecto, double *vecfrom)
	{
	vecto[0]=vecfrom[0];
	vecto[1]=vecfrom[1];
	vecto[2]=vecfrom[2];
	}

inline double sqr(double n)
{	
	return n*n;
}
void norm(double vector[],int n)
{	double l;
	l=dot(vector, vector,n);
	l=sqrt(l);
	arraymult(vector, n, 0, 1./l);
}

void norm(double vector[3])
	{
	double l=mod(vector);
	vector[0] /= l;
	vector[1] /= l;
	vector[2] /= l;
	}
inline double mod(double v[3])
	{
	return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	}

void printarray(double array[], int n1, int n2) {
	int c5, c6;
	if (n2==0) {
		for (c5=0 ; c5<n1 ; c5++){
			cout<<array[c5]<<" , "; }}
	else {
		for (c5=0 ; c5<n1 ; c5++) {
			for (c6=0; c6<n2 ; c6++) {
				cout<<array[c5,c6]; }
			cout<<"\n"; }}
	cout<<"\n"; }

double dot(double vector1[], double  vector2[], int dim)
{	int c7;
	double sum;
	sum=0;
	for (c7=0; c7<dim ; c7++) {
		sum += vector1[c7]*vector2[c7];}
	return sum;}

inline double dot(double vector1[3], double vector2[3])
	{
	return vector1[0]*vector2[0] + vector1[1]*vector2[1] + vector1[2]*vector2[2];
	}

void arraymult(double array[], int dim1, int dim2, double mult)
{	int c8, c9;
	if (dim2==0) {
		for (c8=0 ; c8<dim1 ; c8++) {
			array[c8]=array[c8]*mult; }}
	else {
		for (c8=0 ; c8<dim1 ; c8++) {
			for (c9=0 ; c9<dim2 ; c9++) {
				array[c8,c9]=array[c8,9]*mult; }}}}

int intpow(int toraise, int power)
	{
	int i, result=1;
	for (i=0; i<power; i++)
		{
		result*=toraise;
		}
	return result;
	}

void cross(double vector1[3], double vector2[3], double output[3]){
	output[2]=vector1[0]*vector2[1]-vector1[1]*vector2[0];
	output[0]=vector1[1]*vector2[2]-vector1[2]*vector2[1];
	output[1]=vector1[2]*vector2[0]-vector1[0]*vector2[2];}


// GETS AN UNSIGNED LONG INTERGER RESPONSE   
int getint (char **pmessages)
	{
	int answer, count=0;
	char characters[maxlength];
	while (count<timeout && count>=0)
		{
		cin.ignore(cin.rdbuf()->in_avail());
		cin.getline(characters, maxlength-1,'\n');
		if (atoi(characters) != 0 || strcmp(characters,"0")==0)
			{
			answer = atoi(characters);
			count -= timeout;
			}
		else
			{
			cout<<*(pmessages+15)<<endl;
			count++;
			}
		}
	if (count >= timeout) {return -666;} else {return answer;}
}



// ASKS FOR A YES/NO ANSWER   
int getyn (char **pmessages)
	{
	int answer=2;
	char character;
	while ((answer >1) && (answer < timeout))
		{
		character=cin.get();
		cin.ignore(cin.rdbuf()->in_avail());
		switch (character)
			{
			case 'n':
				answer=0;
				break;
			case 'y':
				answer=1;
				break;
			case 'a':
				answer=-1;
				break;
			default:
				cout<<*(pmessages+22)<<endl;
				answer++;
				break;
			}
		}
	return answer;
	}

// GETS AN INTEGER IN A SPECIFIED RANGE   
int getnumber(int N0, int Nmax, char abort)
	{
	int answer=N0-1, count=0;
	char character[maxlength];
	while (( answer < N0 || answer > Nmax) && count<timeout)
		{
		answer=N0-1;
		if (count!=0)
			{
			cout<<"Please enter a digit between "<<N0<<" and "<<Nmax<<": , or q to exit"<<endl;
			}
		cin.getline(character,maxlength-1);
		cin.ignore(cin.rdbuf()->in_avail());
		if (character[0]==abort && character[1]=='\0') {break;}

		if (atoi(character) != 0 || (N0 <= 0 && Nmax >= 0 && strcmp(character,"0")==0))
			{
			answer = atoi(character);
			}
		count++;
		}
	return answer;
	}
	
// GETS A DOUBLE-PRECISION VALUE   
double getdouble(char **pmessages, double max)
	{
	double answer=-1;
	int count=0, ok=0;
	char characters[maxlength];
	while (ok == 0 && count<timeout)
		{
		cin.ignore(cin.rdbuf()->in_avail());
		cin.getline(characters, maxlength-1,'\n');
		if (atof(characters)>max && max>0)
			{
			cout<<*(pmessages+19)<<max<<endl;
			}
		else
			{
			if (atof(characters) != 0 || characters[0] == '0')
				{
				answer=atof(characters);
				ok=1;
				}
				else
				{
				cout<<*(pmessages+15)<<endl;
				}
			}
		count++;
		}
	return answer;
}


// GETS AN UNSIGNED LONG INTERGER RESPONSE   
long int getlong (char **pmessages)
	{
	long int answer=0;
	int count=0;
	char characters[maxlength];
	while (answer==0 && count<timeout)
		{
		cin.ignore(cin.rdbuf()->in_avail());
		cin.getline(characters, maxlength-1,'\n');
		if (atoi(characters)>0)
			{
			answer=atoi(characters);
			}
		else
			{
			cout<<*(pmessages+15)<<endl;
			}
		count++;
		}
	return answer;
}

double total_angle(double theta, double phi)
	{
	double total;
	total=acos(cos(theta)*cos(phi));
	return total;
	}



// ###########  CREATED FOR TELESCOPE.CPP, BUT COULD BE USEFUL ##############


/*  This function is intended to return coefficients of the electric field squared in each frequency band. The intention is to first calculate a linear interpolation by which the fluxes over each band are functions of the frequencies specified, then simply work out the flux of signal in each band. This can then be processed according to some 'is it detected' function, or simply plotted vs threshold.
	
	The routine must be passed a single frequency band, and also a set of frequencies. It outputs an array which holds the weightings.
*/
void gen_weights(double *frequencies, double *fweights, double band_min, double band_max, int n_freqs)
	{
	double minf,maxf,ref1,ref2, weights[2];
	int i,j,k, up_i, low_i;

/* This returns where the band mins and max's; j gives the highest f the band is less than, i the lowest it is more than */	
	i=0;
	j=0;
	while (i<n_freqs && frequencies[i]<band_min) {i++;}
	while (j<n_freqs && frequencies[j]<band_max) {j++;}
	

// simple weightings for ranges completely covered by the frequency band.	
	for (k=i; k<(j-1); k++)
		{
		minf=frequencies[k];
		maxf=frequencies[k+1];
		ref1=frequencies[k];
		ref2=frequencies[k+1];
		weightings(minf,maxf,ref1,ref2,weights);
		fweights[k] += weights[0];
		fweights[k+1] += weights[1];
		}
//	edge factors - calculates the parts where i is outside the edge of the band.	
	if (i<j)
		{
		/* the lower part */		
		minf=band_min;
		maxf=frequencies[i];
		if (i==0) {up_i=1; low_i=0;} else {up_i=i; low_i=i-1;}
		weightings(minf,maxf,frequencies[low_i],frequencies[up_i],weights);
		fweights[low_i] += weights[0];
		fweights[up_i] += weights[1];
		
		/* the upper part */
		minf=frequencies[j-1];
		maxf=band_max;
		if (j==n_freqs) {up_i=n_freqs-1; low_i=n_freqs-2;} else {up_i=j; low_i=j-1;}
		weightings(minf,maxf,frequencies[low_i],frequencies[up_i],weights);
		fweights[low_i] += weights[0];
		fweights[up_i] += weights[1];
		}
	if (i==j)
		{
		minf=band_min;
		maxf=band_max;
		if (i==0) {up_i=1; low_i=0;}
		if (i==n_freqs) {up_i=i-1; low_i=i-2;}
		if (i != 0 && i != n_freqs) {up_i=i; low_i=i-1;}
		weightings (minf, maxf, frequencies[low_i], frequencies[up_i],weights);
		fweights[low_i] += weights[0];
		fweights[up_i] += weights[1];
		}
	}

inline double inv_gauss(double result, double std_dev)
	{
	return std_dev*sqrt(-log(result)*2);
	}


// standard cosine rule solver for angle - c^2 = a^2 + b^2 - 2abcos(theta)
inline double cosine_angle(double a, double b, double c)
	{
	return acos((a*a+b*b-c*c)/(2*a*b));
	}

inline double cosine_rule(double a, double b, double c)
	{
	return acos((a*a+b*b-c*c)/(2*a*b));
	}

// standard cosine rule which solves for distance
inline double cosine_dist(double a,double b, double theta)
	{
	return sqrt(a*a+b*b-2*a*b*cos(theta));
	}

// if all values estimated from ref1 and ref2 via linear interpolation, then int f(x) from min to max = weight[0]*ref1 + weight[1]*ref2
void weightings (double min,double max, double ref1, double ref2,double weights[2])
	{
	double r;
	r=((min+max)/2. - ref1)/(ref2-ref1);
	weights[0]=(1-r)*(max-min);
	weights[1]=r*(max-min);
	}

// a non-normalised gaussian - the peak is at one
inline double gaussian(double mean, double x, double std_dev)
	{
	return exp(-sqr((x-mean)/std_dev)/2.);
	}

inline double gaussian(double offset, double std_dev)
	{
	return exp(-sqr(offset/std_dev)/2.);
	}

void readline(double *vector, int ninputs)
	{
	int i, max=256;
	char tempvector[max];
	for (i=0; i<ninputs-1; i++)
		{
		cin.getline(tempvector,max-1,' ');
		vector[i]=atof(tempvector);
		}
	cin.getline(tempvector, max-1,'\n');
	vector[ninputs-1]=atof(tempvector);
	}

