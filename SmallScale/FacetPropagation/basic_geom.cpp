

void cross(long double v1[3], long double v2[3], long double v3[3])
	{
	v3[2] = v1[0]*v2[1] - v2[0]*v1[1];
	v3[0] = v1[1]*v2[2]-v2[1]*v1[2];
	v3[1] = v1[2]*v2[0]-v2[2]*v1[0];
	normalise(v3);
	}

long double my_atan_xy(long double x, long double y)
	{
	long double result;
	if (x==0)
		{
		if (y > 0.) {result=piontwo;} else {result = -piontwo;}
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
	return result;
	}


// contains geometric routines

// calculates the distance fro,
// normal must be normalised, D=0 always
// solve A*x+By+Cz+t=0
long double dtoplane(long double *point, long double *normal)
	{
	long double t=0;
	int i;
	for (i=0; i<3; i++) {t -= point[i]*normal[i];}
	return t;
	}

// gets a normalised vector of p2-p1
void construct_norm_vector(long double p1[3], long double p2[3], long double vec[3])
	{
	int i;
	for (i=0; i<3; i++)
		{
		vec[i]=p2[i]-p1[i];
		}
	normalise(vec);
	}
// gets a normalised vector of p2-p1
void construct_vector(long double p1[3], long double p2[3], long double vec[3])
	{
	int i;
	for (i=0; i<3; i++)
		{
		vec[i]=p2[i]-p1[i];
		}
	}
void normalise(long double vec[3])
	{
	long double norm=0.;
	int i;
	for (i=0; i<3; i++)
		{
		norm += vec[i]*vec[i];
		}
	norm=pow(norm,0.5l);
	for (i=0; i<3; i++)
		{
		vec[i]/=norm;
		}
	}


long double calc_dist(long double x1, long double y1, long double z1, long double x2, long double y2, long double z2)
	{
	return pow(pow(x1-x2,2)+pow(y1-y2,2) + pow(z1-z2,2),0.5l);
	}

long double calc_dist(long double p1[3], long double p2[3])
	{
	long double dist=0;
	int i;
	for (i=0; i<3; i++)
		{
		dist += pow(p1[i]-p2[i],2);
		}
	dist = pow(dist,0.5l);
	return dist;
	}

// calculates distance between two points, and returns a vector between them (normalised)
long double calc_dist(long double p1[4], long double p2[3], long double vec[3])
	{
	long double dist=0;
	int i;
	for (i=0; i<3; i++)
		{
		vec[i]=p2[i]-p1[i];
		dist += pow(vec[i],2);
		}
	dist = pow(dist,0.5l);
	for (i=0; i<3; i++) {vec[i]/=dist;}
	return dist;
	}


// simple dot product
long double dot(long double *vec1, long double *vec2)
	{
	return vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
	}
	
// distance from a point to a plane with Ax+By+Cz=0 and point a,b,c
long double find_dist_zero_D(long double A, long double B, long double C, long double a, long double b, long double c)
	{
	long double t;
	t=-A*a-B*b-C*c;
	return t;
	}
// not necessarily xs!!!
long double get_height(long double x0, long double x1, long double z0, long double z1, long double x)
	{
	long double k,z;
	k=(x-x0)/(x1-x0);
	z=k*z1+(1.-k)*z0;
	return z;
	}
	
void my_itoa(int number, int n, char *string)
	{
	long double temp;
	int remainders[n];
	int i;
	
	for (i=0; i<n; i++)
		{
		remainders[n-i-1]=number % 10;
		number -= remainders[n-i-1];
		number /= 10;
		string[n-i-1]=48+remainders[n-i-1];
		}
	string[n]='\0';
	return;
	}
