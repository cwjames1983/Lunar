#include <fstream>
#include <iostream>

using namespace std;

/*
Code to generate a trivial (completely flat) surface file.
Compile with g++ make_trivial_surface.cpp -o mts.exe
The grid size and spacing are hard-coded.
The Hurst exponent, H, is meaningles, since Srms is zero.
*/

int main()
	{
	int nx=1024, ny=1024;
	long double ds=0.2, Srms=0., H=0.78;
	char *opfile="trivial.dat";
	ofstream out(opfile);
	int i,j;
	
	out<<"# "<<nx<<" "<<ny<<" "<<H<<" "<<ds<<" "<<Srms<<"\n";
	for (i=0; i<nx; i++)
		{
		for (j=0; j<ny; j++)
			{
			out<<"0 ";
			}
		out<<"\n";
		}
	out.close();
	}
