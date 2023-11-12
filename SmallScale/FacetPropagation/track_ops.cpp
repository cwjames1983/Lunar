void gensubtrack(track &atrack, long double ntdiv, int nthdiv, track &subtrack)
	{
	long double k1,k2;
	int i;
	k1=nthdiv/ntdiv;
	k2=(nthdiv+1)/ntdiv;
	
	subtrack.x0=k1*atrack.x1+(1-k1)*atrack.x0;
	subtrack.x1=k2*atrack.x1+(1-k2)*atrack.x0;
	subtrack.centre[0]=(subtrack.x0+subtrack.x1)/2.;
	subtrack.y0=k1*atrack.y1+(1-k1)*atrack.y0;
	subtrack.y1=k2*atrack.y1+(1-k2)*atrack.y0;
	subtrack.centre[1]=(subtrack.y0+subtrack.y1)/2.;
	subtrack.z0=k1*atrack.z1+(1-k1)*atrack.z0;
	subtrack.z1=k2*atrack.z1+(1-k2)*atrack.z0;
	subtrack.centre[2]=(subtrack.z0+subtrack.z1)/2.;
	
	subtrack.t0=k1*atrack.t1+(1-k1)*atrack.t0;
	subtrack.t1=k2*atrack.t1+(1-k2)*atrack.t0;
	subtrack.centre[3]=(subtrack.t0+subtrack.t1)/2.;
	//ts[0]=(ts[1]+ts[2])/2.;
	
	subtrack.delt=atrack.delt/ntdiv;
	
	for (i=0; i<3; i++) {subtrack.vec[i]=atrack.vec[i];}
	subtrack.length=atrack.length/ntdiv;
	subtrack.hl=subtrack.length/2.;
	
	subtrack.p[0][0]=subtrack.x0;
	subtrack.p[0][1]=subtrack.y0;
	subtrack.p[0][2]=subtrack.z0;
	subtrack.p[1][0]=subtrack.x1;
	subtrack.p[1][1]=subtrack.y1;
	subtrack.p[1][2]=subtrack.z1;
	
	subtrack.q=atrack.q;
	
	subtrack.beta=atrack.beta;
	}



// assumes the tracks are such that z is the shower axis, and x/y are perpendicular to this. We want a 'standard' 0-transform though to correspond to the cascade being parallel to the surface, i.e. the shower axis in the xy plane. So I swap the cascade x and z. // theta positive is into the surface, -ve is upwards. phi +ve is towards +ve y direction; 0 is x direction
void transform_tracks(long double theta, long double phi, long double dx, long double dy, long double dz, int ntracks, track *tracks)
	{
	long double costheta, sintheta, cosphi, sinphi;
	long double x,y,z;
	long double xprime, yprime, zprime;
	int i;
	
	costheta=cos(theta);
	sintheta=sin(theta);
	cosphi=cos(phi);
	sinphi=sin(phi);
	
	for (i=0; i<ntracks; i++)
		{
		// swaps x and z, so now the cascade is travelling in the x-direction.
		x=tracks[i].z0;
		y=tracks[i].y0;
		z=tracks[i].x0;
		
		//offsets the cascade
		
		// rotation about y (in x-z plane) -- lowers the cascade
		xprime=costheta*x + sintheta*z;
		zprime=costheta*z - sintheta*x;
		x=xprime;
		z=zprime;
		
		// rotation about z (in x-y plane) -- rotates it. phi=0 leaves everything unchanged, +ve phi turn x=1,y=0 into positive y-direction
		xprime=(cosphi*x-sinphi*y);
		yprime=(sinphi*x+cosphi*y);
		x=xprime;
		z=zprime;
		
		x += dx+X2;
		y += dy+Y2;
		z += dz;
		
		cout<<"Track "<<i<<" p0 = ("<<x<<","<<y<<","<<z<<")   ";
		
		// reassigns the shower
		tracks[i].x0=x;
		tracks[i].y0=y;
		tracks[i].z0=z;
		
		// swaps x and z, so now the cascade is travelling in the x-direction.
		x=tracks[i].z1;
		y=tracks[i].y1;
		z=tracks[i].x1;
		
		//offsets the cascade
		
		// rotation about y (in x-z plane) -- lowers the cascade
		xprime=costheta*x + sintheta*z;
		zprime=costheta*z - sintheta*x;
		x=xprime;
		z=zprime;
		
		// rotation about z (in x-y plane) -- rotates it
		xprime=(cosphi*x-sinphi*y);
		yprime=(sinphi*x+cosphi*y);
		x=xprime;
		z=zprime;
		
		x += dx+X2;
		y += dy+Y2;
		z += dz;
		
		cout<<"  p1 = ("<<x<<","<<y<<","<<z<<")\n";
		
		// reassigns the shower
		tracks[i].x1=x;
		tracks[i].y1=y;
		tracks[i].z1=z;
		}
	}


void print_tracks(char *filename, track *tracks, int ntracks, int dowhat)
	{
	ofstream out;
	int i;
	if (dowhat == -1)
		{
		for (i=0; i<ntracks; i++)
			{
			cout<<"Track "<<i<<": "<<tracks[i].x0<<" "<<tracks[i].y0<<" "<<tracks[i].z0<<" "<<tracks[i].x1<<" "<<tracks[i].y1<<" "<<tracks[i].z1<<"\n";
			}
		}
	else
		{
		if (dowhat==0)
			{
			out.open(filename);
			out<<"# Track   X0     Y0     Z0    X1    Y1     Z1\n";
			}
		else {out.open(filename, ios::app);}
		for (i=0; i<ntracks; i++)
			{
			out<<i<<": "<<tracks[i].x0<<" "<<tracks[i].y0<<" "<<tracks[i].z0<<" "<<tracks[i].x1<<" "<<tracks[i].y1<<" "<<tracks[i].z1<<"\n";
			}
		}
	return;
	}
	
	
// this routines just gets tracks from the tracks file
// the format for x,y,z is just in m. The assumption is that the main cascade
// started at 0,0 and moves in the z-direction
// the two times are: start in seconds, end in seconds
int read_tracks(char *track_file, track *&tracks, int &ntracks)
	{
	ifstream in;
	int i;
	
	in.open(track_file);
	if (!in.is_open())
		{
		return 0;
		}
	in>>ntracks;
	tracks = new track[ntracks];
	for (i=0; i<ntracks; i++)
		{
		in>>tracks[i].t0>> tracks[i].x0>> tracks[i].y0>> tracks[i].z0>> tracks[i].beta>>tracks[i].x1>> tracks[i].y1>> tracks[i].z1>>tracks[i].q;
		}
	return 1;
	}
	
void init_track(track &atrack)
	{
	atrack.p[0][0]=atrack.x0;
	atrack.p[0][1]=atrack.y0;
	atrack.p[0][2]=atrack.z0;
	atrack.p[1][0]=atrack.x1;
	atrack.p[1][1]=atrack.y1;
	atrack.p[1][2]=atrack.z1;
	atrack.length=calc_dist(atrack.p[0], atrack.p[1]);
	atrack.t1=atrack.t0+atrack.length/(atrack.beta*c_light);
	atrack.vec[0]=atrack.x1-atrack.x0;
	atrack.vec[1]=atrack.y1-atrack.y0;
	atrack.vec[2]=atrack.z1-atrack.z0;
	atrack.centre[0]=(atrack.x1+atrack.x0)/2.;
	atrack.centre[1]=(atrack.y1+atrack.y0)/2.;
	atrack.centre[2]=(atrack.z1+atrack.z0)/2.;
	atrack.centre[3]=(atrack.t1+atrack.t0)/2.;
	
	atrack.hl=atrack.length/2.;
	atrack.delt=atrack.t1-atrack.t0;
	normalise(atrack.vec);
//	atrack.beta = min(1.,atrack.length/atrack.delt/c_light);
	TRACKS_C[0] += atrack.centre[0]*atrack.q*atrack.length;
	TRACKS_C[1] += atrack.centre[1]*atrack.q*atrack.length;
	TRACKS_C[2] += atrack.centre[2]*atrack.q*atrack.length;
	TRACKS_SUM += atrack.q*atrack.length;
	}

