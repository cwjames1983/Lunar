/***** 
Program to calculate the Askaryan radiation from a cascade under a rough surface

The cascade is represented by a number of discrete particle tracks.
The surface is represented by a regular grid of height functions.
The radiation is represented by complex numbers at discrete
frequencies and exit angles.

The program contains many loops; these are:
1: Loop over individual tracks (1D)
  2: Loop over frequency divisions (1D)
    3: Loop over facets (2D)
      4: Loop over subfacets and subtracks (3D)
        5: Loop over output angles (2D)
	  6: Loop over frequencies (remainder of dimension in '2').

That is, we have nine dimensions of nested loops. This program might take a long time!
The aim is to represent a cascade through only a few tracks.


The geometry for this program assumes a centre at 0,0, so the surface is from -X2,-Y2 to X2, Y2 
However, check this! We can use this for the surface, or change the inputs

******/

int FMODE=1;
int NIX, NIY;

#include "rough_header.h"

int main()
	{
	
// generic and input variables
	int i,j,k, OK, ix, iy;
	char surf_file[mfl], track_file[mfl], modsurf_file[mfl], stats_file[mfl], tempstring[mfl];
	char freq_file[mfl], theta_file[mfl], phi_file[mfl], time_file[mfl]; // to contain outputs
	ifstream in;
	ofstream out;
	int timeyn, phiyn, thetayn, freqyn, writemodyn;
	
// variables relating to the tracks
	int ntracks;
	track *tracks;
	long double theta,phi,dx,dy,depth,dz;
	char modtrackfile[mfl];
	
// variables connected with the surface
	float *surface, *rsurface;
	int naveraged, nosx, nosy, surf_file_type, Nconcat;
	long double rX,rY,rdd,rX2,rY2, H, Srms;
	int rnx, rny;
	long double alengthonghz;
	
// relating to outputs
	long double dphi, phi0, dtheta, theta_0, theta_1, f0, f1;
	int hnphi, logyn;
	long double *phis, *freqs, *thetas;
	long double ***normals; // contains normals to the outgoing planes for theta/phi.

// misc temporary variables
	float fmemory; // memory used for facets
	int sqrtn, reportn;
// ######## program starts ##########	
	
	cin>>USE_WHICH>>ERROR_MARGIN>>SIN_APPROX>>MAX_SPLITS>>FMODE>>SIN_APPROX_2>>TDIV_METH>>TDIV_CONST>>FDIV_METH>>FDIV_FACTOR; cin.getline(tempstring,mfl);
	
	cout<<USE_WHICH<<" "<<ERROR_MARGIN<<" "<<SIN_APPROX<<" "<<MAX_SPLITS<<" "<<FMODE<<" "<<SIN_APPROX_2<<" "<<TDIV_METH<<" "<<FDIV_METH<<" "<<FDIV_FACTOR<<"\n";
	
// ######## handles the surface #########
	TIMES[0]=time(0);
	cin>>surf_file_type>>surf_file; cin.getline(tempstring,mfl);
	cout<<"Reading in surface from file "<<surf_file<<endl;//"\n";
	if (surf_file_type == 1)
		{
		cout<<"reading ascii surface..."<<"\n";
		read_surface(surf_file, rsurface, rnx, rny, rdd, H, Srms);
		}
	else
		{
		cout<<"Reading binary file..."<<"\n";
		read_float_surface(surf_file, rsurface, rnx, rny, rdd, H, Srms);
		}
	cout<<"This surface has "<<rnx<<" x "<<rny<<" facets of size "<<rdd<<"^2 m^2.\n";
	
// calculates a new surface based on taking the mean of the generated surface
	cin>>naveraged>>nosx>>nosy>>Nconcat; cin.getline(tempstring,mfl);
	cout<<"Creating a new surface, averaged over "<<naveraged<<" steps of the previous surface, ";
	cout<<"starting at "<<nosx<<", "<<nosy<<", ";
	cout<<" and concatenating the outer "<<Nconcat<<" facets:"<<endl;
	
	mod_surface(rsurface, rnx, rny, naveraged, nosx, nosy, surface, rdd, Nconcat);
	
	delete [] rsurface; rsurface=0; // frees up much-needed memory by unitialising the 'real' surface!
	Srms *= pow(naveraged,H-1.); // modifies the Srms according to the change of scale
	cout<<"The modified surface is "<<NX<<"x"<<NY<<", with total dimensions "<<X<<"x"<<Y<<" m."<<endl;
	
	cin>>writemodyn>>modsurf_file; cin.getline(tempstring,mfl);
	if (writemodyn == 1)
		{
		cout<<"Writing the modified surface to "<<modsurf_file<<endl;
		print_surface(modsurf_file, surface, NX, NY, H, DD, Srms, 2766);
		}
	else
		{
		cout<<"Not writing the modified surface to "<<modsurf_file<<endl;
		}
	TIMES[0]=time(0)-TIMES[0];
	
// ######## handles tracks #########
	TIMES[2]=time(0);
	cin>>track_file; cin.getline(tempstring,mfl);
	cout<<"Reading data on shower tracks from "<<track_file<<endl;
	OK=read_tracks(track_file, tracks, ntracks);
	if (OK ==0)
		{
		cout<<"Error in reading track file: returning..."<<endl;
		return 0;
		}
	cout<<"There are "<<ntracks<<" tracks in the track file."<<endl;
	
//	print_tracks(0,tracks, ntracks,-1);
	
	// gets two angles with which to transform the tracks, and an xyz position for cascade starting position relative to the 0,0,0 of the surface. The only fancy thing is the use of 'depth', so the z-coordinate is relative to the surface.
	cin>>theta>>phi>>dx>>dy>>depth; cin.getline(tempstring,mfl);
	cout<<"The tracks are to be rotated by "<<theta<<" degrees about the y-axis in the x-z plane, and rotated "<<phi<<" degrees in the x-y plane from x towards y.\n";
	
	theta *= torad;
	phi *= torad;
	
	if (fabs(dx) >= X2 || fabs(dy) >= Y2) {cout<<"Specified offsets of the cascade are outside the surface bounds. Resetting to dx=dy=0.\n"; dx=0.; dy=0.;}
	
	depth_to_z(dx,dy,depth,dz,surface);
	cout<<"The zero-point of the cascade is ("<<dx<<","<<dy<<","<<dz<<"), i.e. "<<-depth<<" below the local surface."<<endl;
	
	transform_tracks(theta,phi,dx,dy,dz, ntracks, tracks);
	for (i=0; i<3; i++) {TRACKS_C[i]=0.;}
	TRACKS_SUM=0.;
	NTRACKS=ntracks;
	for (i=0; i<NTRACKS; i++) {init_track(tracks[i]);}
	for (i=0; i<3; i++) {TRACKS_C[i]/= TRACKS_SUM;}
	cout<<"Tracks initialised!"<<endl;
	
//####
//  Note that here, we expect the start point of a cascade to be 0,0,0 and travel in the z-direction. But this gets transformed to be at 0,0,0 and travel in the x-direction, i.e. start 0,0,0, travel 1,0,0.
//####	
	
	strcpy(modtrackfile,"modded_");
	strcat(modtrackfile, track_file);
	print_tracks(modtrackfile,tracks, ntracks,0);
	cout<<"The modded tracks have been written to "<<modtrackfile<<"\n";
	TIMES[2]=time(0)-TIMES[2];
	
// ####### Finished with tracks. Now getting medium properties #######
	
	TIMES[3]=time(0);
	// entering medium properties
	cin>>REFNi>>REFNt>>alengthonghz; cin.getline(tempstring,mfl);
	ALENGTHONHZ=alengthonghz*1e9;
	cout<<"Refractive indices are "<<REFNi<<" and "<<REFNt<<" for the cascade / detection media respectively; the attenuation length at one GHz is "<<alengthonghz<<" m."<<endl;
	
/*	
#####
At this point we have some tracks, and a surface. Now we want to specify at what frequencies, and how many directions etc, we want to calculate the radiation. This is all specified relative to the angle of rotation phi in the x-y plane. So the definition of phi becomes the same, albeit it re-zeroed at the previous phi, while theta is redifined to be zero parallel to the surface, and positive above the surface
####
*/
	// now we talk about theta and phi in terms of angles above the surface:
	cin>>phi0>>hnphi>>dphi; cin.getline(tempstring,mfl); // nphi = 2*hnphi+1 (ensures symmetry)
	cin>>NTHETAS>>theta_0>>theta_1; cin.getline(tempstring,mfl);
	cin>>NFREQS>>f0>>f1>>logyn; cin.getline(tempstring,mfl);// specifies the frequencies we want to calculate stuff at
	
	// determins what sort of outputs are required
	cin>>stats_file; cin.getline(tempstring,mfl);
	cin>>freqyn>>freq_file; cin.getline(tempstring,mfl); // write columns const theta, rows const phi, tables const freq
	cout<<"Writing freq file "<<freq_file<<"? "<<freqyn<<endl;
	cin>>phiyn>>phi_file; cin.getline(tempstring,mfl); // writes columns const freq, rows const theta, tables const phi
	cout<<"Writing phi file "<<phi_file<<"? "<<phiyn<<endl;
	cin>>thetayn>>theta_file; cin.getline(tempstring,mfl);
	cout<<"Writing theta file "<<theta_file<<"? "<<thetayn<<endl;
	cin>>timeyn>>time_file; cin.getline(tempstring,mfl);
	cout<<"Writing time file "<<time_file<<"? "<<timeyn<<endl;
//	in>>LU_nsize>>LU_ntheta>>LU_nphi;
/*	in>>SURF_nsurfs;
	
	SURF_intracks = new int[SURF_nsurfs];
	SURF_infreqs = new int[SURF_nsurfs];
	SURF_surfnames = new char[mfl*SURF_nsurfs];
	for (i=0; i<SURF_nsurfs; i++)
		{
		in>>SURF_intracks[i]>>SURF_infreqs[i]>>SURF_surfnames+i*mfl;
		}
	*/
//	in.close();
	
// ###### Now calculating output angles, frequencies, and the like ######
	
	NPHIS=2*hnphi+1;
	NOPS=NPHIS*NTHETAS;
	dphi *= torad; theta_0 *= torad; theta_1 *= torad; phi0 *= torad;
	
	gen_phi(phi0, hnphi, dphi);
	gen_theta(theta_0,theta_1);
	gen_fs(f0,f1, logyn);
	
	cout<<"Performing initialisation (mostly of shortcuts)..."<<endl;
	fmemory=sizeof(facet)*(NX-1);
	fmemory *= (NY-1)/1024./1024.;
	cout<<"The memory required for facet initialisation is "<<fmemory<<" MB\n";
	if (fmemory > MAX_SYS_MB && FMODE != 2)
		{
		FMODE=2;
		cout<<"Changing FMODE to allow enough memory to be used;"<<endl;
		}
	OK=init(out, surface); // initialises lots of stuff, mostly shortcuts
	cout<<"    ....initialisation completed."<<endl;
	
	if (OK != 1)
		{
		cout<<"Failed to allocate memory. Exiting..."<<endl;
		exit(1);
		}
	
	TIMES[3]=time(0)-TIMES[3];
// ########## FINISHED INITIALISATION ##########

	TIMES[4]=time(0);

// ######## Program starts for real ########	
	
	if (FMODE == 1) // i.e. facets an inner loop
		{
		NIX=1;
		NIY=1;
		}
	else if (FMODE == 2) // facets an outer loop
		{
		NIX=NX-1;
		NIY=NY-1;
		sqrtn=(int) pow((double) NIX,0.5);
		}
	else // does not make sense, and program will not do anything
		{
		NIX=0;
		NIY=0;
		}
	cout<<"Beginning main computation loop: there will be "<<NIX<<" "<<NIY<<" inner x-y loops, each over "<<ntracks<<" tracks, and "<<NFDIVS<<" separate frequency divisions, for a total number of "<<NIX*NIY*ntracks*NFDIVS<<" calls. Now, go get some coffee..."<<endl;
	reportn=pow((double) (NIX*NIY*ntracks*NFDIVS),0.5);
	for (ix=0; ix<NIX; ix++)
		{
		if (FMODE == 2 && ix%sqrtn==0)
			{
			cout<<"Starting row "<<ix<<endl;
			}
		for (iy=0; iy<NIY; iy++)
			{
			if (FMODE == 2)
				{
				init_facet(surface, NX, NY, ix, iy);
				}
			// loops through tracks
			for (i=0; i<ntracks; i++)
				{
				// does this for all facets, unless skipped; begins with facets directly above
				for (j=0; j<NFDIVS; j++)
					{
				//	cout<<"looping facets for track "<<i<<" and frequency division "<<j<<".\n";
					LOOP_COUNT=0;
					WRITESURF=0;
				/*	for (k=0; k<SURF_nsurfs; k++)
					{
				//	cout<<i<<" "<<j<<" "<<SURF_intracks[k]<<" "<<FDIV_STARTS[j]<<" "<<SURF_infreqs[k]<<" "<<FDIV_STARTS[j]+NPERFDIV[j]<<endl;
					if (i==SURF_intracks[k] && FDIV_STARTS[j]<=SURF_infreqs[k] && FDIV_STARTS[j]+NPERFDIV[j] > SURF_infreqs[k]) {WRITESURF=1; SURF_file.open(SURF_surfnames+k*mfl); break;}
					} */
					if (FMODE == 1)
						{
						loop_facets(tracks[i], FDIV_STARTS[j], NPERFDIV[j], surface);
						}
					else if (FMODE==2)
						{
						div_facet_track(tracks[i], FACETS[0][0], FDIV_STARTS[j], NPERFDIV[j], ix, iy, surface);
						}
					if (WRITESURF==1) {SURF_file.close();}
					if (FMODE == 1) {cout<<"   ... looped through "<<LOOP_COUNT<<" loops."<<endl;}
					}
				}
			}
		}
	
	TIMES[4] = time(0)-TIMES[4];
	cout<<"Total time for looping: "<<TIMES[4]<<endl;
	
	cout<<"Maximum divisions...\n...of tracks: "<<MAX_TDIV<<"\n...of facets: "<<MAX_FDIV<<endl;
	cout<<"We had a total of "<<N_VISIBLE<<" visible subfacets, and "<<N_INVISIBLE<<" where the output path was blocked"<<endl;
	cout<<"Facets were split "<<NSPLITCALLS<<" times.\n";
	
// ########## Program has finished - write info ###########

	write_stats(stats_file, rnx, rny);
	write_all_outputs(freqyn, freq_file, thetayn, theta_file, phiyn, phi_file, timeyn, time_file);
	}










// GAP










void write_stats(char *file, int rnx, int rny)
	{
	ofstream out;
	int i,j,k;
	cout<<"writing stats to file "<<file<<endl;
	out.open(file);

	//out<<"WHY ISN'T THIS WORKING???\n"	
	out<<"\n# statistics!\n";
	out<<"\n#The "<<NPHIS<<" phis are: ";
	for (i=0; i<NPHIS; i++)
		{
		out<<PHIS[i]*todeg<<" ";
		}
	
	out<<"\n\n#The "<<NTHETAS<<" thetas are: ";
	for (i=0; i<NTHETAS; i++)
		{
		out<<THETAS[i]*todeg<<" ";
		}
	
	
	out<<"\n\n#The "<<NFREQS<<" frequencies in GHz are: ";
	for (i=0; i<NFREQS; i++)
		{
		out<<FREQS[i]/1e9<<" ";
		}
	out<<"\n";

	//Calculating total time
	float time_total = TIMES[0] + TIMES[4] + TIMES[2] + TIMES[3];	
	float total_time_seconds = TIMES[0] + TIMES[2] + TIMES[3] + TIMES[4];
	float total_time_minutes = total_time_seconds / 60;
	float total_time_hours = total_time_minutes / 60;	
	out<<"#times:\n";
//	TIMES[4]=time(0)-TIMES[4];
	out<<"The size of this surface was "<<rnx<<"x"<<rny<<"m\n";
	out<<"The maximum number of track subdivisions was: "<<MAX_TDIV<<", and the max number of facet divisions was "<<MAX_FDIV<<"^2\n\n";
	out<<"Surface initialisation took "<<TIMES[0]<<" seconds.\n";
	out<<"Track initialisation took "<<TIMES[2]<<" seconds.\n";
	out<<"Other initialisation took "<<TIMES[3]<<" seconds.\n";
	out<<"Cherenkov calculations took "<<TIMES[4]<<" seconds.\n";
	out<<"\n";

	//Eloquent way to output calculation time. I'm pretty sure there's a more
	//efficient way to do this, but it's above my pay grade.
	if(total_time_hours < 1) {
		if(fmod(total_time_minutes,60)) {
			float min_remainder = total_time_minutes-floor(total_time_minutes);
			float sec_int = min_remainder * 60;	
			out<<"Total simulation took "<<floor(total_time_minutes)<<" minutes and "<<sec_int<<" seconds.\n";
		}
		else {
			out<<"Total simulation took "<<total_time_minutes<<" minutes.\n";
		}
	}
	else {
		if(fmod(total_time_hours,60)) {
			float hour_remainder = total_time_hours-floor(total_time_hours);
			float min_remainder = hour_remainder * 60;
			float min_int = floor(min_remainder);
			float min_remainder2 = min_remainder-floor(min_remainder);
			float sec_int = min_remainder2 * 60;	
			out<<"Total simulation took "<<floor(total_time_hours)<<" hours, "<<min_int<<" minutes and "<<sec_int<<" seconds. Damn.\n";
		}
		else {
			out<<"Total simulation took "<<total_time_hours<<" hours. Sheesh.\n";
		}
	}
	out<<"\n";
	
	out<<"Maximum 'curvature' MAX_CURVE = "<<MAX_CURVE<<"\n";
	out<<"Max num track divs = "<<MAX_TDIV<<"\n";
	out<<"Max num facet divs = "<<MAX_FDIV<<"\n";
	out<<"Facets were split "<<NSPLITCALLS<<" times.\n";
	
	out<<"\n\n\n\n\n#Visibility stats:\n";
	out<<"# Theta    phi   Good    Badcos    Blocked\n";
	for (i=0; i<NTHETAS; i++)
		{
		for (j=0; j<NPHIS; j++)
			{
			out<<THETAS[i]<<"  "<<PHIS[j];
			for (k=0; k<3; k++)
				{
				out<<"      "<<VIS_STATS[i][j][k];
				}
			out<<"\n";
			}
		//out<<"\n"; ///APRAENTLY THIS DOESN'T DO ANYTHING?
		}
	out.close();
	cout<<"done"<<endl;
	}

void write_all_outputs(int freqyn, char *freq_file, int thetayn, char *theta_file, int phiyn, char *phi_file, int timeyn, char *time_file)
	{
	int i,j,k;
	ofstream out;
	long double t1,t2, I, Q, U, V;
	long double ****stokes;
	char suffix[5][10]={"_I.out", "_Q.out", "_U.out", "_V.out", "_XYZ.out"};
	char tempfreq[mfl], tempphi[mfl], temptheta[mfl];
	cout<<"writing outputs "<<freqyn<<thetayn<<phiyn<<timeyn<<endl;
	if (timeyn ==1)
		{
		cout<<"Generating time-domain pulses..."<<endl;
		gen_time_domain(time_file);
		cout<<"    ...finished."<<endl;
		}
if (thetayn == 1 || phiyn == 1 || freqyn == 1)
	{
	stokes = new long double***[NTHETAS];
	for (i=0; i<NTHETAS; i++)
		{
		stokes[i]=new long double **[NPHIS];
		for (j=0; j<NPHIS; j++)
			{
			stokes[i][j] = new long double *[NFREQS];
			for (k=0; k<NFREQS; k++)
				{
				stokes[i][j][k] = new long double[4];
// ####### Stokes IQUV ###########
				t1=EZ[i][j][2*k]*EZ[i][j][2*k] + EZ[i][j][2*k+1]*EZ[i][j][2*k+1];
				t2=EXY[i][j][2*k]*EXY[i][j][2*k] + EXY[i][j][2*k+1]*EXY[i][j][2*k+1];
				I = t1+t2;
				Q= (t1-t2)/I; // all others are fractional
				U= 2*(EZ[i][j][2*k]*EXY[i][j][2*k] + EZ[i][j][2*k+1]*EXY[i][j][2*k+1])/I;
				V= 2*(EZ[i][j][2*k+1]*EXY[i][j][2*k] - EZ[i][j][2*k]*EXY[i][j][2*k+1])/I;
				// note: V=0 implies linearity.
				stokes[i][j][k][0]=I;
				stokes[i][j][k][1]=Q;
				stokes[i][j][k][2]=U;
				stokes[i][j][k][3]=V;
				}
			}
		}
	}
	
	for (i=0; i<4; i++)
		{
		strcpy(tempfreq, freq_file);
		strcat(tempfreq, suffix[i]);
		strcpy(tempphi, phi_file);
		strcat(tempphi, suffix[i]);
		strcpy(temptheta, theta_file);
		strcat(temptheta, suffix[i]);
		if (freqyn == 1) {write_freqs(tempfreq, stokes, i);}
		if (thetayn == 1) {write_thetas(temptheta, stokes, i);}
		if (phiyn == 1) {write_phis(tempphi, stokes, i);}
		}
	i=4;
		{
		strcpy(tempfreq, freq_file);
		strcat(tempfreq, suffix[i]);
		strcpy(tempphi, phi_file);
		strcat(tempphi, suffix[i]);
		strcpy(temptheta, theta_file);
		strcat(temptheta, suffix[i]);
		if (freqyn == 1) {write_freqs(tempfreq, EZ, EXY);}
		if (thetayn == 1) {write_thetas(temptheta, EZ, EXY);}
		if (phiyn == 1) {write_phis(tempphi, EZ, EXY);}
		}
	}
		
void write_freqs(char *file, long double ****stokes, int which)
	{
	ofstream out;
	out.open(file);
	int i,j,k;
	
	out<<"# Outputting data by frequency for Stokes "<<which<<"\n";
	out<<"# We have "<<NFREQS<<" frequency-specific tables, each of\n";
	out<<"# which contains nphi = "<<NPHIS<<" rows by ntheta= "<<NTHETAS<<" columns;";
	out<<"# The phi are: ";
	for (i=0; i<NPHIS; i++) {out<<PHIS[i]*todeg<<" ";}
	out<<"\n# The theta are: ";
	for (i=0; i<NTHETAS; i++) {out<<THETAS[i]*todeg<<" ";}
	out<<"\n\n\n\n";
	
	for (i=0; i<NFREQS; i++)
		{
		out<<"# All theta - phi for frequency "<<FREQS[i]<<"\n";
		for (j=0; j<NPHIS; j++)
			{
			for (k=0; k<NTHETAS; k++)
				{
				out<<stokes[k][j][i][which]<<" ";
				}
			out<<"\n";
			}
		out<<"\n\n";
		}
	out.close();
	}
	
void write_freqs(char *file, long double ***EZ, long double ***EXY)
	{
	ofstream out;
	out.open(file);
	int i,j,k;
	
	out<<"# Outputting data by frequency for R(EZ), I(EZ), R(EXY), I(EXY)\n";
	out<<"# We have "<<NFREQS<<" frequency-specific tables, each of\n";
	out<<"# which contains nphi = "<<NPHIS<<" rows by ntheta= "<<NTHETAS<<" columns;";
	out<<"# The phi are: ";
	for (i=0; i<NPHIS; i++) {out<<PHIS[i]*todeg<<" ";}
	out<<"\n# The theta are: ";
	for (i=0; i<NTHETAS; i++) {out<<THETAS[i]*todeg<<" ";}
	out<<"\n\n\n\n";
	
	for (i=0; i<NFREQS; i++)
		{
		out<<"# All theta - phi for frequency "<<FREQS[i]<<"\n";
		for (j=0; j<NPHIS; j++)
			{
			for (k=0; k<NTHETAS; k++)
				{
				out<<EZ[k][j][2*i]<<" "<<EZ[k][j][2*i+1]<<" "<<EXY[k][j][2*i]<<" "<<EXY[k][j][2*i+1]<<"         ";
				}
			out<<"\n";
			}
		out<<"\n\n";
		}
	out.close();
	}
	
void write_thetas(char *file, long double ****stokes, int which)
	{
	ofstream out;
	out.open(file);
	int i,j,k;
	
	out<<"# Outputting data by theta for Stokes "<<which<<"\n";
	out<<"# We have "<<NTHETAS<<" theta-specific tables, each of\n";
	out<<"# which contains nphi = "<<NFREQS<<" rows by nfreqs = "<<NPHIS<<" columns;";
	out<<"# The phi are: ";
	for (i=0; i<NPHIS; i++) {out<<PHIS[i]*todeg<<" ";}
	out<<"\n# The freqs are: ";
	for (i=0; i<NFREQS; i++) {out<<FREQS[i]/1e6<<" ";}
	out<<" (MHz)\n\n\n\n";
	
	for (i=0; i<NTHETAS; i++)
		{
		out<<"# All phi-freq for theta "<<THETAS[i]*todeg<<"\n";
		for (j=0; j<NFREQS; j++)
			{
			out<<FREQS[j]<<" ";
			for (k=0; k<NPHIS; k++)
				{
				out<<" "<<stokes[i][k][j][which];
				}
			out<<"\n";
			}
		out<<"\n\n";
		}
	out.close();
	}

void write_thetas(char *file, long double ***EZ, long double ***EXY)
	{
	ofstream out;
	out.open(file);
	int i,j,k;
	
	out<<"# Outputting data by theta for R(EZ), I(EZ), R(EXY), I(EXY)\n";
	out<<"# We have "<<NTHETAS<<" theta-specific tables, each of\n";
	out<<"# which contains nphi = "<<NFREQS<<" rows by nfreqs = "<<NPHIS<<" columns;";
	out<<"# The phi are: ";
	for (i=0; i<NPHIS; i++) {out<<PHIS[i]*todeg<<" ";}
	out<<"\n# The freqs are: ";
	for (i=0; i<NFREQS; i++) {out<<FREQS[i]/1e6<<" ";}
	out<<" (MHz)\n\n\n\n";
	
	for (i=0; i<NTHETAS; i++)
		{
		out<<"# All phi-freq for theta "<<THETAS[i]*todeg<<"\n";
		for (j=0; j<NFREQS; j++)
			{
			out<<FREQS[j]<<" ";
			for (k=0; k<NPHIS; k++)
				{
				out<<EZ[i][k][2*j]<<" "<<EZ[i][k][2*j+1]<<" "<<EXY[i][k][2*j]<<" "<<EXY[i][k][2*j+1]<<"         ";
				}
			out<<"\n";
			}
		out<<"\n\n";
		}
	out.close();
	}
	
void write_phis(char *file, long double ****stokes, int which)
	{
	ofstream out;
	out.open(file);
	int i,j,k;
	
	out<<"# Outputting data by phi for Stokes "<<which<<"\n";
	out<<"# We have "<<NPHIS<<" phi-specific tables, each of\n";
	out<<"# which contains nthetas = "<<NTHETAS<<" rows by nfreqs = "<<NFREQS<<" columns;";
	out<<"# The thetas are: ";
	for (i=0; i<NTHETAS; i++) {out<<THETAS[i]*todeg<<" ";}
	out<<"\n# The freqs are: ";
	for (i=0; i<NFREQS; i++) {out<<FREQS[i]/1e6<<" ";}
	out<<" (MHz)\n\n\n\n";
	
	for (i=0; i<NPHIS; i++)
		{
		out<<"# All phi-freq for phi "<<PHIS[i]*todeg<<"\n";
		for (j=0; j<NTHETAS; j++)
			{
			out<<THETAS[j]*todeg<<" ";
			for (k=0; k<NFREQS; k++)
				{
				out<<" "<<stokes[j][i][k][which];
				}
			out<<"\n";
			}
		out<<"\n\n\n\n";
		}
	out.close();
	}	

void write_phis(char *file, long double ***EZ, long double ***EXY)
	{
	ofstream out;
	out.open(file);
	int i,j,k;
	
	out<<"# Outputting data by phi for R(EZ), I(EZ), R(EXY), I(EXY)\n";
	out<<"# We have "<<NPHIS<<" phi-specific tables, each of\n";
	out<<"# which contains nthetas = "<<NTHETAS<<" rows by nfreqs = "<<NFREQS<<" columns;";
	out<<"# The thetas are: ";
	for (i=0; i<NTHETAS; i++) {out<<THETAS[i]*todeg<<" ";}
	out<<"\n# The freqs are: ";
	for (i=0; i<NFREQS; i++) {out<<FREQS[i]/1e6<<" ";}
	out<<" (MHz)\n\n\n\n";
	
	for (i=0; i<NPHIS; i++)
		{
		out<<"# All phi-freq for phi "<<PHIS[i]*todeg<<"\n";
		for (j=0; j<NTHETAS; j++)
			{
			out<<THETAS[j]*todeg<<" ";
			for (k=0; k<NFREQS; k++)
				{
				out<<EZ[j][i][2*k]<<" "<<EZ[j][i][2*k+1]<<" "<<EXY[j][i][2*k]<<" "<<EXY[j][i][2*k+1]<<"         ";
				}
			out<<"\n";
			}
		out<<"\n\n";
		}
	out.close();
	}	
	
	

// assumes x is slow and y is fast, i.e. index is [ix*ny+iy]
// naveraged is the step size, rnx and rnx is the offset size
// we also need something that treats only the central nxm or so
void mod_surface(float *rsurface, int rnx, int rny, int naveraged, int nosx, int nosy, float *&surface, long double rdd, int Nconcat)
	{
	int i,j,ix,iy;
	NX = (rnx-1-nosx-2*Nconcat)/naveraged+1;
	NY = (rny-1-nosy-2*Nconcat)/naveraged+1;
	surface = new float[NX*NY];
	
	ix=0;
	for (i=nosx+Nconcat; i<rnx-Nconcat; i+=naveraged)
		{
		iy=0;
		for (j=nosy+Nconcat; j<rny-Nconcat; j += naveraged)
			{
			surface[ix*NY+iy]=rsurface[i*rny+j];
			iy++;
			}
		ix++;
		}
	DD = naveraged*rdd;
	X=(NX-1)*DD;
	Y=(NY-1)*DD;
	X2=X/2.;
	Y2=Y/2.;
	}
	
void gen_theta(long double theta_0, long double theta_1)
	{
	int i;
	long double dtheta;
	THETAS = new long double[NTHETAS];
	if (NTHETAS > 1)
		{
		dtheta=(theta_1-theta_0)/(NTHETAS-1);
		for (i=0; i<NTHETAS; i++)
			{
			THETAS[i]=theta_0+i*dtheta;
			}
		}
	else
		{
		THETAS[0]=(theta_0+theta_1)/2.;
		}
	
	}


void gen_phi(long double phi0,  int hnphi, long double dphi)
	{
	int i;
	PHIS=new long double[NPHIS];
	for (i=1; i<hnphi+1; i++)
		{
		PHIS[hnphi-i]=-dphi*i+phi0;
		PHIS[hnphi+i]=dphi*i+phi0;
		}
	PHIS[hnphi]=phi0;
	
	}

void gen_fs(long double f0, long double f1, int logyn)
	{
	long double f,df;
	int i;
	FREQS = new long double[NFREQS];
	
	if (logyn==1)
		{
		f0=log10(f0);
		f1=log10(f1);
		}
	if (NFREQS > 1)
		{
		df=(f1-f0)/(NFREQS-1);
		f=f0;
		for (i=0; i<NFREQS; i++)
			{
			FREQS[i]=f;
			f += df;
			}
		if (logyn==1)
			{
			for (i=0; i<NFREQS; i++)
				{
				FREQS[i]=pow(10.l,FREQS[i]);
				}
			}
		}
	else
		{
		FREQS[0]=(f0+f1)/2.;
		if (logyn==1) {FREQS[0] = pow(10.l,FREQS[0]);}
		}
	}

// ny nx since the xs vary fastest
void read_surface(char *surf_file, float *&surface, int &nsurfx, int &nsurfy, long double &ddistance, long double &H, long double &Srms)
	{
	ifstream in;
	char bit, tempstring[128];
	int i, size;
	
	in.open(surf_file);
	in>>bit;
	if (bit != '#') {cout<<"Incorrect format!"<<endl;}
	in>>nsurfy>>nsurfx>>H>>ddistance>>Srms;
	size=nsurfx*nsurfy;
	surface = new float[nsurfy*nsurfx];
	for (i=0; i<size; i++)
		{
		in>>surface[i];
		}
	}

// ny nx since the xs vary fastest
void read_float_surface(char *surf_file, float *&surface, int &nsurfx, int &nsurfy, long double &ddistance, long double &H, long double &Srms)
	{
	float temp;
	long int size;
	ifstream in;
	in.open(surf_file);
	in.read((char *) &nsurfx, sizeof(int));
	in.read((char *) &nsurfy, sizeof(int));
	in.read((char *) &temp, sizeof(float));
	H = (long double) temp;
	in.read((char *) &temp, sizeof(float));
	ddistance = (long double) temp;
	in.read((char *) &temp, sizeof(float));
	Srms = (long double) temp;
	
	cout<<"heading bits are: "<<nsurfx<<" "<<nsurfy<<" "<<H<<" "<<ddistance<<" "<<Srms<<"\n";
	
	size = nsurfx;
	size *= nsurfy;
	cout<<"size = "<<size<<endl;
	surface = new float[size];
	if (surface == 0) {cout<<"Cannot initialise surface - not enough memory!"<<endl;}
	else {in.read((char *) surface, sizeof(float)*size); cout<<"Surface read in."<<endl;}
	in.close();
	}


// surface is assumed to have centre 0,0, and dimesions X by Y = nx*dd + ny*dd
long double depth_to_z(long double dx, long double dy, long double depth, long double &dz, float *surf)
	{
	int x0,x1,y0,y1;
	long double h;
	
	x0=(int) ((X2+dx)/DD);
	x1=x0+1;
	y0=(int) ((Y2+dy)/DD);
	y1=y0+1;
	
	h=get_height(x0*DD, y0*DD, x1*DD, y1*DD,(long double) surf[x0*NY+y0], (long double) surf[x0*NY+y1], (long double) surf[x1*NY+y0], (long double) surf[x1*NY+y1],X2+dx, Y2+dy);
	dz=h-depth;
	
	cout<<"height at X="<<x0*DD<<", Y="<<y0*DD<<" is "<<dz<<"\n";
	
	return dz;
	}
	
// first finds points on x-lines, then interpolates in y
long double get_height(long double x0, long double y0, long double x1, long double y1, long double z00, long double z01, long double z10, long double z11, long double x, long double y)
	{
	long double z, kx,ky;
	long double dx10, dy10, kxp, kyp, dzx, dzy;
	dx10=1./(x1-x0);
	dy10=1./(y1-y0);
	kx=(x-x0)*(dx10);
	ky=(y-y0)*(dy10);
	kxp=1-kx;
	kyp=1-ky;
	
	z = kx*ky*z11 + kxp*ky*z10 + kx*kyp*z01 + kxp*kyp*z00; // labels are xy
	
	return z;
	}
	
void print_surface(char *surf_file, float *surface, int nx, int ny, long double H, long double esize, long double Srms, int app)
	{
	char tempstring[256], catstring[4], single[1];
	ofstream out;
	int i,j;
	single[0]='_';
	strcpy(tempstring,surf_file);
/*	strcat(tempstring,single);
	my_itoa(app, 4, catstring);
	strcat(tempstring, catstring); */
	
	out.open(tempstring, ios::app);
	out<<"#"<<nx<<" "<<ny<<" "<<H<<" "<<esize<<" "<<Srms<<"\n";
	for (i=0; i<nx; i++)
		{
		for (j=0; j<ny; j++)
			{
			out<<surface[i*ny+j]<<" ";
			}
		out<<"\n";
		}
	out<<"\n\n\n\n";
	out.close();
	}


// great initialisation routine
int init(ofstream &out, float *surface)
	{
	int i,j,k, size;
	int OK=1;
	long int memcount;
	
	memcount = sizeof(long double)*(NFREQS*7 + NPHIS*2 + NTHETAS*2 + NTHETAS*NPHIS*12 + NTHETAS*NPHIS*4*NFREQS);
	memcount += sizeof(int)*(NTHETAS*NPHIS*3);
	cout<<"The memory from basic initialisation is: "<<(memcount/1024./1024.)<<" MB\n";
	
	if (FMODE == 1) {if(!init_facets(surface, NX, NY)) {cout<<"Insufficient memory. Abort me now please!"<<endl; OK=0;}}
	
	// calculates error factors
	DIST_ERROR = ERROR_MARGIN;
	LAMBDA_ERROR = acos(1.-ERROR_MARGIN)/twopi;
	
	FACET_DIST_ERROR = DIST_ERROR;
	TRACK_DIST_ERROR = DIST_ERROR*32;
	FACET_LAMBDA_ERROR = LAMBDA_ERROR;
	TRACK_LAMBDA_ERROR = LAMBDA_ERROR*32;
	
	
	REFN_iont=REFNi/REFNt;
	REFN_toni=REFNt/REFNi;
	
	OPVECS = new long double**[NTHETAS]; if (OPVECS==0) {OK=0;}
	OPBASISZ = new long double**[NTHETAS]; if (OPBASISZ==0) {OK=0;}
	OPBASISXY = new long double**[NTHETAS]; if (OPBASISXY==0) {OK=0;}
	OPBASISK = new long double**[NTHETAS]; if (OPBASISK==0) {OK=0;}
	LAMBDAis = new long double[NFREQS]; if ( LAMBDAis==0) {OK=0;}
	LAMBDAts = new long double[NFREQS]; if ( LAMBDAts==0) {OK=0;}
	CHERENKOV_Cs = new long double[NFREQS]; if ( CHERENKOV_Cs==0) {OK=0;}
	OMEGAS = new long double[NFREQS]; if ( OMEGAS==0) {OK=0;}
	COSPHIS=new long double[NPHIS]; if ( COSPHIS==0) {OK=0;}
	SINPHIS=new long double[NPHIS]; if ( SINPHIS==0) {OK=0;}
	COSTHETAS=new long double[NTHETAS]; if ( COSTHETAS==0) {OK=0;}
	SINTHETAS=new long double[NTHETAS]; if ( SINTHETAS==0) {OK=0;}
	VIS_STATS = new int**[NTHETAS]; if ( VIS_STATS==0) {OK=0;}
	DS_TO_PHASES = new long double[NFREQS]; if ( DS_TO_PHASES==0) {OK=0;}
	FLENGTHS = new long double[NFREQS]; if ( FLENGTHS==0) {OK=0;}
	EZ = new long double**[NTHETAS]; if ( EZ==0) {OK=0;}
	EXY = new long double**[NTHETAS]; if ( EXY==0) {OK=0;}
	A_COEFS = new long double[NFREQS]; if ( A_COEFS==0) {OK=0;}
	
	for (i=0; i<NFREQS; i++)
		{
		LAMBDAis[i]=c_light/FREQS[i]/REFNi;
		LAMBDAts[i]=c_light/FREQS[i]/REFNt;
		CHERENKOV_Cs[i] = CHERENKOV_C*FREQS[i];
		OMEGAS[i] = FREQS[i]*twopi;
		DS_TO_PHASES[i] = twopi/LAMBDAts[i];
		FLENGTHS[i] = ALENGTHONHZ/FREQS[i];
		A_COEFS[i] = REFNt*FREQS[i]/c_light/2.;
		}
	for (i=0; i<NPHIS; i++)
		{
		COSPHIS[i]=cos(PHIS[i]);
		SINPHIS[i]=sin(PHIS[i]);
		}
	for (i=0; i<NTHETAS; i++)
		{
		OPVECS[i] = new long double*[NPHIS];  if ( OPVECS[i]==0) {OK=0;}
		OPBASISZ[i] = new long double*[NPHIS];  if ( OPBASISZ[i]==0) {OK=0;}
		OPBASISXY[i] = new long double*[NPHIS];  if ( OPBASISXY[i]==0) {OK=0;}
		OPBASISK[i] = new long double*[NPHIS];  if ( OPBASISK[i]==0) {OK=0;}
		VIS_STATS[i] = new int*[NPHIS]; if ( VIS_STATS[i]==0) {OK=0;}
		COSTHETAS[i]=cos(THETAS[i]);
		SINTHETAS[i]=sin(THETAS[i]);
		EZ[i] = new long double *[NPHIS]; if ( EZ[i]==0) {OK=0;}
		EXY[i] = new long double *[NPHIS]; if ( EXY[i]==0) {OK=0;}
		
		for (j=0; j<NPHIS; j++)
			{
			OPVECS[i][j]=new long double[3]; if ( OPVECS[i][j]==0) {OK=0;}
			OPBASISZ[i][j]=new long double[3]; if ( OPBASISZ[i][j]==0) {OK=0;}
			OPBASISXY[i][j]=new long double[3]; if ( OPBASISXY[i][j]==0) {OK=0;}
			OPBASISK[i][j]=new long double[3]; if ( OPBASISK[i][j]==0) {OK=0;}
			VIS_STATS[i][j] = new int[3]; if ( VIS_STATS[i][j]==0) {OK=0;}
			for (k=0; k<3; k++) {VIS_STATS[i][j][k]=0;}
			
			OPVECS[i][j][0]=COSPHIS[j]*COSTHETAS[i];
			OPVECS[i][j][1]=SINPHIS[j]*COSTHETAS[i];
			OPVECS[i][j][2]=SINTHETAS[i];
			
			// polarisation basis vector in XY plane
			OPBASISXY[i][j][0]=SINPHIS[j];
			OPBASISXY[i][j][1]=-COSPHIS[j];
			OPBASISXY[i][j][2]=0.;
			
			// polarisation basis vector in 'z' plane; i.e. maximal upwards pointing
			OPBASISZ[i][j][0]=-COSPHIS[j]*SINTHETAS[i];
			OPBASISZ[i][j][1]=-SINPHIS[j]*SINTHETAS[i];
			OPBASISZ[i][j][2]=COSTHETAS[i];
			
			OPBASISK[i][j][0]=COSPHIS[j]*COSTHETAS[i];
			OPBASISK[i][j][1]=SINPHIS[j]*COSTHETAS[i];
			OPBASISK[i][j][2]=SINTHETAS[i];
			
			EZ[i][j]= new long double[2*NFREQS]; if ( EZ[i][j]==0) {OK=0;}
			EXY[i][j] = new long double[2*NFREQS]; if ( EXY[i][j]==0) {OK=0;}
			for (k=0; k<NFREQS*2; k++)
				{
				EZ[i][j][k]=0.;
				EXY[i][j][k]=0.;
				}
			}
		}
	
	// we define the equation of each plane to be Ax+By+Cz=0, since this makes things easier!
	
	size=NTHETAS*NPHIS*(NX-1)*(NY-1);
	FPDONE = new bool[(NX-1)*(NY-1)]; // tracks whether visibility has yet been calculated
	if (FPDONE==0) {OK=0;}
	FPVISIBLE=new bool[size]; // tracks whether something is visible
	if (FPVISIBLE==0) {OK=0;}
	
	for (i=0; i<NFACETS; i++)
		{
		FPDONE[i]=0;
		}
	
	// divides frequencies into bins with ranges no greater than factors of 2.
// this means not too much time is wasted dividing tracks and facets
// this assumes of course that the limiting calculation is that of
// the Cherenkov radiation - if not, then this will be costing time!
	if (REFNt < REFNi)
		{
		TIR=asin(REFNt/REFNi);
		}
	else
		{
		TIR=piontwo;
		}
	
	if (f_division(out) != 1) {OK=0;}
	return OK;
	}


// assumes that the frequencies are ordered from lowest to highest
// the idea is that the facets will be subdivided according to near-field tests
// we want to minimise the subdivision at low frequencies if possible
// note that this could easily be ignored at a cost of only maybe a factor of 2 in run time.
int f_division(ofstream &out)
	{
	int i, i0, j, div, count, OK=1;
	long double last_f;
	
	NFDIVS=1;
	last_f=FREQS[0];
	count=1;
	for (i=1; i<NFREQS; i++)
		{
		if ((FREQS[i] > FDIV_FACTOR*last_f))
			{
			NFDIVS++;
			last_f = FREQS[i];
			count=1;
			}
		else {count++;}
		}
	
	NPERFDIV=new int[NFDIVS]; if (NPERFDIV == 0) {OK=0; return OK;}
	FDIV_STARTS=new int[NFDIVS]; if (FDIV_STARTS == 0) {OK=0; return OK;}
	last_f=FREQS[0];
	i0=0;
	div=0;
	count=1;
	FDIV_STARTS[0]=0;
	for (i=1; i<NFREQS; i++)
		{
		if ((FREQS[i] > FDIV_FACTOR*last_f))
			{
			NPERFDIV[div]=i-i0;
			last_f = FREQS[i];
			i0=i;
			div++;
			FDIV_STARTS[div]=i;
			count=1;
			}
		else {count++;}
		}
	i=NFREQS-1;
	NPERFDIV[div]=i-i0+1;
	if (NFREQS == 1) {NPERFDIV[0]=1;}
	cout<<"Frequencies are divided as follows:\n";
	for (i=0; i<NFDIVS; i++)
		{
		out<<"#Division "<<i<<" with "<<NPERFDIV[i]<<" frequencies:\n    ";
		cout<<"    Division "<<i<<" with "<<NPERFDIV[i]<<" frequencies.\n";
		for (j=0; j<NPERFDIV[i]; j++)
			{
			out<<" "<<FREQS[FDIV_STARTS[i]+j];
			}
		out<<"\n";
		}
	return OK;
	}

void gen_time_domain(char *time_file)
	{
	fftw_plan planxy, planz;
	fftw_complex *fz, *fxy;
	double *tz, *txy, dt, norm;
	int i, j, k, nfreq, ntime;
	ofstream out;
	cout<<"Beginning time-domain transform"<<endl;
	dt = 1./(FREQS[NFREQS-1]-FREQS[0]);
	nfreq=NFREQS;
	norm = dt;
	ntime = (NFREQS-1)*2;
	
	tz = (double *) fftw_malloc(sizeof(double)*ntime);
	txy = (double *) fftw_malloc(sizeof(double)*ntime);
	fz = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nfreq);
	fxy = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nfreq);
	
	planz = fftw_plan_dft_c2r_1d(nfreq, fz, tz, FFTW_ESTIMATE);
	planxy = fftw_plan_dft_c2r_1d(nfreq, fxy, txy, FFTW_ESTIMATE);
	
	out.open(time_file);
	out<<"# Gives the time-domain voltages at every calculated angle\n";
	
	for (i=0; i<NTHETAS; i++)
		{
		for (j=0; j<NPHIS; j++)
			{
			for (k=0; k<NFREQS; k++)
				{
				fz[k][0] = (double) (EZ[i][j][2*k]);
				fz[k][1] = (double) (EZ[i][j][2*k+1]);
				fxy[k][0] = (double) (EXY[i][j][2*k]);
				fxy[k][1] = (double) (EXY[i][j][2*k+1]);
				}
			fftw_execute(planz);
			fftw_execute(planxy);
			out<<"\n\n\n\n";
			out<<"# time-domain signal for theta "<<i<<"(="<<THETAS[i]*todeg<<") and PHI "<<j<<" (="<<PHIS[j]*todeg<<").\n";
			out<<"# time,   Ez,    Exy \n";
			for (k=0; k<NFREQS; k++)
				{
				out<<k*dt<<" "<<tz[k]*norm<<" "<<txy[k]*norm<<"\n";
				}
			}
		}
	out.close();
	cout<<"done"<<endl;
	}
