#pragma once

/*

We must loop through:
facets (2)
tracks (1)
exit angles (2)
frequencies (1)

If we loop through angles first, then we can:
1: pick a facet
2: calculate for this facet at
3: work 'uphill' to the max facet
4: work back 'downhill' until some limit is reached
5: repeat for other angles


If we loop through facets, we get the benefit of keeping all the geometry,
but it's harder to make shortcuts on the limit of calculations.
Oh well, bugger it...
*/

#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <complex>

#include "../../Utilities/long_constants.h"
#include "/usr/local/include/fftw3.h"

using namespace std;
struct track
	{
	long double x0,y0,z0,t0,x1,y1,z1,t1, q;
	long double p[2][3];
	long double vec[3], delt;
	long double centre[4];
	long double length, hl;
	long double beta;
	};

// holds all the info about a facet you could ever want to know!
struct facet
	{
	long double p[4][3];
	long double centre[3];
	long double x0,x1,y0,y1,z00,z10,z01,z11;
	long double norm[3], xvec[3], yvec[3];
	long double area, dx,dy, sqrtA;
	long double xbasis[3], ybasis[3], xproj, yproj;
	};

long double XBASIS[3]={1,0,0}, YBASIS[3]={0,1,0}, ZBASIS[3]={0,0,1};

// which basis to use correctly and which to screw up
int USE_WHICH=2;
long double SIN_APPROX=1.e-7, SIN_APPROX_2=1e-7;
int FDIV_METH=0, TDIV_METH=0;
bool TDIV_CONST=false;

//############## Constants #############

const float MAX_SYS_MB=1024.; // max MB allowed in system memory (approx)
const int mfl=256; // max number of characters for file length
int MAX_SPLITS=2; // max number of times to split a subfacet due to being near the Cherenkov angle
long double ERROR_MARGIN=0.01; // approximate global error margin
long double LAMBDA_ERROR, DIST_ERROR; // error values for phase and distance to facets
long double LIMIT_EST=1e-5; // approximates sin theta/theta = 1 if theta < this
long double FACET_LAMBDA_ERROR, TRACK_LAMBDA_ERROR, FACET_DIST_ERROR, TRACK_DIST_ERROR; // applied to each of facet and track seperately

// Output variables
ofstream SURF_file;
int *SURF_intracks, *SURF_infreqs, SURF_nsurfs;
char *SURF_surfnames;
bool WRITESURF;

// ######### Global Variables ##########

long double X,Y,X2,Y2,DD;
int NX, NY;
long double REFNi, REFNt, ALENGTHONHZ, TIR; // to be entered in by the user later
long double REFN_iont, REFN_toni;

// shortcuts related to theta and phi
long double *THETAS, *PHIS, *COSTHETAS, *SINTHETAS, *COSPHIS, *SINPHIS;
int NTHETAS, NPHIS, NOPS;
long double ***OPVECS;
long double *DS_TO_PHASES;
long double ***OPBASISZ, ***OPBASISXY, ***OPBASISK; // basic vectors for output polarisations

// outputs!!!!
long double ***EZ;
long double ***EXY;

bool *FPDONE, *FPVISIBLE;

long double CHERENKOV_C=(e_charge)/(epsilon_0 *c_light*c_light);
long double *CHERENKOV_Cs;

// these are for fresnel zones
int ***FXCENTRES, ***FYCENTRES;
// related to trackss
long double TRACKS_C[3], TRACKS_SUM;
int NTRACKS;

// related to frequency
long double *FLENGTHS, *FREQS, *LAMBDAis, *LAMBDAts, *OMEGAS, *A_COEFS;
int NFREQS, NFDIVS, *FDIV_STARTS, *NPERFDIV;

// related to facets
int NFX, NFY, NFACETS;
facet **FACETS;
int ***VIS_STATS;

// for lookup table
/*
long double ****LU_TABLE;
long double LU_dtheta, LU_dphi, LU_dsize;
int LU_nsize, LU_ntheta, LU_nphi;
long double LU_MAX_SIZE, LU_MIN_SIZE, LU_MAX_THETA;
const int LU_NTHETA_INT=50, LU_NPHI_INT=50;
long int NLULOOPS;
long double LU_MAX_SDIF, LU_MAX_PDIF, LU_MAX_TDIF;
long double LOG_LU_MINS, LOG_LU_MAXS;
*/

// NOT INCCLUDED!!!!
// for determining fresnel zone
// calculated for each frequency range / output direction / facet
//bool ****INFZONE;
long double MAX_FT_DIST=512.;
bool in_zone(track atrack, facet afacet, long double fmax);

// ########## Facet splitting #########

long double DNTIRDIV=2.;
int NTIRDIV=2, NTIRDIV2=4;
int DEFAULT_FACETDIV=0;
facet *SUBFACETS=0;

// ########## Global Counters #########

int MAX_TDIV=0;
int MAX_FDIV=0;
long long int LOOP_COUNT=0;
long int NSPLITCALLS=0;
long double MAX_CURVE=0.;

//Charlie's stuff
long NTRACKDIV_COUNT = 0;

//switches: for output
int domags=0;
long double FDIV_FACTOR=1000.;
long long int N_VISIBLE=0, N_INVISIBLE=0;

int NTIMES=5, TIMES[5];

// ####### FUNCTION DECLARATIONS #######


// ##### from lookup_table.cpp ########

int gen_lookup_table();
long double get_norm_factors(long double phi_i, long double theta_i, long double dxyonlambda, long double &paranorm, long double &perpnorm);
long double calc_norm_factor_t(long double theta_i, long double phi_i, long double n_i, long double n_t, int ntheta, int nphi, long double dxyonlambda, long double &sum_para, long double &sum_perp);
long double calc_norm_factor_r(long double theta_i, long double phi_i, long double n_i, long double n_t, int ntheta, int nphi, long double dxyonlambda, long double &sum_para, long double &sum_perp);

// ###### from 'cherenkov_facets.cpp' #######

void gen_para_perp(long double Epol[3], long double normal[3], long double &fpara, long double &fperp, long double sin_theta_i);
void gen_Epol(long double tvec[3], long double propvec[3], long double polvec[3]);
long double calc_min_theta(long double x, long double y, int x0, int x1, int y0, int y1, long double pz, float *surf, long double phi, long double cosphi, long double sinphi, int nx, int ny);
void init_fp(float *surface, int isx, int isy, int ifacet);
bool is_plane_visible(int ifx, int ify, int itheta, int iphi, float *surface);
bool facet_to_plane(facet &afacet, int theta_i, int phi_j, long double &cos_theta_t, long double &dist, int ifx, int ify, float *surface);
void gen_vec3_norm_to_vec2_from_vec1(long double vec1[3], long double vec2[3], long double vec3[3]);
void gen_random_perp(long double ve1[3], long double vec2[3], long double vec3[3]);
void split_facet(track & subtrack, facet &subfacet, int if0, int nfs, int ifx, int ify, float *surface, int nsplits);
void cherenkov(int if0, int nfs, facet &subfacet, track &subtrack, int ifx, int ify, float *surface, int nsplits);

// #### functions from 'facet_ops.cpp' ####
void fill_facet_points(facet &afacet);
void get_facet_normal(facet &afacet);
void facet_normal(long double *p1, long double *p2, long double *p3, long double norm[3]);
void gen_subfacet(struct facet &facet, long double nxdivs, long double nydivs, int xth, int yth, struct facet &asubfacet);
//void init_facets(facet **&facets, float *surface, int nx, int ny);
int init_facets(float *surface, int nx, int ny);
void init_facet(float *surface, int nx, int ny, int i, int j);

// ### functions from 'track_ops.cpp' ####

int read_tracks(char *track_file, track *&tracks, int &ntracks);
void transform_tracks(long double theta, long double phi, long double dx, long double dy, long double dz, int ntracks, track *tracks);
void init_track(track &atrack);
void print_tracks(char *filename, track *tracks, int ntracks, int dowhat);
void gensubtrack(track &atrack, long double ntdiv, int nthdiv, track &subtrack);

// #### handling functions: from 'rough.cpp' #### (the ligaments of the program???)




// #### functions from 'basic_geom.cpp' ######
void cross(long double v1[3], long double v2[3], long double v3[3]);
long double get_height(long double x0, long double x1, long double z0, long double z1, long double x);
bool dtoplane(facet &afacet, int theta_i, int phi_j, long double &angle, long double &dist);
long double dtoplane(long double *point, long double *normal);
long double find_dist_zero_D(long double A, long double B, long double C, long double a, long double b, long double c);
long double dot(long double *vec1, long double *vec2);
long double calc_dist(long double p1[3], long double p2[3]);
long double calc_dist(long double x1, long double y1, long double z1, long double x2, long double y2, long double z2);
long double calc_dist(long double p1[3], long double p2[3], long double vec[3]);
void my_itoa(int number, int n, char *string);
void construct_norm_vector(long double p1[3], long double p2[3], long double vec[3]);
void construct_vector(long double p1[3], long double p2[3], long double vec[3]);
void normalise(long double vec[3]);
long double tan_theta(long double thisz, long double pz, long double hdist);
long double calc_dist(long double p1[4], long double p2[3], long double vec[3]);
long double my_atan_xy(long double y, long double x);
// initialisation functions



// ##### functions from 'rough.cpp' #######

int init(ofstream &out, float *surface);
void write_stats(char *file, int rnx, int rny);
void write_all_outputs(int freqyn, char *freq_file, int thetayn, char *theta_file, int phiyn, char *phi_file, int timeyn, char *time_file);
void write_freqs(char *file, long double ****stokes, int which);
void write_thetas(char *file, long double ****stokes, int which);
void write_phis(char *file, long double ****stokes, int which);
void write_freqs(char *file, long double ***EZ, long double ***EXY);
void write_thetas(char *file, long double ***EZ, long double ***EXY);
void write_phis(char *file, long double ***EZ, long double ***EXY);
//void init_theta_phi(long double ***&dirs, int ntheta, int nphi, long double *thetas, long double *phis);
void gen_time_domain(char *time_file);

void get_facet_0(track &atrack, int &xi, int &yi);
void loop_facets(track &track, int if0, int nfs, float *surface);
int f_division(ofstream &out);
void div_facet_track(track &atrack, facet &afacet, int if0, int nfs, int ifx, int ify, float *surface);
void read_surface(char *surf_file, float *&surface, int &nsurfx, int &nsurfy, long double &ddistance, long double &H, long double &Srms);
void read_float_surface(char *surf_file, float *&surface, int &nsurfx, int &nsurfy, long double &ddistance, long double &H, long double &Srms);
long double depth_to_z(long double dx, long double dy, long double depth, long double &dz, float *surf);
long double get_height(long double x0, long double y0, long double x1, long double y1, long double z00, long double z01, long double z10, long double z11, long double x, long double y);
void mod_surface(float *rsurface, int rnx, int nry, int naveraged, int nosx, int nosy, float *&surface, long double dd, int Nconcat);


void print_surface(char *surf_file, float *surface, int nx, int ny, long double H, long double esize, long double Srms, int app);
void my_itoa(int number, int n, char *string);
void gen_fs(long double f0, long double f1, int logyn);
void gen_phi(long double phi0, int hnphi, long double dphi);
void gen_theta(long double theta_0, long double theta_1);


#include "cherenkov_facets.cpp"
#include "facet_ops.cpp"
#include "track_ops.cpp"
#include "basic_geom.cpp"
// #include "lookup_table.cpp"