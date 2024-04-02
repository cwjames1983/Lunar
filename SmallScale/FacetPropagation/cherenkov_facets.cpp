#include "omp.h"

// just gets the facet immediately above the track - obviously, this is pretty stupid
// in terms of a starting facet, but oh well!
void get_facet_0(track &atrack, int &ix0, int &iy0)
	{
	ix0 = (int) ((atrack.centre[0])/DD);
	iy0 = (int) ((atrack.centre[1])/DD);
	}


// loop this within a frequency bin, since there will likely be different cuts
// as a function of f;
void loop_facets(track &atrack, int if0, int nfs, float *surface)
	{
	int iy,ix, ix0, iy0;
	int sqrf, count;
	get_facet_0(atrack, ix0, iy0);
	
	sqrf = (int) (pow((long double) NFX*NFY, 0.5l)*3.12532l)+1;
	count=0;
	
	for (iy=0; iy<NFY; iy++)
		{
		for (ix=0; ix<NFX; ix++)
			{
			count++;
			if (count%sqrf == 0) {cout<<"       Up to facet "<<ix<<" "<<iy<<" (centre [x,y] = ["<<FACETS[ix][iy].centre[0]<<" , "<<FACETS[ix][iy].centre[1]<<"])"<<endl;}
			div_facet_track(atrack, FACETS[ix][iy], if0, nfs, ix, iy, surface);
			}
		}
	}


void split_facet(track &subtrack, facet &subfacet, int if0, int nfs, int ifx, int ify, float *surface, int nsplits)
	{
	int i,j,k,l,m;
	facet *subfacets;
	
	NSPLITCALLS++;
	

	subfacets = new facet[NTIRDIV2];
	for (i=0; i<NTIRDIV; i++)
		{
		for (j=0; j<NTIRDIV; j++)
			{
			gen_subfacet(subfacet, DNTIRDIV, DNTIRDIV, i,j, subfacets[i*NTIRDIV+j]);
			}
		}
	
	for (j=0; j<NTIRDIV; j++)
		{
		for (k=0; k<NTIRDIV; k++)
			{
			if (WRITESURF)
				{
				for (l=0; l<4; l++)
					{
					for (m=0; m<3; m++)
						{
						SURF_file<<subfacets[j*NTIRDIV+k].p[l][m]<<" ";
						}
					}
				SURF_file<<"\n";
				}
			// put Cherenkov calculation as a function of frequencies here
			cherenkov(if0, nfs, subfacets[j*NTIRDIV+k], subtrack, ifx, ify, surface, nsplits+1);
			}
		}
	delete [] subfacets; subfacets=0;
	}


//void get_max_divs(facet **facets, int nfx, int nfy, long double maxf, track *tracks, int ntracks


// the conditions imposed are:
// the maximum wavelength offsets 
void div_facet_track(track &atrack, facet &afacet, int if0, int nfs, int ifx, int ify, float *surface)
	{
	long double maxf, eff_lambda, base_distance, wave_curvature;
	long double vtemp[3], v0[3], eff_wavelength, dnfacetdiv, dntrackdiv;
	int i,j,k,l,m;
	int ntrackdiv, nfacetdiv;
	long double fxs[3], fys[3], fzs[5], txs[3], tys[3], tzs[3], tts[3];
	facet *subfacets=0;
	track subtrack, *usetrack;
	long double dist, ddist, dlambda, dr; //CMW 04-03-24 okay, maybe ddist is delta d.
	long double mindot=1e10, ddot;
	
	maxf=FREQS[if0+nfs-1];
	eff_wavelength=c_light/maxf/REFNi;
	
	base_distance = calc_dist(afacet.centre, atrack.centre);
	construct_norm_vector(afacet.centre, atrack.centre, v0);
	// from tracks extrema to facet centre
	ddist=0.;
	for (i=0; i<2; i++)
		{
		construct_vector(afacet.centre, atrack.p[i], vtemp);
		dist=pow(dot(vtemp, vtemp),0.5l);
		ddot=fabs(dot(v0,vtemp))/dist;
		if (ddot < mindot) {mindot=ddot;}
		if (fabs(dist - base_distance) > ddist) {ddist = fabs(dist-base_distance);} //CMW 04-03-24 yeah this looks like delta d. This line just makes ddist an absolute value.
		}
	dr=ddist/base_distance; // clearly delta r
	wave_curvature=(1.-mindot)*base_distance/eff_wavelength;
	ntrackdiv = (int) (max(dr/TRACK_DIST_ERROR, wave_curvature/TRACK_LAMBDA_ERROR))+1;
	
	if (TDIV_METH > 0) {ntrackdiv=TDIV_METH;} // sets number of track devisions: for hard-coded testing
	
	mindot = 1e9;
	// for facet extrema to track centre
	ddist=0.;
	
	for (i=0; i<4; i++)
		{
		construct_vector(atrack.centre, afacet.p[i], vtemp);
		dist=pow(dot(vtemp,vtemp),0.5l);
		ddot=fabs(dot(v0,vtemp))/dist;
		if (ddot < mindot) {mindot=ddot;}
		if (fabs(dist - base_distance) > ddist) {ddist = fabs(dist-base_distance);}
		}
	dr=ddist/base_distance;
	wave_curvature=(1.-mindot)*base_distance/eff_wavelength;
	nfacetdiv = (int) (max(dr/FACET_DIST_ERROR, wave_curvature/FACET_LAMBDA_ERROR))+1;
	
	if (FDIV_METH > 0) {nfacetdiv = FDIV_METH;}
	
	dnfacetdiv=(long double) nfacetdiv;
	dntrackdiv=(long double) ntrackdiv;
	
	// logs greatest number of divisions
	if ( TDIV_CONST )
		{
		
		ntrackdiv = MAX_TDIV;
		if( ntrackdiv == 0) {ntrackdiv = 1;}
			// line added to print the value of ntrackdiv 12/11/2023:
			//std::cout<<"ntrackdiv is now "<<ntrackdiv<<std::endl;
		}
	else if ( ntrackdiv > MAX_TDIV ) { 
		MAX_TDIV=ntrackdiv; 
			//std::cout<<"ntrackdiv was  > MAX_TDIV and ntrackdiv is now "<<ntrackdiv<<std::endl;
		}
	if (nfacetdiv > MAX_FDIV) {MAX_FDIV=nfacetdiv;}
	
	
	// ensures we don't re-create memory unless we have to
	// i.e. only increases number of subfacets if this is larger than previous
	if (nfacetdiv > DEFAULT_FACETDIV || DEFAULT_FACETDIV==0)
		{
		if (SUBFACETS) {delete [] SUBFACETS;} // delete old subfacets
		DEFAULT_FACETDIV = max(nfacetdiv, 2*DEFAULT_FACETDIV);
		SUBFACETS = new facet[DEFAULT_FACETDIV*DEFAULT_FACETDIV];
		cout<<"Increased size of DEFAULT_FACETDIV to "<<DEFAULT_FACETDIV<<"\n";
		}
	
	// This forces the size of every (sub)facet to be equal. We must of course be careful later to calculate the incident E-fields according to projected area (this is a real weighting).
	
	if (nfacetdiv == 1) {subfacets=&afacet;}
	else
	 	{
		subfacets = SUBFACETS;
		for (i=0; i<nfacetdiv; i++)
			{
			for (j=0; j<nfacetdiv; j++)
				{
				gen_subfacet(afacet, dnfacetdiv, dnfacetdiv, i,j, subfacets[i*nfacetdiv+j]);
				}
			}
		}
//CMW 02-04-2024 making changes for parallelisation here
#pragma omp parallel
{
    int thread_id =omp_get_thread_num();
	for (i=0; i<ntrackdiv; i++)
		{
		if (ntrackdiv == 1)
			{
			usetrack = &atrack;
			}
		else
			{
			gensubtrack(atrack, dntrackdiv, i, subtrack);
			usetrack = &subtrack;
			NTRACKDIV_COUNT++;
            std::cout<<"Hello from thread "<<thread_id<<"!\n";
			}
		for (j=0; j<nfacetdiv; j++)
			{
			for (k=0; k<nfacetdiv; k++)
				{
			/*	if (WRITESURF && i==0)
					{
					for (l=0; l<4; l++)
						{
						for (m=0; m<3; m++)
							{
							SURF_file<<subfacets[j*nfacetdiv+k].p[l][m]<<" ";
							}
						}
					SURF_file<<"\n";
					} */
				// put Cherenkov calculation as a function of frequencies here
				cherenkov(if0, nfs, subfacets[j*nfacetdiv+k], *usetrack, ifx, ify, surface,0);
				}
			}
		}
}
//	if (nfacetdiv != 1) {delete [] subfacets; subfacets=0;}
	}

// creates a random vector to use as a basis for calculating a perpendicular.
void gen_random_perp(long double vec1[3], long double vec2[3], long double vec3[3])
	{
	long double temp1[3];
	int i;
	for (i=0; i<3; i++)
		{
		temp1[i]=vec1[i];
		}
	temp1[0] += 0.28373683;
	temp1[1] -= 0.028272727;
	temp1[2] += 0.20938722;
	normalise(temp1);
	gen_vec3_norm_to_vec2_from_vec1(temp1,vec2,vec3);
	}

// constructs a vector in vec1, vec2 plane which is normal to vec2
// i.e.
void gen_vec3_norm_to_vec2_from_vec1(long double vec1[3], long double vec2[3], long double vec3[3])
	{
	long double t;
	int i;
	t=dot(vec1,vec2);
	if (t == 0.)
		{
		//cross(vec1,vec2,vec3);
		for (i=0; i<3; i++) {vec3[i]=vec1[i];}
		}
	else if (t ==1.)
		{
		gen_random_perp(vec1,vec2,vec3);
		}
	else
		{
		t=-1./t;
		for (i=0; i<3; i++)
			{
			vec3[i]=vec2[i] + t*vec1[i];
			}
		normalise(vec3);
		}
	}

/* this is the heart of the program.


Q:	How do we weight for solid angles and the like?
A:	It is implicit.
	We already project the outgoing radiation vectors as a function of angle.
	We weight each area element .dA by a cos(theta_i) cos(theta_t) factor,
		to give the projected areas both incoming and outgoing. While
		this seems like a long double-compensation, in fact we integrate
		over the unprojected area, which cancels one effect.
	The field strength is reduced by 1/R, where R is the distance to the plane.
		Since the outgoing radiation is generally proportional to the
		illuminated sqrt area, this cancels the 1/R effect.

*/



/* TOMORROW, HOW TO TREAT THIS? [wtf are these comments???]

#1: Use 'logic' to work out projection effects - use fresnel to get the field on the other side, then project!
#2: do it with n1=n2. Then fresnel coefs = 1. We still require projection factors of the fields though, both for incoming, and outgoing.
#3: divide by wavelength as the normalising factor!
#4: DON'T normalise for area - it should be an integral over this, so that it is essentially smooth between facets (i.e. integrals add)


*/








//#######################################################################################################################################################################################
//												CHERENKOV FUNCTION


void cherenkov(int if0, int nfs, facet &subfacet, track &subtrack, int ifx, int ify, float *surface, int nsplits)
	{
	// utility variables
	int i,j,k,l;
	
	// track-to-plane variables - independent of outgoing radiation
	long double dtf, tfvec[3]; // distance track to facet
	long double cos_theta_i, sin_theta_i, theta_i;
	
	// Cherenkov variables
	long double phase_constant, phase_constant_delt_on_2; // variables for the Cherenkov calculation
	long double omega_phase_bit, omega_phase_bit_delt_on_2; // variables for the Cherenkov calculation
	long double cos_tangle, sin_tangle; // sin and cos of angle to shower axis
	long double incident_imag_Es[NFREQS];
	long double decoherence; // not a phase actually, but decoherence factor w/in a track
	long double app_velocity; // apparent tracklength
	long double attenuation; // attenuation of the radiation in the medium
	
	// refraction variables
	long double dfp; // distance facet to plane
	long double polvec[3], slap_vec[3], poke_vec[3];
	long double fslap, fpoke;
	long double sin_alpha, cos_alpha; // electric fields on the transmitted side; specular!
	long double t_para_fresnel, t_perp_fresnel; // refraction coefficients for Fresnel refraction
	long double *spec_slap_vec, spec_poke_vec[3];
	long double a,b; // used for constructing spec_poke_vec
	long double phi_i, cos_phi_i, sin_phi_i; // angle of incident radiation wrst x and y
	
	// innermost loop variables
	long double cos_theta_t, sin_theta_t, cos_phi_t, sin_phi_t, sin_beta, cos_beta; // sin and cos of angle wrst transmission basis vectors
	long double phi_hat[3], theta_hat[3]; // basis vectors for electric field calculations
	long double A; // decoherence factor during A. = effective height x interference
	long double phase0, sin_phase0, cos_phase0;
	long double Cxy, Cz, Cptz, Cptxy, Cppz, Cppxy, Cstz, Cstxy, Cspz, Cspxy; // coefs tracing from poke/slap, to theta_phi, to z/xy
	long double proj_pz, proj_pxy, proj_tz, proj_txy; // projections of theta-phi basis to z-xy basis
	long double theta_poke, theta_slap, phi_poke, phi_slap; //coeffs from Stratton-Chu of poke-slap to Etheta-Ephi;
	long double ctcap1, capct, ktix, ktixf, ktiy, ktiyf, slap_bit, poke_bit; // miscellaneous shortcuts
	
	int ntir; // counts number of corners where radiation is totally internally reflected
	
// ##### We now calculate quantities relevant to the radiation emitted by the track #######
	ntir=0;
	if (nsplits < MAX_SPLITS)
		{
		// checks to see how many corners are totally internally reflected
		for (i=0; i<4; i++)
			{
			dtf=calc_dist(subtrack.centre, subfacet.p[i], tfvec);
			cos_theta_i=dot(tfvec, subfacet.norm);
			theta_i=acos(cos_theta_i);
			if (theta_i > TIR) {ntir++;}
			}
		// if they all are, then throw away the surface (doesn't work for very big facets and near-surface cascades and extremely long wavelengths - but who the *?*! would calculate that????. Umm... for cosmic rays... but this should already have been handled by div-facet-track
		if (ntir == 4)
			{
			return;
			}
		else if (ntir > 0)
			{
			split_facet(subtrack, subfacet,if0, nfs, ifx, ify, surface, nsplits);
			return;
			}
		else {/* we continue! */}
		}
	
	// gets the distance from the centre of the track to the facet. tfvec IS normalised!!!
	dtf=calc_dist(subtrack.centre, subfacet.centre, tfvec);
	cos_theta_i=dot(tfvec, subfacet.norm);
	theta_i=acos(cos_theta_i);
	if (theta_i > TIR) {if (nsplits < MAX_SPLITS) {cout<<"why wasn't this picked up earlier????"<<endl;} return;} // since then nothing will come out!
	
	 // tangle is track angle, i.e. theta in most Cherenkov parameterisations.
	cos_tangle=dot(subtrack.vec,tfvec);
	sin_tangle=pow(1.-cos_tangle*cos_tangle,0.5l);
	
	// generates the phase constant in the Cherenkov parameterisations
	phase_constant=(1.-cos_tangle*REFNi*subtrack.beta); // the omega minus k dot v
	phase_constant_delt_on_2=phase_constant*subtrack.delt/2.; // multiplied by delt/2.
	app_velocity=subtrack.beta*sin_tangle*c_light;
	
// ##### Quantities related to the A of radiation through the facet #####
	
	// incident radiation
	
	if (theta_i == 0.)
		{
		sin_theta_i = 0.;
		phi_i = 0.; // arbitrary!
		cos_phi_i=1.;
		sin_phi_i=0.;
		}
	else
		{
		sin_theta_i=pow(1.-cos_theta_i*cos_theta_i,0.5l);
		// angle relative to the projections of the x-plane and y-plane onto the facet
		phi_i = my_atan_xy(dot(tfvec,subfacet.xbasis), dot(tfvec, subfacet.ybasis));
		cos_phi_i=cos(phi_i);
		sin_phi_i=sin(phi_i);
		}
	
	// generates the polarisation vector of the incident radiation
	gen_vec3_norm_to_vec2_from_vec1(subtrack.vec, tfvec,polvec);
	
	if (dot(subtrack.vec, polvec) < 0.) // we define polarisation as being in the direction of the track (again, arbitrary)
		{
	//	cout<<"pol was -ve!!!"<<endl;
		for (i=0; i<3; i++) {polvec[i] *= -1.;}
		}
	
	// generates the incident 'pokey' vector for incoming radiation
	gen_vec3_norm_to_vec2_from_vec1(subfacet.norm, tfvec, poke_vec);
	if (dot(subfacet.norm, poke_vec) > 0.)
		{
	//	cout<<"poke was +ve!!!"<<endl; // should be -ve since this is the way stuff is derived
		for (i=0; i<3; i++) {poke_vec[i] *= -1.;}
		}
	
	// generates the 'slappy' basis vector for both incoming and outgoing radiation
	cross(tfvec, poke_vec, slap_vec);
	
	// generates projections of the polarisation vector for the slappy and pokey cases of incoming radiation
	// these have fpoke^2+fslap^2=1, i.e. are simply orthogonal basis components. It has been checked this
	// always squares to 1.
	fslap=dot(polvec,slap_vec);
	fpoke=dot(polvec,poke_vec);
	
	sin_alpha = sin_theta_i*REFN_iont; // sin of specular angle
	cos_alpha=pow(1.-sin_alpha*sin_alpha,0.5l); // guess what this is then?
	
	// the relative x-y angle does not change upon refraction
	sin_beta=sin_phi_i;
	cos_beta=cos_phi_i;
	
	// standard Fresnel coefficients for A, giving the electric field on the outgoing
	// side -- this is the field which gets integrated over.
	t_para_fresnel = 2 * cos_theta_i / (REFN_toni*cos_theta_i + cos_alpha);
	t_perp_fresnel = 2 * cos_theta_i / (cos_theta_i + REFN_toni*cos_alpha);
	
// ####### genrating the specular pokey and slappy vectors #########

	spec_slap_vec=slap_vec; // they are the same!
	
	if (sin_theta_i <= SIN_APPROX_2) // allows a small deviation to prevent problems later
		{
		cross(spec_slap_vec, subfacet.norm, spec_poke_vec); // outputs to spec_poke_vec
		}
	else
		{
		a=cos_alpha/sin_theta_i;
		b=-sin_alpha-a*cos_theta_i; // have made b -ve now
		for (i=0; i<3; i++) {spec_poke_vec[i]=a*tfvec[i]+b*subfacet.norm[i];}
	//	if (dot(subfacet.norm, spec_poke_vec) > 0.) {cout<<"dot > 0!\n";} it is never > 0! hurray!
		// the dot between the spec poke and normal should be -ve
		if (fabs(dot(spec_poke_vec,spec_poke_vec) - 1.) > SIN_APPROX_2)
			{
			cout<<"Arghh!!!!!!!!!!! SIN APPROX 2 checking... "<<fabs(dot(spec_poke_vec,spec_poke_vec)-1.)<< " "<<sin_theta_i<<endl;
			normalise(spec_poke_vec);
			}
		}
// #### Now we generate frequency-specific shortcuts, since frequency is the innermost loop of the E-field calculation ####
	for (k=if0; k<if0+nfs; k++)
		{
		attenuation = exp(-dtf/FLENGTHS[k]); // attenuation in the interaction medium
		
		if (fabs(phase_constant*OMEGAS[k]) < LIMIT_EST)
			{
			decoherence=subtrack.delt; // limit of formula as phase -> 0
			}
		else
			{
			omega_phase_bit=OMEGAS[k]*phase_constant;
			omega_phase_bit_delt_on_2 = OMEGAS[k]*phase_constant_delt_on_2;
			decoherence = 2*sin(omega_phase_bit_delt_on_2)/omega_phase_bit;
			}
		// includes: decoherence along the track, apparent tracklength, Cherenkov coefficients, 1/R factor to the surface, and attenuation in the medium, and projection factor sin tangle. Shoukd there be a cos theta i here?
		
	//	incident_imag_Es[k] = subtrack.q*CHERENKOV_Cs[k]*attenuation/dtf; // spherical wave of appropriate magnitude
		incident_imag_Es[k]=subtrack.q*decoherence*CHERENKOV_Cs[k]*app_velocity*attenuation/dtf; // true values!
		}
	
// ####### calculates these phase factors for all frequencies ##############
	for (i=0; i<NTHETAS; i++)
		{
		for (j=0; j<NPHIS; j++)
			{
			// determines if the radiation from the facet gets to the plane unimpeded
			// it will also return the distance from the facet to plane.
			if (facet_to_plane(subfacet,i,j, cos_theta_t, dfp, ifx, ify, surface))
				{		
// phase when it gets there, i.e. includes propagation contributions to phase are: distance facet to exit plane dfp as a function of wavelength on A, distance track to facet dtf as function of wavelength on incidence, and the initial time tts[1] due to the centre of the track (phases relative to this).
				N_VISIBLE++;
				
				// angle of the normal to the output vector
			//	cos_theta_t = dot(subfacet.norm, OPBASISK[i][j]);
				if (cos_theta_t < 0.) {cout<<"something is wrong!!!!!\n"<<endl;}
				if (cos_theta_t == 1.) {sin_theta_t = 0.;}
				else {sin_theta_t = pow(1.-cos_theta_t*cos_theta_t,0.5l);}
				
				// if the angle doesn't matter, then we just set them both to be sensible;
				if (sin_theta_t == 0.)
					{
					cos_phi_t=1.;
					sin_phi_t=0.;
					}
				else
					{
					// get phi from angle of non-z-normal to facet.xbasis
					cos_phi_t=dot(subfacet.xbasis, OPBASISK[i][j]);
					if (cos_phi_t != 0.) {cos_phi_t /=sin_theta_t;}
					sin_phi_t=dot(subfacet.ybasis, OPBASISK[i][j]);
					if (sin_phi_t != 0.) {sin_phi_t /=sin_theta_t;}
					}
				// constructs basis vector phi_hat and theta hat
				for (l=0; l<3; l++)
					{
					phi_hat[l] = -sin_phi_t*subfacet.xbasis[l] + cos_phi_t*subfacet.ybasis[l];
					theta_hat[l] = cos_theta_t*cos_phi_t*subfacet.xbasis[l] + cos_theta_t*sin_phi_t*subfacet.ybasis[l] - sin_theta_t*subfacet.norm[l];
					}
				
				// we calculate A only for n_i = n_t, since transmission is handled sperately	
				
				// optical lengths of x- and y-paths in the intergral over the phases
				// note that we use the refn_t since refnt*sin_alpha = sin_theta_t * refni
				ktix = 0.5*(sin_theta_t*cos_phi_t - sin_alpha*cos_phi_i);
				ktiy = 0.5*(sin_theta_t*sin_phi_t - sin_alpha*sin_phi_i);
				
				// we calculate based on basic vectors of \hat{phi} and \hat{theta}. But we need to eventually transform these into vectors of basis z and xy. Hence we must calculate coefficients which project theta hat and phi hat onto these.
				
				// we want single coefficients which trace incoming pokey and slappy vectors (2; p & s) to theta and phi vectors (2 also (t & p)) and then onto the output basic vectors (again, 2! (z & xy)). This generates 8 coefficients in total.
				poke_bit = t_para_fresnel*fpoke;
				slap_bit = t_perp_fresnel*fslap;
				
				ctcap1 = cos_theta_t*cos_alpha + 1.;
				capct = cos_theta_t+cos_alpha;
				
				// note that phi_i = beta here
				theta_slap = ctcap1*(sin_phi_t*cos_phi_i - cos_phi_t*sin_phi_i);
				phi_slap = capct*(sin_phi_t*sin_phi_i + cos_phi_t*cos_phi_i);
				theta_poke = phi_slap;
				phi_poke = - theta_slap;
				
				
				proj_pz = dot(phi_hat, OPBASISZ[i][j]);
				proj_pxy = dot(phi_hat, OPBASISXY[i][j]);
				proj_tz = dot(theta_hat, OPBASISZ[i][j]);
				proj_txy = dot(theta_hat, OPBASISXY[i][j]);
				
				// the labelling goes: p/s (pokey/slappy) t/p (theta/phi coefs) xy/z (projecting onto output directions)
				Cptz = poke_bit * theta_poke * proj_tz;
				Cptxy = poke_bit * theta_poke * proj_txy;
				Cppz = poke_bit * phi_poke * proj_pz;
				Cppxy = poke_bit * phi_poke * proj_pxy;
				Cstz = slap_bit * theta_slap * proj_tz;
				Cstxy = slap_bit * theta_slap * proj_txy;
				Cspz = slap_bit * phi_slap * proj_pz;
				Cspxy = slap_bit * phi_slap * proj_pxy;
				
				// adds the components of Cz and Cxy to get the correct magnitudes from incoming to outgoing.
				Cxy = Cptxy + Cppxy + Cstxy + Cspxy;
				Cz = Cptz + Cppz + Cstz + Cspz;
				
				
				// now looping through frequencies
				for (k=if0; k<if0+nfs; k++)
					{
					// total phase at this frequency. larger start time = earlier phase means omega t component should be +ve, not e^ikx-wt, and phase for track centre should be negative. If emitted late then it has less far to go.
					phase0 = twopi * (subtrack.centre[3]*FREQS[k] + dtf/LAMBDAis[k] + dfp/LAMBDAts[k]);
					LOOP_COUNT++;
					
					// now a phase, rather than an optical length
					ktixf = ktix*DS_TO_PHASES[k];// = twopi/LAMBDAtS[k]; (
					ktiyf = ktiy*DS_TO_PHASES[k]; //twopi/LAMBDAtS[k];
		
// this gives the magnitude, no projection effects: this is an E-field coefficient!!!
// (integrals over two exponentials are seperable; ddon2 is half the distance)
// note that as the phases become small, this scales with ddon2^2, i.e. the area. dd is to convert from basic electric field to a unit length. (field per area)
					if (fabs(ktixf*subfacet.dx) <= SIN_APPROX)
						{
						if (fabs(ktiyf*subfacet.dy) <= SIN_APPROX)
							{
							A = subfacet.area;
							}
						else
							{
							A = subfacet.dx*sin(ktiyf*subfacet.dy)/ktiyf;
							}
						}
					else
						{
						if (fabs(ktiyf*subfacet.dy) <= SIN_APPROX)
							{
							A = subfacet.dy*sin(ktixf*subfacet.dx)/ktixf;
							}
						else
							{
							// note that this is the only one which will ever occur, so this has the only shortcuts
							A = sin(ktixf*subfacet.dx) * sin(ktiyf*subfacet.dy) / (ktixf*ktiyf);
							}
						}
					
					sin_phase0=sin(phase0);
					cos_phase0=cos(phase0);
					
					A *= incident_imag_Es[k]*A_COEFS[k]; // A_COEFS are simply k on four pi
					
					
					// Go phases! They start off imaginary, then evolve with time to...
					EZ[i][j][2*k] += A * sin_phase0 * Cz; // real
					EXY[i][j][2*k] += A * sin_phase0 * Cxy;
					EZ[i][j][2*k+1] += A * cos_phase0 * Cz; // imaginary
					EXY[i][j][2*k+1] += A * cos_phase0 * Cxy;
					if (isnan(EZ[i][j][2*k]) || isnan(EZ[i][j][2*k+1]) || isnan(EXY[i][j][2*k]) || isnan(EXY[i][j][2*k+1]))
						{
						cout<<"isnan!"<<endl;
						}
					}
				}
			else
				{
				N_INVISIBLE++;
				}
			}
		}
	}

// ########### SET OF ROUTINES FOR CALCULATING VISIBILITY ###########


// calculates the distace from a facet to the output planes; only does a very simple "is it detectable?"
// ifx and ify identify which original facet (NOT subfacet!)
// angle is the angle from the facet normal to the plane, i.e. it is the theta_t
bool facet_to_plane(facet &afacet, int theta_i, int phi_j, long double &cos_theta_t, long double &dist, int ifx, int ify, float *surface)
	{
	bool OK;
	cos_theta_t=dot(afacet.norm, OPVECS[theta_i][phi_j]); // operates on a sub-plane level
	if (cos_theta_t < 0)
		{
		OK=0;
		VIS_STATS[theta_i][phi_j][1]++;
		}
	else
		{
		OK=is_plane_visible(ifx, ify, theta_i, phi_j, surface); // only on a facet level; *could* be subfacet, but unlikely (VERY long!)
		if (OK)
			{
			VIS_STATS[theta_i][phi_j][0]++;
			dist=dtoplane(afacet.centre, OPVECS[theta_i][phi_j]); // negative is OK.
			}
		else {VIS_STATS[theta_i][phi_j][2]++;}
		}
	return OK;
	}

// determines if the details of facet-plane distances have been calculated yet,
// and does so if not
bool is_plane_visible(int ifx, int ify, int itheta, int iphi, float *surface)
	{
	int index, ifacet;
	ifacet=ifx*NFY+ify;
	
	index=ifacet*NOPS+iphi*NTHETAS+itheta;
	if (!FPDONE[ifacet])
		{
		FPDONE[ifacet]=1;
		init_fp(surface, ifx, ify, ifacet);
		}
	return FPVISIBLE[index];
	}

// calculates it for all theta and phi for that facet
// takes it as yes or no regardless of the subfacet - ie based on actual facet
void init_fp(float *surface, int isx, int isy, int ifacet)
	{
	int i,j, ifx, ify;
	long double mintheta;
	if (FMODE == 1) {ifx = isx; ify=isy;} else {ifx=0; ify=0;}
	
	for (i=0; i<NPHIS; i++)
		{
		// we always take the mid-point
		mintheta=calc_min_theta(isy+0.5, isx+0.5, isy, isy+1, isx, isx+1, FACETS[ifx][ify].centre[2], surface, PHIS[i], COSPHIS[i], SINPHIS[i], NY, NX);
		for (j=0; j<NTHETAS; j++)
			{
			if (THETAS[j] < mintheta)
				{
				FPVISIBLE[ifacet*NOPS+i*NTHETAS+j]=0;
				}
			else
				{
				FPVISIBLE[ifacet*NOPS+i*NTHETAS+j]=1;
				}
			}
		
		}
	}

// BE WARNED!!! The definitions of x and y are swapped around here, and everything operates in integer units, with a correction only being made right at the end. Oh, the joy of importing a routine!!! Even trickier, cosphi is stil cos phi, and because we use a matrix with indices swapped, it ALL WORKS!!!
long double calc_min_theta(long double x, long double y, int x0, int x1, int y0, int y1, long double pz, float *surf, long double phi, long double cosphi, long double sinphi, int nx, int ny)
	{
	int dx, dy, temp;
	long double dxdist, dydist, delx, dely, xdist, ydist, thisz, dist;
	long double maxtantheta=-1000., tantheta, theta;
	
	dist=0.;
	
	while (phi > twopi) {phi -= twopi;}
	while (phi < 0.) {phi += twopi;}
	
	// this ensures that x1 and y1 are the 'leading' points, while x0,y0 are the lagging
	
	if (phi < piontwo && phi > -piontwo) {dy=1;} else {dy = -1; temp=y0; y0=y1; y1=temp;}
	if (phi < 0.) {dx=-1; temp=x0; x0=x1; x1=temp;} else {dx=1;} 
	

// note these dels and distances may be negative!
	if (cosphi==0) {dxdist=9e9;} else {dxdist=1./cosphi;}
	if (sinphi == 0) {dydist=9e9;} else {dydist=1./sinphi;}
	
	delx=x1-x;
	dely=y1-y;
	
	// calculates distance to next x-line and y-line
	xdist=delx/cosphi*dx;
	ydist=dely/sinphi*dy;
	
	while (y1 < ny && y1 >= 0 && x1 < nx && x1 >= 0)
		{
		
		if (fabs(ydist) < fabs(xdist))
			{
			thisz=get_height(x0,x1,surf[y1*nx+x0], surf[y1*nx+x1],x+ydist*cosphi*dy);
		//	tantheta=tan_theta(thisz, pz, ydist);
			tantheta=(thisz-pz)/fabs(ydist);
			if (tantheta > maxtantheta)
				{
				maxtantheta=tantheta;
				dist=fabs(ydist);
				}
			y1 += dy;
			y0 += dy;
			ydist += dydist;
			}
		else
			{
			thisz=get_height(y0,y1,surf[y0*nx+x1], surf[y1*nx+x1], y+xdist*sinphi*dx);
		//	tantheta=tan_theta(thisz, pz, xdist);
			tantheta=(thisz-pz)/fabs(ydist);
			if (tantheta > maxtantheta)
				{
				maxtantheta=tantheta;
				dist=fabs(xdist);
				}
			x1 += dx;
			x0 += dx;
			xdist += dxdist;
			}
		}
	theta=atan(maxtantheta/DD);
	return theta;
	
	}


/*
// for a given output direction, calculates the semi-major axes
// of the fresnel zone
void get_zone(track atrack, long double f_min, int thetai, int phij, long double *oc
	{
	long double cos_theta_i,
	
// first calculates offset of track to surface.
// this gives the central point of the ellipse

	cos_theta_t = OPBASISK[thetai][phij][2];
	sin_theta_t = pow(1.-cos_theta_t*cos_theta*t,0.5);
	sin_theta_i = sin_theta_t*REFNt/REFNi;
	
	magxy=pow(1.-OPBASISK[i][j][2]*OPBASISK[i][j][2],0.5);
	oc[0] = atrack.centre[2]/cos_theta_i*OPBASISK[thetai][phij][0]/magxy;
	oc[1] = atrack.centre[2]/cos_theta_i*OPBASISL[thetai][phij][1]/magxy;
	oc[2] = 0.;
	
// now we get the axes...
	// in perp direction, d_dist
	
	}
*/
/*
bool in_zone(track atrack, facet afacet, long double fmax)
	{
	bool OK;
	if (calc_dist(afacet.centre, atrack.centre) >  MAX_FT_DIST)
		{
		OK=0;
		}
	else
		{
		OK=1.;
		}
	return OK;
	}
*/

