
// physical constants
const double c_light= 2.99792458e8; // speed of light in vacuum, m/s
const double  k_b=1.38066e-23; // Boltzman constant
const double Z_nought=119.9169832*M_PI; // impedence of free space
const double h_planck = 6.626068e-34;


// lunar constants
double R_m=1.7500e6; // Radius of the moon, m
const double d_E_m=3.844e8; // Earth-moon distance,m
const double d_E_m_2=d_E_m*d_E_m; // distance squared
double moon_temp=225.; //temperature of the moon in kelvin
const double Zmoon=14.295948; // from ZHS: this is the mean by weight I think
const double Amoon=Zmoon*2.; // only approximate
//const double Zmoon=11.;
//const double Amoon=22.;

//conversion factors
const double rad_to_deg = 180./M_PI;
const double todeg = rad_to_deg, to_deg=todeg;
const double deg_to_rad = M_PI/180;
const double torad=deg_to_rad, to_rad=torad;

//short fractions
const double twopi=2.*M_PI;
const double fourpi=4.*M_PI;
const double pion2=M_PI/2.;
const double piontwo=M_PI/2.;
const double one_third=1./3., two_thirds=2./3.,four_thirds=4./3.;
const double threepiontwo=1.5*M_PI;
const double root_two=pow(2,0.5);
const double exp_e=2.7182818284;
const double sqrt_exp_e=pow(exp_e,0.5);

// miscellaneous physical constants
const double N_A=6.023e23; // Avogadro's number
const double N_A_26=6.023e26; // Avogadro's number per kilogram
const double alpha_fine=7.29735257e-3;
const double Lambda_e=3.8616e-11; // compton wavelength of electron in cm
const double M_mu=105.658369e6, M_e=510998.918, M_tau=1.777e9; // lepton masses in eV
const double e_charge=1.60217646e-19; // C: charge of the electron in Coulombs
const double epsilon_0=8.85418782e-12; //m^-3 kg^-1 s^2 C^2
const double mu_0=fourpi*1e-7; // permeability free space
