#pragma once

// hysical constants
const long double c_light= 2.99792458e8l; // speed of light in vacuum, m/s
const long double  k_b=1.3806504e-23; // Boltzman constant
const long double Z_nought=376.73031; // impedence of free space


// lunar constants
long double R_m=1.7500e6; // Radius of the moon, m
const long double d_E_m=3.844e8; // Earth-moon distance,m
const long double d_E_m_2=d_E_m*d_E_m; // distance squared
long double moon_temp=225.; //temperature of the moon in kelvin
const long double Zmoon=14.295948; // from ZHS: this is the mean by weight I think
const long double Amoon=Zmoon*2.; // only approximate
//const long double Zmoon=11.;
//const long double Amoon=22.;

//conversion factors
const long double rad_to_deg = 180.l/M_PI;
const long double todeg = rad_to_deg, to_deg=todeg;
const long double deg_to_rad = M_PI/180.l;
const long double torad=deg_to_rad, to_rad=torad;

//short fractions
const long double twopi=2.l*M_PI;
const long double fourpi=4.l*M_PI;
const long double pion2=M_PI/2.l;
const long double piontwo=M_PI/2.l;
const long double one_third=1.l/3.l, two_thirds=2.l/3.l,four_thirds=4.l/3.l;
const long double threepiontwo=1.5l*M_PI;
const long double root_two=pow(2,0.5l);
const long double exp_e=2.71828182845904523536;
const long double sqrt_exp_e=pow(exp_e,0.5l);

// miscellaneous physical constants
const long double N_A=6.023e23; // Avogadro's number
const long double N_A_26=6.023e26; // Avogadro's number per kilogram
const long double alpha_fine=7.29735257e-3;
const long double Lambda_e=3.8616e-11; // compton wavelength of electron in cm
const long double M_mu=105.658369e6, M_e=510998.918, M_tau=1.777e9; // lepton masses in eV
const long double e_charge=1.60217646e-19; // C: charge of the electron in Coulombs
const long double epsilon_0=8.85418782e-12; //m^-3 kg^-1 s^2 C^2
