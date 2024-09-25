"""

The program performs the following funcitons:
    - Calculates a cascade profile for the regolith
        give a cosmic energy E
    - Divides that cascade up into N tracks of varying
        magnitudes and equal widths
    - Writes out N track files for a cosmic ray impact
        at some angle theta to the surface normal

Author: Clancy W. James

Uses scientific results from:
    ZHS: Zas, Halzen, Stanev, PhysRevD 45 (1992) 362
        "Electromagnetic pulses from high-energy showers:
            Implications for neutrino detection"
        https://ui.adsabs.harvard.edu/abs/1992PhRvD..45..362Z/abstract

    AMZ: Alvarez-Muniz & Zas, PhysLettB 434 (1998) 396
        "The LPM effect for EeV hadronic showers in ice:
            implications for radio detection of neutrinos"
        https://ui.adsabs.harvard.edu/abs/1998PhLB..434..396A/abstract
    
    AMMVZ: Alvarez-Muniz et al, PhysRevD 74 (2006) 023007 
        "Coherent radio pulses from showers in different media:
            A unified parametrization"
        https://ui.adsabs.harvard.edu/abs/2006PhRvD..74b3007A/abstract
"""

import numpy as np
from matplotlib import pyplot as plt

global rhoIce,rhoReg,chi0Reg,chi0Ice,RmollierReg,EcritIce,EcritReg,\
    RefNReg,RefNIce,ElpmReg,kEReg,kEIce

# constant medium properties
rhoIce = 0.92
rhoReg = 1.8 #g/cm2
chi0Reg = 22.59 #g/cm2
chi0Ice = 36.08
RmollierReg = 11.7 #g/cm3
EcritIce = 73.0
EcritReg = 40.0e6 #eV
RefNReg = 1.73
RefNIce = 1.78
ElpmReg = 770e12 # eV
kEReg = 3.3e-22 # constant of EM emission. 1e-6 less due to nu in Hz, not Mhz
kEIce = 4.79e-22



# other constants I should probably import from scipy
global clight,echarge,eps0
clight = 299792458
echarge = 1.602176634e-19
eps0 = 8.85418782E-12    

    
def main():
    """
    Main program:
    
    
    Required Arguments:
        E (float, eV): cosmic ray energy in eV (e.g. 1e18)
    
        N (int): Number of pieces in which to divide the cascade
    
        Theta (float, deg): IGNORED Angle of impact relative to the surface
            normal (0 is down). Always impacts in the x-z plane
        
    Secondary arguments:
    
    """
    from argparse import ArgumentParser

    a = ArgumentParser()
    a.add_argument('-E','--energy', type=float)
    a.add_argument('-N','--ntracks', type=int)
    a.add_argument('-T','--theta', type=float)
    a.add_argument('-O','--outprefix', type=str)
    args = a.parse_args()
    
    energy, ntracks, theta, prefix=check_inputs(args)
    
    Nmax,XmaxChi_m,X0,Lambda_m=getHillasParams(energy)
    
    # splits the tracks into N segements
    # each has ampitude, start, and stop positions
    xslims,ys = splitHillas(ntracks,Nmax,XmaxChi_m,X0,Lambda_m,plot=prefix)
    
    write_files(ntracks,xslims,ys,theta,prefix)

def write_files(N,lslims,mags,theta,outname = "trackfile"):
    """
    ls here is the longitudinal direction
    x,y,z refers to simulation Cartesian coordinates
    """
    global clight
    
    theta_rad = theta * np.pi/180.
    ctheta = np.cos(theta_rad)
    stheta = np.sin(theta_rad)
    
    allfile=outname + "_alltracks.dat"
    allout = open(allfile,'w')
    allout.write(str(N)+"\n")
    
    print("lslims is ",lslims)
    for i in np.arange(N):
        l0=lslims[i] # m
        l1=lslims[i+1] #m
        mag=mags[i] # this is in units of fundamental charge
        #print(l0,l1,mag)
        
        x0 = 0.
        x1 = 0.
        y0=0.
        y1=0.
        
        z0 = l0
        z1 = l1
        
        t0 = l0/clight
        t1 = l1/clight
        
        beta=1.
        
        # writes output file
        outfile = outname + "_" + str(i) + ".dat"
        print("i value = ", i)
        out = open(outfile, 'w')
        out.write("1\n") # single track
        string="{0:1.5e} {1:6.3f} {2:6.3f} {3:6.3f} {4:1.3f} {5:5.3f} {6:5.3f} {7:5.3f} {8:1.5e}\n".format(t0,x0,y0,z0,beta,x1,y1,z1,mag)
        #print(string)
        out.write(string)
        out.close()
        
        allout.write(string)
    allout.close()

def splitHillas(Nsplit,Nmax,XmaxChi_m,X0,Lambda_m,thresh=1e-3,plot=None):
    """
    Splits up the GaisserHillas progile into N segments
    
    Thresh is threshold below which we don't care
    """
    intmin=0.
    intmax=XmaxChi_m*5.
    xs = np.linspace(intmin,intmax,1000)
    ys = GaiserHillas(xs,Nmax,XmaxChi_m,X0,Lambda_m)
    
    themax = np.max(ys)
    Nthresh = themax * thresh
    OK = np.where(ys > Nthresh)[0]
    newxs = xs[OK]
    newys = ys[OK]
    
    # now calculates the function at the midpoint of each segment
    # this is a gross simpification!
    xmin = newxs[0]
    xmax = newxs[-1]
    dx = (xmax-xmin)/Nsplit
    splitxs = np.linspace(xmin + dx/2.,xmax - dx/2.,Nsplit)
    splitxs_limits = np.linspace(xmin,xmax,Nsplit+1)
    splitys = GaiserHillas(splitxs,Nmax,XmaxChi_m,X0,Lambda_m)
    
    if plot is not None:
        xs = np.linspace(intmin,intmax,100)
        ys = GaiserHillas(xs,Nmax,XmaxChi_m,X0,Lambda_m)
        plt.figure()
        plt.plot(xs,ys,label='cascade')
        #plt.plot(splitxs,splitys)
        
        # uses dumb routines to make this plot
        plt.hist(splitxs,bins=splitxs_limits,weights=splitys,alpha=0.5,label='track distribution')
        
        ymax=np.max(ys)
        ymax=10**(int(np.log10(ymax))+1)
        ymin=ymax/1e6
        plt.ylim(ymin,ymax)
        
        plt.yscale('log')
        plt.xlabel('x [m]')
        plt.ylabel('N(x)')
        
        plt.legend()
        plt.tight_layout()
        plt.savefig(plot+'_cascade.pdf')
        plt.close()
    
    return splitxs_limits,splitys
    

def get_excess_tracklength(EeV):
    """
    Implements Eqs 1-3 from AMZ.
    Gives total excess tracklength
    
    Note: a higher density causes pi0 to interact before decaying,
        thus reducing EM content; but it also causes pi+- to
        do so, increasing EM content (produces more pi0, less mu)
        Hence, no modification here is made for density.
        The only scaling is by radiation length.
        
    Input: Eev (float): energy in eV
    
    Output: tptReg: total projected excess tracklength in m
    """
    global rhoIce,rhoReg,chi0Reg,chi0Ice,RmollierReg,EcritIce,EcritReg,\
        RefNReg,RefNIce,ElpmReg,kEReg,kEIce
    
    # epsilon factor: eq 2
    eps = np.log10(EeV/1e12)
    eps3 = eps+3 # why not defined at 1 GeV?
    
    feps = -1.27e-2 - 4.76e-2 * eps3 - 2.07e-3*eps3**2 + 0.52*eps3**0.5
    
    # total projected tracklength in m
    tIce = 6.25 * feps * (EeV/1e9)
    
    # total projected tracklength
    tptIce = 0.21*tIce
    
    # now scales to regolith, according to Unified Parameterisation Eq 11
    
    # since kE defined at the Cherenkov angle, we must account for this also
    sinthetaCIce = (1.-RefNIce**-2)**0.5
    sinthetaCReg = (1.-RefNReg**-2)**0.5
    tptReg = tptIce * (chi0Reg/chi0Ice) * (rhoIce/rhoReg) * \
        (kEReg/kEIce) * (sinthetaCReg/sinthetaCIce)
    
    return tptReg
    
    

def getHillasParams(energy,plot=None):
    """
    Gets parameters of a Gaiser-Hillas function for a cosmic ray
    cascade profile.
    
    We require:
        Cascade length: Based on Charlie Morley-Wong's fits to results from  
            J. Alvarez-Muniz, E. Zas, PHYs.Lett.B 434 (1998) 396
            Scaled to regolith according to density
        
        Total excess tracklength: Uses the results of 
            Alvarez-Muniz et al, PhysRevD 74 (2006) 023007,
            (which gives frequency spectrum), to derive total
            excess tracklength, then scaled to hadronic cascades
            using the "electromagnetic factor" of ...
        
    Scaled to the regolith by C.W. James
    """
    global rhoIce,rhoReg,chi0Reg,chi0Ice,RmollierReg,EcritIce,EcritReg,\
        RefNReg,RefNIce,ElpmReg,kEReg,kEIce
    
    log10E = np.log10(energy)
    
    # Xmax and lambda in units of the radiation length
    
    Lambda = -0.125*log10E + 3.177
    X0=0. # always start at 0
    
    # 14 rad lengths at 10 TeV
    # 18 at 100 PeV
    # (this parameter wasn;t fully fitted)
    XmaxChi = 1. + log10E # approximate
    
    # we take the overall length as due to the hadronic shower core
    # Hence, we first get actual length parameters in cm,
    # then normalise by density (proportional to nucleon density)
    
    ScaleIce = chi0Ice/rhoIce # cascade in cm lengths
    
    # now we convert to regolith
    ScaleReg = chi0Ice/rhoReg
    
    # performs the scaling to regolith units of m
    Lambda_m = Lambda * ScaleReg/100.
    XmaxChi_m = XmaxChi * ScaleReg/100.
    
    
    # we require Nmax to be fitted according to the integral of the profile
    # it must equal the total excess tracklength for the regolith
    from scipy.integrate import quad
    intmin=0.
    intmax=XmaxChi_m*5.
    Nmax=1. # leading constant - peak Nparticles
    X0=0. # starting point of cascade
    
    # total excess tracklength in m for normalisation of unity
    C = quad(GaiserHillas,intmin,intmax, args=(Nmax,XmaxChi_m,X0,Lambda_m))
    C = C[0]
    
    # we require a constant for the total tracklength
    req_norm = get_excess_tracklength(energy)
    Nmax = req_norm/C
    
    if plot is not None:
        xs = np.linspace(intmin,intmax,100)
        ys = GaiserHillas(xs,Nmax,XmaxChi_m,X0,Lambda_m)
        
        ymax=np.max(ys)
        ymax=10**(int(np.log10(ymax))+1)
        ymin=ymax/1e6
        plt.ylim(ymin,ymax)
        
        plt.figure()
        plt.plot(xs,ys)
        plt.yscale('log')
        plt.xlabel('x [m]')
        plt.ylabel('N(x)')
        plt.tight_layout()
        plt.savefig(plot+'_cascade.pdf')
        plt.close()
    
    return Nmax,XmaxChi_m,X0,Lambda_m
    
def failed_get_excess_tracklength(EeV):
    """
    Calculates normalisation for total excess tracklength.
    
    This routine is incorrect and I don't know why...
    
    Ensures that a total excess tracklength matches the electric field
    
    We get this by solving Eq 10 from AMZ:
        rE = kE * (E/Ec) * (chi0/rho) * nu * np.sin(theta)
    
        and ZHS Eq 13: normalisation from ZHS formula;
        rE = echarge * mur * omega / (2*pi*eps0*c^2) * vperp * t
        
        hence vt = kE (E/Ec)* (chi0/rho) * nu /
            echarge * mur * omega / (2*pi*eps0*c^2)
        hence vt = kE (E/Ec)* (chi0/rho) *(eps0*c^2) /
            echarge)
            
        Possible 2 pi difference due to FT definitions
            
    Returns: vt in units of total excess tracklength
        (units: fundamental charge * m)
    
    """
    global rhoReg, chi0Reg,RmollierReg,EcritReg,RefNReg,ElpmReg,RefNReg
    global eps0, kEReg, clight, echarge
    
    vt_cm = kEReg * (EeV/EcritReg) * (chi0Reg/rhoReg) * eps0\
        * clight**2 / echarge
    vt_m = vt_cm/100 # conversion from cm to m
    
    print("WARNING: the total normalisation needs to be revised")
    print("Currently re-divided by 10^4")
    print("Please consult with Jaime...")
    vt_m /= 1e4
    
    return vt_m 

def GaiserHillas(X,Nmax,Xmax,X0,Lambda):
    """
    Gaisser-Hillas function:
        - X (float): cascade length (grams/cm3)
        - Xmax (float): cascade peak (g/cm3)
        - X0 (float): shower start (g/cm3)
        - Lambda (float): shape parameter
    
    Returns:
        N(x): number of particles as function of cascade length
    
    """
   
    N = Nmax * ((X-X0)/(Xmax-X0))**((Xmax-X0)/Lambda)\
        * np.exp((Xmax-X/Lambda))
    return N

def check_inputs(a):
    """
    a is the argparse objects
    """
    
    err=False
    
    # checks inputs
    if not a.energy:
        print("No energy provided. Please use -E [energy] in eV")
        err=True
    elif a.energy <= 0:
        print("Energy ",a.energy," too small - must be positive")
        err=True
    else:
        energy=a.energy
    
    if not a.ntracks:
        print("No ntracks provided. Please use -N [ntracks]")
        err=True
    elif a.ntracks <=0:
        print("Ntracks value ",a.ntracks," negative: must be positive integer")
        err=True
    else:
        ntracks = a.ntracks
    
    
    if not a.theta:
        print("Theta not provided, Please use -T [theta]")
        err=True
    if a.theta > 90:
        print("Theta value ",a.theta," too big: 0 =< theta =< 90")
        err=True
    elif a.theta < 0.:
        print("Theta value ",a.theta," too small: 0 =< theta =< 90")
        err=True
    else:
        theta=a.theta
    
    if not a.outprefix:
        print("prefix set to 'example'")
        prefix="example"
    else:
        prefix=a.outprefix
    
    if err:
        print("Errors occurred, exiting...")
        exit()
    return energy,ntracks,theta,prefix
    
main()
