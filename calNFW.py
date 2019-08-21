'''
This function calculates the enclosed mass of an NFW profile, based on its redshift,
reffp, and meffp, where meffp is the enclosed mass (in Mssun) of the DM halo at reffp
(in kpc).
The virial over-density and virial radius are defined in Bryan G. L., Norman M. L., 1998
The concentration mass relation can be found in Aaron A. Dutton  Andrea V. Macciò, 2014
(Eqs. 12 & 13)
'''
import numpy as np


def calNFW(rlist,reffp,meffp,zr=0): 
    rhocrit=127*(0.27*(1+zr)*(1+zr)*(1+zr)+0.73)
    xLam=-0.73/(0.27*(1+zr)*(1+zr)*(1+zr)+0.73)
    delta=18.0*np.pi*np.pi+82.0*(xLam)-39*(xLam)*(xLam)
    samlogMv=np.arange(7.0,14.0,0.05)
    logc=0.537+(1.025-0.537)*np.exp(-0.718*np.power(zr,1.08))+(0.024*zr-0.097)*(samlogMv-np.log10(0.7)-np.log10(1e12))
    cont=np.power(10,logc)
    Rvirn=np.power(np.power(10,samlogMv)*3.0/4.0/np.pi/rhocrit/delta,1.0/3.0)
#       finding the NFW profile (at z=0) that fits at reff and plot it
    #x=reffp/Rvirn
    #rho = rhocrit*delta/cont/cont/cont/x/(1/cont+x)/(1/cont+x)
    Rs = Rvirn/cont
    cdelta = cont*cont*cont*delta/3./(np.log(1+cont)-cont/(1+cont))
    #enclosed NFW halo mass within reffp with different total halo masses
    menc = 4*np.pi*rhocrit*cdelta*Rs*Rs*Rs*(np.log((Rs+reffp)/Rs)-reffp/(Rs+reffp))
    #infer the total halo mass with meffp
    logMveff = np.interp(meffp, menc, samlogMv)
    #calculate the enclosed mass with the total halo mass above
    logceff=0.537+(1.025-0.537)*np.exp(-0.718*np.power(zr,1.08))+(0.024*zr-0.097)*(logMveff-np.log10(0.7)-np.log10(1e12))
    conteff=np.power(10,logceff)
    Rvirneff=np.power(np.power(10,logMveff)*3.0/4.0/np.pi/rhocrit/delta,1.0/3.0)
    Rseff = Rvirneff/conteff
    cdeltaeff = conteff*conteff*conteff*delta/3./(np.log(1+conteff)-conteff/(1+conteff))
    #output the enclosed NFW halo mass in rlist that matches [reffp,meffp]
    mnfw = np.array(4*np.pi*rhocrit*cdeltaeff*Rseff*Rseff*Rseff*(np.log((Rseff+rlist)/Rseff)-rlist/(Rseff+rlist)))
    return {'mnfw':mnfw}


def main():
    rlist = np.linspace(0.1,20,num=20)
    reffp = 4.7; meffp = 7.0e9;
    info = calNFW(rlist,reffp,meffp,zr=0)
    print 'mnfw', info['mnfw']
    return None


if __name__== "__main__":
      main()
