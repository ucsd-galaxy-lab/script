'''
This is a function to read AHF catalog and halo files.
The path to the catalog/file directory is dir = homedir+maindir+rundir
(set them to fit your need)

Catalog should be in the Pep folder; halo file should be in the halos folder.


comoving =1: then the function will output COMOVING distance (not physical)
if set hubble= (little h; to be read from the snapshot header): the code will remove the h in the unit;
otherwise (e.g. in idealized runs), set hubble=1

read_halo_history_pep: read AHF catalog
should put: snapshottime_FIRE1.txt and snapshot_times_FIRE2.txt into the same folder
(to generate catalog file name)


'lambdalist': spin parameter defined in Bullock+01
'redshift': redshift
'ID': halo ID
'x','y','z': halo coordinate
'xv','yv','zv': halo velocity
'M': halo virial mass
'Ms': stellar mass
'Mg': gas mass
'R': virial radius
'fMhires': percentage of high resolution particles (1=no pollution from low resolution region)
'Lxstar','Lystar','Lzstar': angular momentum of stars
'Lxgas','Lygas','Lzgas': angular momentum of gas

'''
import numpy as np


def snapshotzstr(firever=1):
        if firever==1:
                spmname='snapshottime_FIRE1.txt'
                spmfile=open(spmname,"r")
                dars = spmfile.readlines()
                snapno=[]
                redshift=[]
                for line in dars:
                        xsd = line.split()
                        snapno = np.append(snapno, int(xsd[0]))
                        redshift = np.append(redshift, float(xsd[2]))
                strzlist=np.array(redshift)
                strnolist=np.array(snapno)
        elif firever==2:
                spmname='snapshot_times_FIRE2.txt'
                spmfile=open(spmname,"r")
                spmfile.readline()
                spmfile.readline()
                spmfile.readline()
                dars = spmfile.readlines()

                snapno=[]
                redshift=[]

                for line in dars:
                        xsd = line.split()
                        snapno = np.append(snapno, int(xsd[0]))
                        redshift = np.append(redshift, float(xsd[2]))
                spmfile.close()
                snapno=np.array(snapno)
                redshift=np.array(redshift)
                strzlist=[]
                strnolist=[]
                for icount in range(len(snapno)):
                        if snapno[icount]>-1:
                                strno=str(int(snapno[icount]))
                                if icount <10:
                                        strno =  '00' + strno
                                elif icount<100:
                                        strno = '0'+strno
                                if strno=='250':
                                        strz = '1.187'
                                elif strno=='000':
                                        strz = '99.013'
                                elif strno=='024':
                                        strz = '9.313'
                                else:
                                        strz="%0.3f" % redshift[icount]
                        strzlist = np.append(strzlist,strz)
                        strnolist = np.append(strnolist,strno)
        return strnolist, strzlist


def read_halo_history(rundir, halonostr='00', multifile='n',suffixadd='', hubble=1,comoving=1,\
                      homedir='/home/tkc004/', maindir='scratch', singlesnap=0, atime=1, snumadd=0): #to output comoving distance comoving = 1
        redlist = []
        halolist = []
        xlist = []
        ylist = []
        zlist = []
        xvlist = []
        yvlist = []
        zvlist = []
        mvirlist = []
        mstarlist = []
        mgaslist = []
        rvirlist = []
        fMhireslist = []
        lambdalist = []
        Lxgaslist = []
        Lygaslist = []
        Lzgaslist = []
        Lxstarlist = []
        Lystarlist = []
        Lzstarlist = []
        dirname=homedir+maindir+'/'+rundir
        halofile=open(dirname+'/halos/halo_000'+halonostr+suffixadd+'.dat','r')
        halofile.readline()
        halofile.readline()
        dars = halofile.readlines()
        halofile.close()
        mvirno = 4
        xcenno = 6
        ycenno = 7
        zcenno = 8
        xvcenno = 9
        yvcenno = 10
        zvcenno = 11
        zredno = 0
        halono = 1
        mstarno = 74
        rvirno = 12
        fMhiresno = 38
        mgasno = 54
        lambdano=20
        Lxgasno = 57
        Lygasno = 58
        Lzgasno = 59
        Lxstarno = 77
        Lystarno = 78
        Lzstarno = 79
        if snumadd==1:
                mvirno+= 1
                xcenno+= 1
                ycenno+= 1
                zcenno+= 1
                xvcenno+= 1
                yvcenno+= 1
                zvcenno+= 1
                zredno+= 1
                halono+= 1
                mstarno+= 1
                rvirno+= 1
                fMhiresno+= 1
                mgasno+= 1
                lambdano+= 1
                Lxgasno+= 1
                Lygasno+= 1
                Lzgasno+= 1
                Lxstarno+= 1
                Lystarno+= 1
                Lzstarno+= 1

        for line in dars:
                xsd = line.split()
                mvir = float(xsd[mvirno])
                xcen = float(xsd[xcenno])
                ycen = float(xsd[ycenno])
                zcen = float(xsd[zcenno])
                xvcen = float(xsd[xvcenno])
                yvcen = float(xsd[yvcenno])
                zvcen = float(xsd[zvcenno])
                Lxgas = float(xsd[Lxgasno])
                Lygas = float(xsd[Lygasno])
                Lzgas = float(xsd[Lzgasno])
                Lxstar = float(xsd[Lxstarno])
                Lystar = float(xsd[Lystarno])
                Lzstar = float(xsd[Lzstarno])
                if comoving==1:
                        zred=0.0
                else:
                        zred=float(xsd[zredno])
                redlist=np.append(redlist,float(xsd[zredno]))
                halolist=np.append(halolist,int(xsd[halono]))
                mvirlist=np.append(mvirlist,mvir/hubble)
                mstarlist=np.append(mstarlist,float(xsd[mstarno])/hubble)
                xlist=np.append(xlist,xcen/hubble/(1.0+zred)) # xcen originally in comoving unit (kpc/h)
                ylist=np.append(ylist,ycen/hubble/(1.0+zred))
                zlist=np.append(zlist,zcen/hubble/(1.0+zred))
                xvlist=np.append(xvlist,xvcen) #originally in pecular velocity
                yvlist=np.append(yvlist,yvcen)
                zvlist=np.append(zvlist,zvcen)
                rvirlist=np.append(rvirlist,float(xsd[rvirno])/hubble/(1.0+zred)) # originally in comoving unit (kpc/h)
                fMhireslist= np.append(fMhireslist, float(xsd[fMhiresno]))
                mgaslist = np.append(mgaslist, float(xsd[mgasno])/hubble)
                lambdalist = np.append(lambdalist, float(xsd[lambdano])) #spin parameter defined in Bullock+2001
                Lxgaslist = np.append(Lxgaslist, float(xsd[Lxgasno]))
                Lygaslist = np.append(Lygaslist, float(xsd[Lygasno]))
                Lzgaslist = np.append(Lzgaslist, float(xsd[Lzgasno]))
                Lxstarlist = np.append(Lxstarlist, float(xsd[Lxstarno]))
                Lystarlist = np.append(Lystarlist, float(xsd[Lystarno]))
                Lzstarlist = np.append(Lzstarlist, float(xsd[Lzstarno]))
        redlist=np.array(redlist)
        halolist=np.array(halolist)
        mvirlist=np.array(mvirlist)
        mstarlist=np.array(mstarlist)
        mgaslist=np.array(mgaslist)
        xlist=np.array(xlist)
        ylist=np.array(ylist)
        zlist=np.array(zlist)
        xvlist=np.array(xvlist)
        yvlist=np.array(yvlist)
        zvlist=np.array(zvlist)
        rvirlist=np.array(rvirlist)
        fMhireslist=np.array(fMhireslist)
        lambdalist=np.array(lambdalist)
        Lxgaslist=np.array(Lxgaslist)
        Lygaslist=np.array(Lygaslist)
        Lzgaslist=np.array(Lzgaslist)
        Lxstarlist=np.array(Lxstarlist)
        Lystarlist=np.array(Lystarlist)
        Lzstarlist=np.array(Lzstarlist)
        if snumadd==1:
                redlist=redlist[::-1]
                halolist=halolist[::-1]
                mvirlist=mvirlist[::-1]
                mstarlist=mstarlist[::-1]
                mgaslist=mgaslist[::-1]
                xlist=xlist[::-1]
                ylist=ylist[::-1]
                zlist=zlist[::-1]
                xvlist=xvlist[::-1]
                yvlist=yvlist[::-1]
                zvlist=zvlist[::-1]
                rvirlist=rvirlist[::-1]
                fMhireslist=fMhireslist[::-1]
                lambdalist=lambdalist[::-1]
                Lxgaslist=Lxgaslist[::-1]
                Lygaslist=Lygaslist[::-1]
                Lzgaslist=Lzgaslist[::-1]
                Lxstarlist=Lxstarlist[::-1]
                Lystarlist=Lystarlist[::-1]
                Lzstarlist=Lzstarlist[::-1]
        if singlesnap==0:
                        return {'lambdalist':lambdalist, 'redshift':redlist, 'ID':halolist, 'x':xlist, 'y':ylist, 'z':zlist,\
                        'xv':xvlist, 'yv':yvlist, 'zv':zvlist, 'M':mvirlist, 'Ms':mstarlist, 'Mg':mgaslist, 'R':rvirlist,\
                        'fMhires':fMhireslist, 'Lxstar':Lxstarlist, 'Lystar':Lystarlist, 'Lzstar':Lzstarlist,\
                        'Lxgas':Lxgaslist, 'Lygas':Lygaslist, 'Lzgas':Lzgaslist}
        else:
                aAHF = 1.0/(1.0+redlist)
                rednow = 1.0/atime-1.0
                print 'atime', atime
                xcens = np.interp(atime,aAHF,xlist)
                ycens = np.interp(atime,aAHF,ylist)
                zcens = np.interp(atime,aAHF,zlist)
                xvcens = np.interp(atime,aAHF,xvlist)
                yvcens = np.interp(atime,aAHF,yvlist)
                zvcens = np.interp(atime,aAHF,zvlist)
                MAHF = np.interp(atime,aAHF,mvirlist)
                SmAHF = np.interp(atime,aAHF,mstarlist)
                GmAHF = np.interp(atime,aAHF,mgaslist)
                RAHF = np.interp(atime,aAHF,rvirlist)
                fAHF = np.interp(atime,aAHF,fMhireslist)
                lambdaspin = np.interp(atime,aAHF,lambdalist)
                ids = np.interp(atime,aAHF,halolist)
                Lxgas = np.interp(atime,aAHF,Lxgaslist)
                Lygas = np.interp(atime,aAHF,Lygaslist)
                Lzgas = np.interp(atime,aAHF,Lzgaslist)
                Lxstars = np.interp(atime,aAHF,Lxstarlist)
                Lystars = np.interp(atime,aAHF,Lystarlist)
                Lzstars = np.interp(atime,aAHF,Lzstarlist)
                return {'lambdalist':lambdaspin, 'redshift':rednow, 'ID':ids, 'x':xcens, 'y':ycens, 'z':zcens,\
                'xv':xvcens, 'yv':yvcens, 'zv':zvcens, 'M':MAHF, 'Ms':SmAHF, 'Mg':GmAHF, 'R':RAHF, 'fMhires':fAHF,\
                'Lxstar':Lxstars,'Lystar':Lystars,'Lzstar':Lzstars,'Lxgas':Lxgas,'Lygas':Lygas,'Lzgas':Lzgas}


def read_halo_history_pep(rundir, finalno, Sheaform=0, snapsep=1, singlesnap=1, beginno=100,halonostr='00',\
                          multifile='n', hubble=1,comoving=1, maindir='scratch',firever=2, outputallhalos=0,\
                          homedir='/home/tkc004/'): #to output comoving distance comoving = 1): #to output comoving distance comoving = 1
        redlist = []
        idlist = []
        xlist = []
        ylist = []
        zlist = []
        xvlist = []
        yvlist = []
        zvlist = []
        Mlist = []
        Mslist = []
        Mglist = []
        Rvlist = []
        fMlist = []
        lambdalist = []
        Lxgaslist=[]
        Lygaslist=[]
        Lzgaslist=[]
        Lxstarlist=[]
        Lystarlist=[]
        Lzstarlist=[]
        halono=int(halonostr)
        strnolist, zstrlist = snapshotzstr(firever=firever)
        if (singlesnap==1): beginno=finalno
        print 'beginno', beginno
        for Nsnap in range(beginno,finalno+1,snapsep):
                redshiftstring = str(zstrlist[Nsnap])
                #Nsnapstring = str(strnolist[Nsnap]) 
                Nsnapstring = str(Nsnap)
                hcat=read_halo_catalog(rundir, Nsnapstring, redshiftstring, usembp='n',\
                                       maindir=maindir,Sheaform=Sheaform,homedir=homedir)
                idl = hcat['ID']
                mvirl = hcat['M']
                xl = hcat['x']
                yl = hcat['y']
                zl = hcat['z']
                xvl = hcat['vx']
                yvl = hcat['vy']
                zvl = hcat['vz']
                Msl = hcat['Mstar']
                Mgl = hcat['Mgas']
                Rl = hcat['Rvir']
                fMl = hcat['fMhires']
                lambdal =hcat['lambdano']
                Lxgasl = hcat['Lxgas']
                Lygasl = hcat['Lygas']
                Lzgasl = hcat['Lzgas']
                Lxstarl = hcat['Lxstar']
                Lystarl = hcat['Lystar']
                Lzstarl = hcat['Lzstar']
                redshift = float(redshiftstring)
                if comoving==1:
                        zred=0.0
                else:
                        zred=redshift
                haloid = halono
                mvirl = mvirl/hubble
                Msl = Msl/hubble
                Mgl = Mgl/hubble
                xl = xl/hubble/(1.0+zred)
                yl = yl/hubble/(1.0+zred)
                zl = zl/hubble/(1.0+zred)
                Rl = Rl/hubble/(1.0+zred)
                if outputallhalos==1 and singlesnap==1:
                        print 'output all halos'
                        return {'fMhires':fMl, 'redshift':redshift, 'ID':idl, 'x':xl, 'y':yl, 'z':zl,\
                        'xv':xvl, 'yv':yvl, 'zv':zvl, 'M':mvirl, 'Ms':Msl, 'Mg':Mgl, 'R':Rl,\
                        'Lxgas':Lxgasl, 'Lygas':Lygasl, 'Lzgas':Lzgasl,\
                        'lambdalist':lambdal,'Lxstar':Lxstarl, 'Lystar':Lystarl, 'Lzstar':Lzstarl}
                mvir = float(mvirl[halono])
                mstar = float(Msl[halono])
                mgas = float(Mgl[halono])
                xcen = float(xl[halono]) # xcen originally in comoving unit (kpc/h)
                ycen = float(yl[halono])
                zcen = float(zl[halono])
                xvcen = float(xvl[halono])
                yvcen = float(yvl[halono])
                zvcen = float(zvl[halono])
                rvir = float(Rl[halono])  # originally in comoving unit (kpc/h)
                lambdaB = float(lambdal[halono])
                Lxgas = float(Lxgasl[halono])
                Lygas = float(Lygasl[halono])
                Lzgas = float(Lzgasl[halono])
                Lxstar = float(Lxstarl[halono])
                Lystar = float(Lystarl[halono])
                Lzstar = float(Lzstarl[halono])
                xlist = np.append(xlist,xcen)
                ylist = np.append(ylist,ycen)
                zlist = np.append(zlist,zcen)
                xvlist = np.append(xvlist,xvcen)
                yvlist = np.append(yvlist,yvcen)
                zvlist = np.append(zvlist,zvcen)
                idlist = np.append(idlist,haloid)
                redlist = np.append(redlist,float(zstrlist[Nsnap]))
                Mlist = np.append(Mlist,mvir)
                Mslist = np.append(Mslist,mstar)
                Mglist = np.append(Mglist,mgas)
                Rvlist = np.append(Rvlist,rvir)
                fMlist = np.append(fMlist,fMl)
                lambdalist = np.append(lambdalist,lambdaB)
                Lxgaslist = np.append(Lxgaslist,Lxgas)
                Lygaslist = np.append(Lygaslist,Lygas)
                Lzgaslist = np.append(Lzgaslist,Lzgas)
                Lxstarlist = np.append(Lxstarlist,Lxstar)
                Lystarlist = np.append(Lystarlist,Lystar)
                Lzstarlist = np.append(Lzstarlist,Lzstar)
        xlist = np.array(xlist)
        ylist = np.array(ylist)
        zlist = np.array(zlist)
        xvlist = np.array(xvlist)
        yvlist = np.array(yvlist)
        zvlist = np.array(zvlist)
        idlist = np.array(idlist)
        redlist = np.array(redlist)
        Mlist = np.array(Mlist)
        Mslist = np.array(Mslist)
        Mglist = np.array(Mglist)
        Rvlist = np.array(Rvlist)
        fMlist = np.array(fMlist)
        Lxgaslist = np.array(Lxgaslist)
        Lygaslist = np.array(Lygaslist)
        Lzgaslist = np.array(Lzgaslist)
        Lxstarlist = np.array(Lxstarlist)
        Lystarlist = np.array(Lystarlist)
        Lzstarlist = np.array(Lzstarlist)
        lambdalist = np.array(lambdalist)
        if (singlesnap==1):
                return {'lambdalist':lambdaB, 'fMhires':fMl, 'redshift':redshift, 'ID':haloid,\
                'x':xcen, 'y':ycen, 'z':zcen, 'xv':xvcen, 'yv':yvcen, 'zv':zvcen,\
                'M':mvir, 'Ms':mstar, 'Mg':mgas, 'R':rvir, 'Lxstar':Lxstar, 'Lystar':Lystar,'Lzstar':Lzstar,\
                'Lxgas':Lxgas, 'Lygas':Lygas,'Lzgas':Lzgas}
        else:
                return {'lambdalist':lambdalist, 'fMhires':fMlist, 'redshift':redlist, 'ID':idlist,\
                'x':xlist, 'y':ylist, 'z':zlist, 'xv':xvlist, 'yv':yvlist, 'zv':zvlist,\
                'M':Mlist, 'Ms':Mslist, 'Mg':Mglist, 'R':Rvlist, 'Lxstar':Lxstarlist, 'Lystar':Lystarlist,'Lzstar':Lzstarlist,\
                'Lxgas':Lxgaslist, 'Lygas':Lygaslist,'Lzgas':Lzgaslist}


def read_halo_catalog(rundir, Nsnapstring, redshiftstring, usembp='n',maindir='scratch',Sheaform=0,homedir='/home/tkc004/'):
        finname = homedir+maindir+'/'+rundir+'/Pep/snap'+ Nsnapstring+ 'RPep.z'+redshiftstring+'.AHF_halos' #or whatever prefix
        f = open(finname)
        C = np.loadtxt(f)
        halID = C[:,0]
        x = C[:,5]
        y = C[:,6]
        z = C[:,7]
        vx = C[:,8]
        vy = C[:,9]
        vz = C[:,10]
        M = C[:,3]
        Vmax = C[:,16]
        Vsig = C[:,18]
        lambdano = C[:,19] #Bullock+01 definition
        Rvir = C[:,11]
        Lx = C[:,21]
        Ly = C[:,22]
        Lz = C[:,23]
        Mgas = C[:,53]
        if Sheaform==1:
                Mstar = C[:,64]
        else:
                Mstar = C[:,73]
        fMhires = C[:,37]
        Lxgas = C[:,56]
        Lygas = C[:,57]
        Lzgas = C[:,58]
        Lxstar = C[:,76]
        Lystar = C[:,77]
        Lzstar = C[:,78]
        #if (use_mbp=='y'): not even going to bother writing this for now - it isnt very accurate for halo center in Amiga. May want to use COM instead though

        return {'Lxgas':Lxgas, 'Lygas':Lygas, 'Lzgas':Lzgas,'Lxstar':Lxstar, 'Lystar':Lystar, 'Lzstar':Lzstar, 'Lx':Lx,'Ly':Ly,'Lz':Lz,'fMhires':fMhires,'ID':halID, 'x':x, 'y':y, 'z':z,\
 'vx':vx, 'vy':vy, 'vz':vz, 'M':M, 'Vmax':Vmax, 'Vsig':Vsig, 'lambdano':lambdano,\
'Rvir':Rvir, 'Mgas':Mgas, 'Mstar':Mstar}


def main():
    rundir = 'FIRE_2_0_or_h61'
    #info = read_halo_history(rundir, halonostr='00', multifile='n',suffixadd='', hubble=0.702,comoving=0,\
    #                  homedir='/home/tkc004/', maindir='oasis', snumadd=0)
    finalno = 600 #the snapshot you want to read
    firever = 2 #FIRE version: 1-> FIRE 1; 2-> FIRE 2;
    info = read_halo_history_pep(rundir, finalno, Sheaform=0, snapsep=1, singlesnap=1, beginno=100,halonostr='00',\
                          multifile='n', hubble=0.702,comoving=0, maindir='oasis',firever=firever, outputallhalos=0,\
                         homedir='/home/tkc004/')
    print 'virial mass', info['M']
    return None


if __name__== "__main__":
      main()
