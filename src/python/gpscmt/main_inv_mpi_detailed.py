#this is the main inversion file
#output detailed inversion result(e.g. dhat)
#load all GFs from *.npy (from zrt2enz.py)
from __future__ import division #for python2
import numpy as np
import matplotlib.pyplot as plt
import glob
import math
import datetime
from joblib import Parallel, delayed
import multiprocessing

#########for Chichi_yu
#name_col=0
#LL_col=[1,2]
#ENU_col=[3,4,5]
#sENU_col=[6,7,8]
#comp_INV=[0,1,2]
#scale_of_obs=1 # data*scale_of_obs=meter
#comp_INV=[0,1]
##################

#########for Chichi_Wu
#name_col=0
#LL_col=[1,2]
#ENU_col=[3,4,5]
#sENU_col=[6,7,8]
#comp_INV=[0,1,2]
#scale_of_obs=0.01 # data*scale_of_obs=meter

#########for CosDisp_Eq20180204_v2_gamit.txt
#name_col=8
#LL_col=[0,1]
#ENU_col=[4,2,6]
#sENU_col=[5,3,7]
#comp_INV=[0,1,2]
#comp_INV=[0,1]
#scale_of_obs=1000 # data*scale_of_obs=meter
##################

###########for fakeEQ###############
#name_col=0
#LL_col=[1,2]
#ENU_col=[3,4,5]
#sENU_col=[6,7,8]
#comp_INV=[0,1,2]
#scale_of_obs=1 
#####################################




def simple_dist(lon1,lat1,lon2,lat2):
    return ((lon1-lon2)**2+(lat1-lat2)**2)**0.5

def makeG(ngrid,use_idx,use_comp=[0,1,2]):
    #make G from .npy file
    #ngrid:name of grid file
    #use_idx: index of the G used
    #use_comp=[0,1,2] 0:E,1:N,2:U component you want to use in the inversion
    G=np.load(ngrid)
    n_sta=int(len(G)/3)
    all_idx_comps=np.array([],dtype='int')
    for icomp in use_comp:
        all_idx_comps=np.hstack([all_idx_comps,np.array(use_idx)+n_sta*icomp])
    G=G[all_idx_comps]
    return G

def makeG_dict(G,use_idx,use_comp=[0,1,2]):
    #make G from dictionary
    #G:G already loaded
    #use_idx: index of the G used
    #use_comp=[0,1,2] 0:E,1:N,2:U component you want to use in the inversion
    n_sta=int(len(G)/3)
    all_idx_comps=np.array([],dtype='int')
    for icomp in use_comp:
        all_idx_comps=np.hstack([all_idx_comps,np.array(use_idx)+n_sta*icomp])
    G=G[all_idx_comps]
    return G

def makeD(new_stations,Loc_col=[0,1],use_Dcol=[2,3,4],err=[5,6,7]):
    #arrange D to an array(1-D matrix) for the inversion
    #new_stations: station matrix
    #Loc_col=column of lon. lat.
    #use_Dcol:Data column order should be same as in use_comp in the (makeG) function
    #err:data error of use_col, put [] if direct inversion
    all_comp=[]
    for sta in new_stations:
        n_comp=[] #if use_col=ENU(2,3,4) in this example, n_comp=[E,N,U]
        n_D=[]
        n_sD=[]
        for LL in Loc_col:
            n_comp.append(sta[LL])
        for comp in use_Dcol:
            n_comp.append(sta[comp])
        for nerr in err:
            n_comp.append(sta[nerr])
        all_comp.append(n_comp)
    all_comp=np.array(all_comp) #now it should be [[lon,lat,E,N,U,sE,sN,sU],[lon,lat,E,N,U,sE,sN,sU],[lon,lat,E,N,U,sE,sN,sU]...]
    #stack the data
    stacked_D=[]
    for ncp in range(len(use_Dcol)):
        stacked_D=np.hstack([stacked_D,all_comp[:,2+ncp]]) #2 is the len(Loc_col)
    stacked_sD=[]
    for ncp in range(len(err)):
        stacked_sD=np.hstack([stacked_sD,all_comp[:,2+len(use_Dcol)+ncp]])
    return all_comp,stacked_D,stacked_sD

def GMD_solve(G,D):
    #calculate inversion here
    GT=np.transpose(G)
    M=np.dot(np.dot(np.linalg.inv(np.dot(GT,G)),np.transpose(G)),D)
    res=(np.sum((D-np.dot(G,M))**2)/len(D))**0.5 #RMS error
    return M,res

def find_STAidx(STA,stations,col_name,col_stlon,col_stlat,Grid,dist_filt):
    #observation station is the n idx in the STA
    #considering the situation of same STAname but different location, call (simple_dist) to check also
    #STA dictionary
    #stations:a matrix with observation data
        #for example
        #122.934600 24.456000 0.816667 -0.550000 4.243333 5.246613 4.493286 18.729340 0 0499
        #121.125680 23.132360 -1.056667 -1.476667 4.230000 7.534462 7.231036 33.314867 0 CHUL
    #col_name: #index of name column, in this example, 9
    #col_stlon:in this example, 0
    #col_stlat:in this example, 1
    #Grid: location of the Grid that you are working
    #dist_filt:filt the stations greater than this distance
    #OUT idx for G matrix, new stations
    all_idx=[]
    new_stations=[]
    for line in stations:
        try:
            stainfo=STA[line[col_name].decode()]
            idx=stainfo[0]
            check_dist=simple_dist(stainfo[1],stainfo[2],line[col_stlon],line[col_stlat])
            if check_dist<=0.05:
                check_dist2=simple_dist(stainfo[1],stainfo[2],Grid[0],Grid[1])
                if check_dist2<=dist_filt:
                    all_idx.append(idx)
                    new_stations.append(line)
                    continue
            else:
                #print(stainfo[1],stainfo[2],line[col_stlon],line[col_stlat],'dist=',check_dist)
                #print(line[col_name].decode(),'Too far away, maybe not the same station')
                continue
        except:
            #print('Cannot find',line[col_name].decode())
            continue #cannot find the name!
    return all_idx,new_stations



def Mij2sdr(Mij):
    if len(Mij)==6:
        Mij=np.matrix([[Mij[0],Mij[1],Mij[2]], [Mij[1],Mij[3],Mij[4]], [Mij[2],Mij[4],Mij[5]]])
    diag_m,vector=np.linalg.eig(Mij)
    diag_m=np.matrix([[diag_m[0],0,0],[0,diag_m[1],0],[0,0,diag_m[2]]])
    values_total=diag_m.copy()
    #--------------------------------------------------------------------------
    # isotropic part of the moment tensor
    #--------------------------------------------------------------------------
    volumetric = (1./3)*np.sum(values_total); #np.sum(diag_m) is the trace
    m_volumetric = volumetric*np.matrix([[1,0,0],[0,1,0],[0,0,1]])
    #--------------------------------------------------------------------------
    # deviatoric part of the moment tensor
    #--------------------------------------------------------------------------
    m_deviatoric = diag_m - m_volumetric
    #--------------------------------------------------------------------------
    # eigenvalues of the deviatoric part of the moment tensor
    #--------------------------------------------------------------------------
    value_dev,junk = np.linalg.eig(m_deviatoric)
    #--------------------------------------------------------------------------
    # calculation of a fault normal and slip from the moment tensor
    #--------------------------------------------------------------------------
    j=np.argsort(value_dev) #the sort index
    values=value_dev[j] #the sort value
    n1=vector[:,j[2]]+vector[:,j[0]] #max+min vector
    n1=n1/np.dot(n1.transpose(),n1)[0,0]**0.5 #n1=n1/norm(n1)
    if n1[2]>0:
        n1=n1*-1.0 #vertical component is always negative!(?
    u1=vector[:,j[2]]-vector[:,j[0]] #max-min vector
    u1=u1/np.dot(u1.transpose(),u1)[0,0]**0.5 #u1=u1/norm(u1)
    if n1.transpose()*Mij*u1+u1.transpose()*Mij*n1 < 0:
        u1=-1.0*u1
    n2 = u1.copy()
    u2 = n1.copy()
    if n2[2]>0:
        n2=-1.0*n2 #vertical component is always negative
    #--------------------------------------------------------------------------
    # 1st solution
    #--------------------------------------------------------------------------
    dip=math.acos(-n1[2])*180/np.pi
    strike=math.asin(-n1[0]/(n1[0,0]**2+n1[1,0]**2)**0.5)*180/np.pi
    #determination of a quadrant
    if (n1[1]<0):
        strike=180-strike
    rake=math.asin(-u1[2]/np.sin(dip*np.pi/180))*180/np.pi
    #determination of a quadrant
    cos_rake = u1[0]*np.cos(strike*np.pi/180)+u1[1]*np.sin(strike*np.pi/180)
    if cos_rake<0:
        rake=180-rake
    if strike<0:
        strike=strike+360
    if rake<-180:
        rake=rake+360
    if rake>180:
        rake=rake-360 #make rake -180~180
    strike1 = strike; dip1 = dip; rake1 = rake #this is the first solution
    #--------------------------------------------------------------------------
    # 2nd solution
    #--------------------------------------------------------------------------
    dip=math.acos(-n2[2])*180/np.pi
    strike=math.asin(-n2[0]/(n2[0,0]**2+n2[1,0]**2)**0.5)*180/np.pi
    #determination of a quadrant
    if (n2[1]<0):
        strike=180-strike
    rake=math.asin(-u2[2]/np.sin(dip*np.pi/180))*180/np.pi
    #determination of a quadrant
    cos_rake = u2[0]*np.cos(strike*np.pi/180)+u2[1]*np.sin(strike*np.pi/180)
    if cos_rake<0:
        rake=180-rake
    if strike<0:
        strike=strike+360
    if rake<-180:
        rake=rake+360
    if rake>180:
        rake=rake-360 #make rake -180~180
    strike2 = strike; dip2 = dip; rake2 = rake;
    return strike1,dip1,rake1,strike2,dip2,rake2


def sdr2Mij(strike,dip,rake,return_type=1):
    #->"['IN':'degree','OUT:','moment tensor','type=1':'M1-M6','type=2':'matrix']"
    delta=dip*np.pi/180
    lamda=rake*np.pi/180
    phi=strike*np.pi/180
    M11=-(np.sin(delta)*np.cos(lamda)*np.sin(2*phi)+np.sin(2*delta)*np.sin(lamda)*np.sin(phi)**2)
    M12=np.sin(delta)*np.cos(lamda)*np.cos(2*phi)+0.5*np.sin(2*delta)*np.sin(lamda)*np.sin(2*phi)
    M13=-(np.cos(delta)*np.cos(lamda)*np.cos(phi)+np.cos(2*delta)*np.sin(lamda)*np.sin(phi))
    M22=np.sin(delta)*np.cos(lamda)*np.sin(2*phi)-np.sin(2*delta)*np.sin(lamda)*np.cos(phi)**2
    M23=-(np.cos(delta)*np.cos(lamda)*np.sin(phi)-np.cos(2*delta)*np.sin(lamda)*np.cos(phi))
    M33=-(M11+M22)
    if return_type==1:
        return M11,M12,M13,M22,M23,M33
    elif return_type==2:
        return np.matrix([[M11,M12,M13],[M12,M22,M23],[M13,M23,M33]])
    else:
        print('return_type can only be 1 or 2')

#--------------------------START----------------------------------------------#

def load_G(ENUdir,GRDfile):
    #load GFs
    #all_grids=glob.glob(ENUdir+'/'+'*ENU.npy')
    GRD=np.genfromtxt(GRDfile)
    if GRD.ndim==1:
        all_grids=[ENUdir+'/'+'SRC%05d_ENU.npy'%(i+1) for i in range(1)]
    else:
        all_grids=[ENUdir+'/'+'SRC%05d_ENU.npy'%(i+1) for i in range(len(GRD))]
    #pre-load the GFs will save time!
    all_G={}
    print('Pre-loading GFs...')
    for i,Gname in enumerate(all_grids):
        all_G[i]=np.load(Gname)
    return all_G

def sta2dict(STAfile):
    #make station check dictionary
    STA={} #initial STA, make it a dictionary, STA={'staname':[idx,stlon,stlat]}
    IN1=open(STAfile,'r')
    for n,line in enumerate(IN1.readlines()):
        staname=line.split()[-1].split('_')[0]
        stainfo=[n,float(line.split()[0]),float(line.split()[1])]
        STA[staname]=stainfo
    IN1.close()
    return STA


def loop_inv(ng,ngrid,all_idx,stacked_D,stacked_sD,return_dhat):
    #ng:index of grid
    #ngrid: GF
    #return_dhat: True/False
    #print('running:',ng)
    #print('ngrid=',ngrid)
    G=makeG_dict(ngrid,all_idx,comp_INV)
    #weighted LSQ
    new_G=[]
    for n,g in enumerate(G):
        new_G.append(g/stacked_sD[n])
    new_D=stacked_D/stacked_sD
    new_G=np.vstack([new_G,np.array([1, 0, 0, 1, 0, 1])*1e-3 ])
    new_D=np.hstack([new_D,0])
    M,res=GMD_solve(new_G,new_D)
    if return_dhat:
        dhat=np.dot(G,np.matrix(M).transpose())
        return ng,res,M,dhat
    else:
        return ng,res,M

def run_inv(ENUdir,all_G,GRD,STA,data_file,outlog,n_cores,outfile):
    '''
        ENUdir: the directory of the *ENU file
        all_grids: pre-loaded G
        GRD: GRD from np.genfromtxt(GRDfile) that generated GFs
        STA: STA in dictionary
        data_file: ENU observations file
        outlog: output log file
        n_cores:number of cores for loop
        outfile:[file name/False] output the best dhat and residuals?
    '''
    #make station check dictionary
    #This should be fine by using the pre-loaded STA dictionary
    #STA={} #initial STA, make it a dictionary, STA={'staname':[idx,stlon,stlat]}
    #IN1=open(STAfile,'r')
    #for n,line in enumerate(IN1.readlines()):
    #    staname=line.split()[-1].split('_')[0]
    #    stainfo=[n,float(line.split()[0]),float(line.split()[1])]
    #    STA[staname]=stainfo
    #IN1.close()
    #load GRD loc file
    #GRD=np.genfromtxt(GRDfile)
    #Use all the stations
    #station file from observation
    stations=np.genfromtxt(data_file,dtype=None)
    all_idx,new_stations=find_STAidx(STA,stations,name_col,LL_col[0],LL_col[1],[121,23],10) #find the id in the G matrix, new_stations:filtered stations
    #make the D array
    #all_comp,stacked_D,stacked_sD=makeD(new_stations,Loc_col=[0,1],use_Dcol=[2,3,4],err=[5,6,7])
    all_comp,stacked_D,stacked_sD=makeD(new_stations,Loc_col=LL_col,use_Dcol=ENU_col,err=sENU_col) #for _v2_gamit.txt
    #plt.quiver(all_comp[:,0],all_comp[:,1],all_comp[:,2],all_comp[:,3])
    #plt.show()
    stacked_D=stacked_D*scale_of_obs #mm to meter
    stacked_sD=stacked_sD*scale_of_obs
    ###Inversion from every grids###
    if not(all_G):
        #load GFs
        #all_grids=glob.glob(ENUdir+'/'+'*ENU.npy')
        if GRD.ndim==1:
            all_grids=[ENUdir+'/'+'SRC%05d_ENU.npy'%(i+1) for i in range(1)]
        else:
            all_grids=[ENUdir+'/'+'SRC%05d_ENU.npy'%(i+1) for i in range(len(GRD))]
        #pre-load the GFs will save time!
        all_G={}
        print('Loading GFs...')
        for i,Gname in enumerate(all_grids):
            all_G[i]=np.load(Gname)
    else:
        if len(all_G)==len(GRD): #quick check if the preload_G input is right
            print('Given preload G...')
        else:
            if GRD.ndim==1:
                all_grids=[ENUdir+'/'+'SRC%05d_ENU.npy'%(i+1) for i in range(1)]
            else:
                all_grids=[ENUdir+'/'+'SRC%05d_ENU.npy'%(i+1) for i in range(len(GRD))]
            #pre-load the GFs will save time!
            all_G={}
            print('Dimension incorrect: Loading GFs...')
            for i,Gname in enumerate(all_grids):
                all_G[i]=np.load(Gname)
    #####run parallel loop#####
    results = Parallel(n_jobs=n_cores)(delayed(loop_inv)(ng,all_G[ng],all_idx,stacked_D,stacked_sD,False) for ng in all_G ) #or all_G.keys(). Note that all_G is dictionary
    #####run parallel loop END#####
    all_res=np.array([result[1] for result in results])
    min_idx=np.where(all_res==np.min(all_res))[0][0] #Check the Parallel keeps the same order!
    #print('min_idx=',min_idx)
    #print('min_results=',results[min_idx])
    if GRD.ndim==1:
        depth=GRD[2]
        if outfile!=False:
            OUTdepth_res=open(outfile+'_%f.res'%(depth),'w')
            OUTdepth_res.write('%f %f %f %d\n'%(GRD[0],GRD[1],all_res[0],0))
            OUTdepth_res.close()
            OUTall_res=open(outfile+'_all.res','w')
            OUTall_res.write('%f %f %f %f %d\n'%(GRD[0],GRD[1],GRD[2],all_res[0],0))
            OUTall_res.close()
            #output dhat in meter
            results_best = Parallel(n_jobs=n_cores)(delayed(loop_inv)(ng,all_G[ng],all_idx,stacked_D,stacked_sD,True) for ng in [0] )
            #print('Best solution=',results_best)
            dhat=results_best[0][3]
            nsta_used=len(all_comp) #number of stations used in the inversion, all_comp:matrix for [lon,lat,E,N,U],[lon,lat,E,N,U]
            print('#of stations=',nsta_used)
            OUTdhat=open(outfile+'_dhat.dat','w')
            for nn,sta_used in enumerate(new_stations):
                OUTdhat.write('%f %f  %f %f %f %f %f %f 0 %s\n'%(sta_used[0],sta_used[1],dhat[nn],dhat[int(nn+nsta_used)],dhat[int(nn+2*nsta_used)],0,0,0,sta_used[name_col].decode() ))
            OUTdhat.close()
            np.savetxt('checkdhat',dhat)
    else:
        depth=GRD[min_idx][2]
        depth_idx=np.where(GRD[:,2]==depth)[0]
        if outfile!=False:
            OUTdepth_res=open(outfile+'_%f.res'%(depth),'w')
            for nres in depth_idx:
                OUTdepth_res.write('%f %f %f %d\n'%(GRD[nres,0],GRD[nres,1],all_res[nres],nres))
            OUTdepth_res.close()
            OUTall_res=open(outfile+'_all.res','w')
            for nres in range(len(all_res)):
                OUTall_res.write('%f %f %f %f %d\n'%(GRD[nres,0],GRD[nres,1],GRD[nres,2],all_res[nres],nres))
            OUTall_res.close()
            #output dhat in meter
            results_best = Parallel(n_jobs=n_cores)(delayed(loop_inv)(ng,all_G[ng],all_idx,stacked_D,stacked_sD,True) for ng in [min_idx] )
            #print('Best solution=',results_best)
            dhat=results_best[0][3]
            nsta_used=len(all_comp) #number of stations used in the inversion, all_comp:matrix for [lon,lat,E,N,U],[lon,lat,E,N,U]
            print('#of stations=',nsta_used)
            OUTdhat=open(outfile+'_dhat.dat','w')
            for nn,sta_used in enumerate(new_stations):
                OUTdhat.write('%f %f  %f %f %f %f %f %f 0 %s\n'%(sta_used[0],sta_used[1],dhat[nn],dhat[int(nn+nsta_used)],dhat[int(nn+2*nsta_used)],0,0,0,sta_used[name_col].decode() ))
            OUTdhat.close()
            np.savetxt('checkdhat',dhat)
    
    sav_M=results[min_idx][2]
    M=sav_M*1e20; #this is the setting when I made GFs
    Mij=np.matrix([[M[0],M[1],M[2]], [M[1],M[3],M[4]], [M[2],M[4],M[5]]])
    strike1,dip1,rake1,strike2,dip2,rake2=Mij2sdr(Mij)
    m0=(M[0]**2+M[3]**2+M[5]**2+2*M[1]**2+2*M[2]**2+2*M[4]**2)**0.5/(2**0.5);
    mw=2/3*(np.log10(m0)-16.1)
    print('----------First solution from all stations----------')
    if GRD.ndim==1:
        print('Loc=',GRD[0])
    else:
        print('Loc=',GRD[min_idx,:])
    print('Mw=%f'%(mw))
    print('SDR=',strike1,dip1,rake1,strike2,dip2,rake2)
    '''
    #if you want to save time, just don't go to the 2nd round
    ###############2nd round Inversion from near-field stations only#################
    #load GFs
    #all_grids=glob.glob(ENUdir+'/'+'*ENU.npy')
    #all_grids=[ENUdir+'/'+'SRC%05d_ENU.npy'%(i+1) for i in range(len(GRD))]
    #sav_res=[]
    #sav_M=[]
    minres=1e20
    print('allG=')
    print('all_G=',all_G)
    for ng,ngrid in enumerate(all_G):
    #for ng,ngrid in enumerate([all_grids[sav_ng]]):
        #if ng%10000==0:
            #print('Second search:',ng,'out of',len(all_grids))
        dist_eq=simple_dist(GRD[sav_ng,0],GRD[sav_ng,1],GRD[ng,0],GRD[ng,1])
        if dist_eq>0.5:
            continue
        #station filtering
        #station file from observation
        stations=np.genfromtxt(data_file,dtype=None)
        if mw<6.5:
            sta_dist_filt=0.6
        elif mw>=6.5 and mw<7.0:
            sta_dist_filt=1.2
        else:
            sta_dist_filt=3.0
        print('dist=',sta_dist_filt)
        #all_idx,new_stations=find_STAidx(STA,stations,9,0,1,[GRD[sav_ng,0],GRD[sav_ng,1]],sta_dist_filt) #find the id in the G matrix, new_stations:new station file without issued stations
        all_idx,new_stations=find_STAidx(STA,stations,name_col,LL_col[0],LL_col[1],[GRD[sav_ng,0],GRD[sav_ng,1]],sta_dist_filt)
        #make the D array
        #all_comp,stacked_D,stacked_sD=makeD(new_stations,Loc_col=[0,1],use_Dcol=[2,3,4],err=[5,6,7])
        #all_comp,stacked_D,stacked_sD=makeD(new_stations,Loc_col=[0,1],use_Dcol=[4,2,6],err=[5,3,7]) #for _v2_gamit.txt
        all_comp,stacked_D,stacked_sD=makeD(new_stations,Loc_col=LL_col,use_Dcol=ENU_col,err=sENU_col) #for _v2_gamit.txt
        stacked_D=stacked_D*scale_of_obs #mm to meter
        stacked_sD=stacked_sD*scale_of_obs
        #G=makeG(ngrid=ngrid,use_idx=all_idx,use_comp=[0,1,2])
        G=makeG_dict(all_G[ng],use_idx=all_idx,use_comp=comp_INV)
        #weighted LSQ
        new_G=[]
        for n,g in enumerate(G):
            new_G.append(g/stacked_sD[n])
        new_D=stacked_D/stacked_sD
        new_G=np.vstack([new_G,np.array([1, 0, 0, 1, 0, 1])*1e0 ])
        new_D=np.hstack([new_D,0])
        M,res=GMD_solve(new_G,new_D)
        if res<minres:
            sav_M2=M.copy()
            sav_ng2=ng
            minres=res

    M2=sav_M2*1e20; #this is the setting when I made GFs
    Mij2=np.matrix([[M2[0],M2[1],M2[2]], [M2[1],M2[3],M2[4]], [M2[2],M2[4],M2[5]]])
    strike1,dip1,rake1,strike2,dip2,rake2=Mij2sdr(Mij2)
    m0=(M2[0]**2+M2[3]**2+M2[5]**2+2*M2[1]**2+2*M2[2]**2+2*M2[4]**2)**0.5/(2**0.5);
    mw=2/3*(np.log10(m0)-16.1)
    print('----------Second solution from close stations----------')
    print('Loc=',GRD[sav_ng2,:])
    print('Mw=%f'%(mw))
    print('SDR=',strike1,dip1,rake1,strike2,dip2,rake2)
    '''
    #write the inversion result
    OUT_inv=open(outlog,'w')
    OUT_inv.write('Inversion made on:%s\n'%(datetime.datetime.now()))
    OUT_inv.write('GFs from: %s\n'%(ENUdir))
    OUT_inv.write('GRDfile from: %s\n'%(GRDfile))
    OUT_inv.write('Data file: %s\n'%(data_file))
    if GRD.ndim==1:
        OUT_inv.write('EQloc: %f %f %f\n'%(GRD[0],GRD[1],GRD[2]))
    else:
        OUT_inv.write('EQloc: %f %f %f\n'%(GRD[min_idx,0],GRD[min_idx,1],GRD[min_idx,2]))
    OUT_inv.write('SDR: %f %f %f %f %f %f\n'%(strike1,dip1,rake1,strike2,dip2,rake2))
    OUT_inv.write('Mw: %f\n'%(mw))
    OUT_inv.close()


#all_G=load_G(ENUdir,GRDfile)
#run_inv(ENUdir,all_G,GRDfile,STAfile,data_file,'testQQQ.log')

#ENUdir='/Users/timlin/Documents/Project/GPSInv/GFs/TW_dense/ENU_out' #GFs for all GRD with respect to STA
#GRDfile='/Users/timlin/Documents/Project/GPSInv/GFs/TW_dense/TW_dense.grid'
#STAfile='/Users/timlin/Documents/Project/GPSInv/GFs/TW_dense/TPN_All_filt.sta'
#data_file='/Users/timlin/Documents/Project/GPSInv/code/Inversion/TEST00000/SRC_sum.enu'

#ENUdir='/home/jtlin/TWGPS/ENU_out' #GFs for all GRD with respect to STA
#GRDfile='/home/jtlin/TWGPS/TW_dense.grid'
#STAfile='/home/jtlin/TWGPS/TPN_All_filt.sta'
#data_files='/home/jtlin/Make_TWEQ/TW_fakeEQM6.0_7.3/TEST%05d/SRC_sum.enu'
###########for fakeEQ###############
#name_col=0
#LL_col=[1,2]
#ENU_col=[3,4,5]
#sENU_col=[6,7,8]
#comp_INV=[0,1,2]
#scale_of_obs=1
#n_cores=4 #parallelizing loop for gridsearch best source location
#####################################



##########These example variable will be changed from outside when you set the right path in the parameters control file##########
#########for Coseis_xxxx.gam
ENUdir='/Users/timlin/Documents/Project/GPSInv/GFs/TW_dense/ENU_out' #GFs for all GRD with respect to STA
GRDfile='/Users/timlin/Documents/Project/GPSInv/GFs/TW_dense/TW_dense.grid'
STAfile='/Users/timlin/Documents/Project/GPSInv/GFs/TW_dense/TPN_All_filt.sta'
#data_file='/Users/timlin/Documents/Project/GPSInv/Invercases/Coseismic_data/Coseis_0602.gam'
data_file='/Users/timlin/Documents/Project/GPSInv/Invercases/Coseismic_data/Coseis_1226.gam'
name_col=9
LL_col=[0,1]
ENU_col=[2,3,4]
sENU_col=[5,6,7]
comp_INV=[0,1,2]
scale_of_obs=0.001 # data*scale_of_obs=meter
n_cores=2
##################



def Main_run():
    all_G=load_G(ENUdir,GRDfile)
    GRD=np.genfromtxt(GRDfile)
    STA=sta2dict(STAfile)
    run_inv(ENUdir,all_G,GRD,STA,data_file,'GPSCMT.log',n_cores,'GPSCMT')



#all_G=load_G(ENUdir,GRDfile)
#GRD=np.genfromtxt(GRDfile)
#STA=sta2dict(STAfile)
#run_inv(ENUdir,all_G,GRD,STA,data_file,'Hengchun1226.log',n_cores,'Hengchun1226')



#import time
#time1=time.time()
#data_file='/Users/timlin/Documents/Project/GPSInv/code/Inversion/TEST00000/SRC_sum.enu'
#run_inv(ENUdir,all_G,GRD,STA,data_file,'test00000.log',n_cores)
#print('Running time=',time.time()-time1)
#time1=time.time()
#data_file='/Users/timlin/Documents/Project/GPSInv/code/Inversion/TEST00001/SRC_sum.enu'
#run_inv(ENUdir,all_G,GRD,STA,data_file,'test00001.log',n_cores)
#print('Running time=',time.time()-time1)
#time1=time.time()
#data_file='/Users/timlin/Documents/Project/GPSInv/code/Inversion/TEST00002/SRC_sum.enu'
#run_inv(ENUdir,all_G,GRD,STA,data_file,'test00002.log',n_cores)
#print('Running time=',time.time()-time1)


