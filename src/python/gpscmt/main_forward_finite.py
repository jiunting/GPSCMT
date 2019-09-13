#make forward ENU at given stations from random finite fault
import numpy as np
import matplotlib.pyplot as plt
import math
import shutil
import glob
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

def makefault(fout,strike,dip,nstrike,dx_dip,dx_strike,epicenter,num_updip,num_downdip,rise_time):
    '''
    Copied from Mudpy
    Make a planar fault
        
    strike - Strike angle (degs)
    dip - Dip angle (degs)200/5
    '''
    from numpy import arange,sin,cos,deg2rad,r_,ones,arctan,rad2deg,zeros,isnan,unique,where,argsort
    import pyproj
    
    proj_angle=180-strike #Angle to use for sin.cos projection (comes from strike)
    y=arange(-nstrike/2+1,nstrike/2+1)*dx_strike
    x=arange(-nstrike/2+1,nstrike/2+1)*dx_strike
    z=ones(x.shape)*epicenter[2]
    y=y*cos(deg2rad(strike))
    x=x*sin(deg2rad(strike))
    #Save teh zero line
    y0=y.copy()
    x0=x.copy()
    z0=z.copy()
    #Initlaize temp for projection up/down dip
    xtemp=x0.copy()
    ytemp=y0.copy()
    ztemp=z0.copy()
    #Get delta h and delta z for up/ddx_dip=1own dip projection
    if num_downdip>0 and num_updip>0:
        dh=dx_dip*cos(deg2rad(dip))
        dz=dx_dip*sin(deg2rad(dip))
        #Project updip lines
        for k in range(num_updip):
            xtemp=xtemp+dh*cos(deg2rad(proj_angle))
            ytemp=ytemp+dh*sin(deg2rad(proj_angle))
            ztemp=ztemp-dz
            x=r_[x,xtemp]
            y=r_[y,ytemp]
            z=r_[z,ztemp]
        #Now downdip lines
        xtemp=x0.copy()
        ytemp=y0.copy()
        ztemp=z0.copy()
        for k in range(num_downdip):
            xtemp=xtemp-dh*cos(deg2rad(proj_angle))
            ytemp=ytemp-dh*sin(deg2rad(proj_angle))
            ztemp=ztemp+dz
            x=r_[x,xtemp]
            y=r_[y,ytemp]
            z=r_[z,ztemp]
    #Now use pyproj to dead reckon anf get lat/lon coordinates of subfaults
    g = pyproj.Geod(ellps='WGS84')
    #first get azimuths of all points, go by quadrant
    az=zeros(x.shape)
    for k in range(len(x)):
        if x[k]>0 and y[k]>0:
            az[k]=rad2deg(arctan(x[k]/y[k]))
        if x[k]<0 and y[k]>0:
            az[k]=360+rad2deg(arctan(x[k]/y[k]))
        if x[k]<0 and y[k]<0:
            az[k]=180+rad2deg(arctan(x[k]/y[k]))
        if x[k]>0 and y[k]<0:
            az[k]=180+rad2deg(arctan(x[k]/y[k]))
    #Quadrant correction
    #Now horizontal distances
    d=((x**2+y**2)**0.5)*1000
    #Now reckon
    lo=zeros(len(d))
    la=zeros(len(d))
    for k in range(len(d)):
        if isnan(az[k]): #No azimuth because I'm on the epicenter
            print('Point on epicenter')
            lo[k]=epicenter[0]
            la[k]=epicenter[1]
        else:
            lo[k],la[k],ba=g.fwd(epicenter[0],epicenter[1],az[k],d[k])
    #Sort them from top right to left along dip
    zunique=np.unique(z)
    for k in range(len(zunique)):
        i=where(z==zunique[k])[0] #This finds all faults at a certain depth
        isort=argsort(la[i]) #This sorths them south to north
        if k==0: #First loop
            laout=la[i][isort]
            loout=lo[i][isort]
            zout=z[i][isort]
        else:
            laout=r_[laout,la[i][isort]]
            loout=r_[loout,lo[i][isort]]
            zout=r_[zout,z[i][isort]]
    #Write to file
    strike=ones(loout.shape)*strike
    dip=ones(loout.shape)*dip
    tw=ones(loout.shape)*0.5
    rise=ones(loout.shape)*rise_time
    L=ones(loout.shape)*dx_strike*1000
    W=ones(loout.shape)*dx_dip*1000
    f=open(fout,'w')
    for k in range(len(x)):
        out='%i\t%.6f\t%.6f\t%.3f\t%.2f\t%.2f\t%.1f\t%.1f\t%.2f\t%.2f\n' % (k+1,loout[k],laout[k],zout[k],strike[k],dip[k],tw[k],rise[k],L[k],W[k])
        f.write(out)
    f.close()

def where_fault(Vmodel,depth):
    accumD=np.cumsum(Vmodel[:,0])
    idx=np.where( depth<=accumD )[0][0] #make sure there is no source depth larger than accumD[-1]!
    return(Vmodel[idx,:])


def gen_random_fault(num_columns,number_updip,number_downdip,dx_strike,dx_dip,target_Mw):
    ##################random finite fault generating#########################
    #using the below defaults
    #Length=40km width=20km
    #Length=30km width=15km
    Inland=np.genfromtxt('Taiwan_inland.dat')
    polygon = Polygon(Inland)
    #num_columns=10
    #number_updip=2
    #number_downdip=2
    #dx_strike=3
    #dx_dip=3
    rise_time=0.0 #rise_time is not important for the static model
    #target_Mw=6.8 #Mw6.5
    target_M0=10**(target_Mw*1.5+9.05)
    while True:
        fault_file='testEQ.fault'
        strike=np.random.rand()*359 #0~359
        dip=np.random.rand()*90 #0~90
        hypocenter=[np.random.rand()*2+120,np.random.rand()*3.5+22,np.random.randint(2,15)] #120~122 22.0~25.5
        #hypocenter=[np.random.rand()*1+120.3,np.random.rand()*0.8+23,np.random.randint(2,15)] #120~122 22.0~25.5
        makefault(fault_file,strike,dip,num_columns,dx_strike,dx_dip,hypocenter,number_updip,number_downdip,rise_time)
        A=np.genfromtxt(fault_file)
        #check if the top is in the air
        if np.min(A[:,3])<0:
            continue
        else:
            #Inland EQ only
            #print('Hypocenter=',hypocenter)
            if polygon.contains(Point(hypocenter[0],hypocenter[1]))==True:
                break
            else:
                continue
    ##################random finite fault geometry done#########################
    ################start generating slip on finite fault ######################
    G=[] #matrix of Laplacian
    for fault_id in A[:,0]:
        ncorners=0
        #make G
        g=[0]*len(A[:,0])
        #test whether the point has any neightbor
        #the simplist idea is to know
        if (fault_id-num_columns) in A[:,0]:
            ncorners+=1 #up point
            g[int(fault_id-num_columns)-1]=1
        if (fault_id+num_columns) in A[:,0]:
            ncorners+=1 #down point
            g[int(fault_id+num_columns)-1]=1
        if ((fault_id+1) in A[:,0]):
            #this should also be in the same row/column
            #get the col row # for new point
            nrow_r=(fault_id+1)//num_columns
            if (fault_id+1)%num_columns !=0:
                nrow_r+=1
            ncol_r=(fault_id+1)%num_columns
            if ncol_r==0:
                ncol_r=num_columns
            #get the #of col/row for current point
            nrow_cur=(fault_id)//num_columns
            if fault_id%num_columns !=0:
                nrow_cur+=1
            ncol_cur=(fault_id)%num_columns
            if ncol_cur==0:
                ncol_cur=num_columns
            if (nrow_cur==nrow_r) or (ncol_cur==ncol_r):
                ncorners+=1 #right point
                g[int(fault_id+1)-1]=1
                #print(fault_id,'has right boundary!')
        if ((fault_id-1) in A[:,0]):
            #this should also be in the same row/column
            #get the col row # for new point
            nrow_r=(fault_id-1)//num_columns
            if (fault_id-1)%num_columns !=0:
                nrow_r+=1
            ncol_r=(fault_id-1)%num_columns
            if ncol_r==0:
                ncol_r=num_columns
            #get the #of col/row for current point
            nrow_cur=(fault_id)//num_columns
            if fault_id%num_columns !=0:
                nrow_cur+=1
            ncol_cur=(fault_id)%num_columns
            if ncol_cur==0:
                ncol_cur=num_columns
            if (nrow_cur==nrow_r) or (ncol_cur==ncol_r):
                ncorners+=1 #right point
                g[int(fault_id-1)-1]=1
                #print(fault_id,'has left boundary!')
        #coef of the fault_id itself
        g[int(fault_id)-1]=ncorners
        G.append(g)
    Gm=np.matrix(G)
    cent_n=len(A[:,0])//2
    SS=[0]*len(A[:,0])
    init_SS=np.random.rand()*2-1 #-1~1
    SS[cent_n]=init_SS
    SS=np.matrix(SS).transpose()
    DS=[0]*len(A[:,0])
    init_DS=np.random.rand()*2-1  #-1~1
    DS[cent_n]=init_DS
    DS=np.matrix(DS).transpose()
    rake=math.atan2(init_DS,init_SS)*180/np.pi
    maxiter=np.max([number_updip,number_downdip,num_columns//2])
    for i in range(maxiter+6):
        SS=Gm*SS
        DS=Gm*DS
    #the rake should be same
    #tmpid1=np.where(np.abs(DS)==np.max(np.abs(DS)))[0][0]
    #tmpid2=np.where(np.abs(SS)==np.max(np.abs(SS)))[0][0]
    #math.atan2(DS[tmpid1],SS[tmpid2])*180/np.pi
    norm=np.max((np.array(SS)**2+np.array(DS)**2)**0.5)
    SS=SS/norm #normalize
    DS=DS/norm
    #plt.subplot(1,2,1)
    #plt.scatter(A[:,1],A[:,2],c=np.array(SS).reshape(-1),cmap=plt.cm.jet)
    #plt.colorbar()
    #plt.title('SS',fontsize=14)
    #plt.subplot(1,2,2)
    #plt.scatter(A[:,1],A[:,2],c=np.array(DS).reshape(-1),cmap=plt.cm.jet)
    #plt.colorbar()
    #plt.title('DS',fontsize=14)
    #plt.savefig('TestEQ_SSDS.png')
    #plt.show()
    #plt.figure()
    GRD=np.genfromtxt('/home/jtlin/TWGPS/TW_dense.grid')
    #plt.plot(GRD[:,0],GRD[:,1],'k.',markersize=1)
    GPSfile='/home/jtlin/TWGPS/TPN_All_filt.sta'
    STA=np.genfromtxt(GPSfile)
    #plt.plot(STA[:,0],STA[:,1],'b^')
    top_fault_idx=np.where(A[:,3]==np.min(A[:,3]))[0]
    #plt.plot(A[top_fault_idx,1],A[top_fault_idx,2],'ro')
    #plt.quiver(A[:,1],A[:,2],np.array(SS).reshape(-1),np.array(DS).reshape(-1))
    #plt.text(np.min(A[:,1]),np.max(A[:,2]),'strike=%5.2f'%(strike))
    #plt.text(np.min(A[:,1]),np.max(A[:,2])-0.05,'dip=%5.2f'%(dip))
    #plt.text(np.min(A[:,1]),np.max(A[:,2])-0.1,'rake=%5.2f'%(rake))
    #plt.savefig('TestEQ_fault.png')
    #plt.show()
    #plt.close() #close figure
    #########input the fault model to fk#############
    #STA=np.genfromtxt('/Users/timlin/Documents/Project/GPSInv/code/Inversion/TPN_All_filt.sta')
    Vmodel=np.genfromtxt('/home/jtlin/TWGPS/vel1D_CWB')
    #calculate Mw
    Area=A[:,8]*A[:,9] #m^2
    vs=np.array([where_fault(Vmodel,d)[1] for d in A[:,3]])
    rho=np.array([where_fault(Vmodel,d)[3] for d in A[:,3]])
    mu=vs**2*rho*1e9 #pa.  Vs=sqrt(mu/rho)
    Disp=(np.array(SS).reshape(-1)**2+np.array(DS).reshape(-1)**2)**0.5 #displacement
    M0=mu*Area*Disp
    ratio=target_M0/np.sum(M0)
    #scale the slip to the right magnitude
    Disp=(np.array(SS*ratio).reshape(-1)**2+np.array(DS*ratio).reshape(-1)**2)**0.5 #displacement rescale to the target Mw
    M0=mu*Area*Disp
    Mw=(np.log10(M0)-9.05)/1.5
    M0_sum=np.sum(mu*Area*Disp)
    Mw_sum=(np.log10(M0_sum)-9.05)/1.5
    #print('Mw=',Mw_sum)
    return strike,dip,rake,Mw,Mw_sum,hypocenter


#########make gridfile, srcfile, and stafile##############
##############set parameters################
#For 6.0: 15km x 10km
#For 6.5: 30km x 10km #this should be regenerated
#For 6.8: 40km x 15km
#For 7.0: 50km x 15km
#For 7.3-7.5: 80km x 15km
'''
target_Mw=np.random.rand() * 1.3  + 6.0
if target_Mw<6.5:    
    num_columns=5
    number_updip=2
    number_downdip=2
    dx_strike=3
    dx_dip=2
elif target_Mw>=6.5 and target_Mw<6.8:
    num_columns=10
    number_updip=2
    number_downdip=2
    dx_strike=3
    dx_dip=2
elif target_Mw>=6.8 and target_Mw<7.0:
    num_columns=10
    number_updip=2
    number_downdip=2
    dx_strike=4
    dx_dip=3
elif target_Mw>=7.0 and target_Mw<7.3:
    num_columns=10
    number_updip=2
    number_downdip=2
    dx_strike=5
    dx_dip=3
elif target_Mw>=7.3:
    num_columns=16
    number_updip=2
    number_downdip=2
    dx_strike=5
    dx_dip=3
'''
rise_time=0.0 #rise_time is not important for the static model
#target_Mw=6.5
#############Generate fault#####################

for i_run in range(10000):
    #generate random Mw
    target_Mw=np.random.rand() * 1.5  + 6.0
    if target_Mw<6.5:
        num_columns=5
        number_updip=2
        number_downdip=2
        dx_strike=3
        dx_dip=2
    elif target_Mw>=6.5 and target_Mw<6.8:
        num_columns=10
        number_updip=2
        number_downdip=2
        dx_strike=3
        dx_dip=2
    elif target_Mw>=6.8 and target_Mw<7.0:
        num_columns=10
        number_updip=2
        number_downdip=2
        dx_strike=4
        dx_dip=3
    elif target_Mw>=7.0 and target_Mw<7.2:
        num_columns=10
        number_updip=2
        number_downdip=2
        dx_strike=5
        dx_dip=3
    elif target_Mw>=7.2:
        num_columns=16
        number_updip=2
        number_downdip=2
        dx_strike=5
        dx_dip=3

    strike,dip,rake,Mw,Mw_sum,hypocenter=gen_random_fault(num_columns,number_updip,number_downdip,dx_strike,dx_dip,target_Mw)
    A=np.genfromtxt('testEQ.fault')

    OUT_grid=open('testEQ.grid','w')
    for nsubf in A:
        OUT_grid.write('%f %f %f\n'%(nsubf[1],nsubf[2],nsubf[3]))
    OUT_grid.close()

    OUT_src=open('testEQ.src','w')
    for i,nsubf in enumerate(A):
        OUT_src.write('%f %f %f %f\n'%(Mw[i],strike,dip,rake))
    OUT_src.close()

    STA=np.genfromtxt('/home/jtlin/TWGPS/TPN_All_filt.sta',dtype=None)
    stlon=np.array([i[0] for i in STA])
    stlat=np.array([i[1] for i in STA])
    stname=np.array([i[2].decode() for i in STA])
    close_idx=np.where( ((stlon-hypocenter[0])**2 + (stlat-hypocenter[1])**2)**0.5 <=1.5 )[0] #filt the stations within 1.5 degree
    OUT_sta=open('testEQ.sta','w')
    for id in close_idx:
        OUT_sta.write('%f %f %s\n'%(stlon[id],stlat[id],stname[id]))
    OUT_sta.close()

    ############run the forward model##############
    import make_forward
    make_forward.home='/home/jtlin/Make_TWEQ'
    make_forward.project_name='TEST%05d'%(i_run)
    make_forward.fk_home='/home/jtlin/fk' #this is the default path
    make_forward.grdfile='testEQ.grid'
    make_forward.srcfile='testEQ.src'
    #make_forward.stafile='/Users/timlin/Documents/Project/GPSInv/code/Inversion/TPN_All_filt.sta'
    make_forward.stafile='testEQ.sta' #only calculate ZRT on the filtered stations
    make_forward.modelpath='/home/jtlin/TWGPS/vel1D_CWB'
    make_forward.nprocess=40
    make_forward.SRCtype=1
    make_forward.run_forward()

    ############read the zrt file and convert to enu data##############
    import zrt2enz
    inpdir=make_forward.home+'/'+make_forward.project_name
    outdir=inpdir
    STA_path='testEQ.sta'
    FMT=0 #output format 1:.npy 2:.txt 3:both or 0:nothing
    sum_data=True #sum the SRCxxxx.zrt from finite fault to a single file (this is the observation)
    SRCtype=1
    zrt2enz.zrt2enz(inpdir=inpdir,outdir=outdir,STA_path=STA_path,FMT=FMT,sum_data=sum_data,SRCtype=SRCtype)

    ##############write summary file#########################
    OUT2=open(outdir+'/testEQ_summary.log','w')
    OUT2.write('Hypocenter: %f %f %f\n'%(hypocenter[0],hypocenter[1],hypocenter[2]))
    OUT2.write('Mw: %f\n'%(Mw_sum))
    OUT2.write('SDR: %f %f %f\n'%(strike,dip,rake))
    #OUT2.write('Number of stations:%d'%(len(stname_filt)))
    OUT2.write('num_columns: %d\n'%(num_columns))
    OUT2.write('number_updip: %d\n'%(number_updip))
    OUT2.write('number_downdip: %d\n'%(number_downdip))
    OUT2.write('dx_strike: %f\n'%(dx_strike))
    OUT2.write('dx_dip: %f\n'%(dx_dip))
    OUT2.close()

    ############move all the files in the directory##############
    shutil.move('testEQ.grid',make_forward.home+'/'+make_forward.project_name+'/testEQ.grid')
    shutil.move('testEQ.src',make_forward.home+'/'+make_forward.project_name+'/testEQ.src')
    shutil.move('testEQ.sta',make_forward.home+'/'+make_forward.project_name+'/testEQ.sta')
    shutil.move('testEQ.fault',make_forward.home+'/'+make_forward.project_name+'/testEQ.fault')
    #shutil.move('TestEQ_SSDS.png',make_forward.home+'/'+make_forward.project_name+'/TestEQ_SSDS.png')
    #shutil.move('TestEQ_fault.png',make_forward.home+'/'+make_forward.project_name+'/TestEQ_fault.png')

