#rotate the .zrt file to enz component and save it to .npy matrix
'''
Input:
 directory of the *.zrt
 format of .zrt:
 dist az Z1 R1 T1 Z2 R2 T2 ... Z6 R6 T6
Output:
 save the file in the format below:
 output format: .npy (or/and) .txt
 *.npy is a matrix with size [3*nsta, 6mt]
 .txt is text file with size the same with input file(this is for checking only)
'''
import numpy as np
import glob
import os

##########Donn't change thise########
inpdir='/Users/timlin/Documents/Project/GPSInv/GFs/TW_dense/GreensFcn/'
outdir='/Users/timlin/Documents/Project/GPSInv/GFs/TW_dense/ENU_out/'
SRCtype=3 #type=1 (Mw/strike/dip/rake) or type=2 (m0/M1~M6) or type3 (generate Green's function for CMT)
FMT=1 #output format 1:.npy 2:.txt 3:both
##########Donn't change thise##########

def run_convert():
    print('Start converting ZRT file to ENZ...')
    SRCs=glob.glob(inpdir+'*.zrt')
    SRCs.sort()
    nSRC=0
    if not(os.path.exists(outdir)):
        os.makedirs(outdir)
    for SRC in SRCs:
        if nSRC%10==0:
            print(nSRC,'out of',len(SRCs))
        A=np.genfromtxt(SRC)
        dists=A[:,0]
        azs=A[:,1]
        azrads=azs*np.pi/180
        sav_E_sta=[]
        sav_N_sta=[]
        sav_U_sta=[]
        if SRCtype==3:
            #zrt is in GF format
            for nsta,azrad in enumerate(azrads):
                Trans=np.matrix([[np.sin(azrad),np.cos(azrad)],[np.cos(azrad),np.sin(azrad)*-1.0]])
                sav_E_M16=[]
                sav_N_M16=[]
                sav_U_M16=[]
                for n_mt in range(6):
                    AA=np.matrix(A[nsta,3+n_mt*3:3+n_mt*3+2])
                    EN=Trans*AA.transpose()
                    #the EN are for different moment tensors
                    sav_E_M16.append(EN[0,0]*1e-2) #cm to meter
                    sav_N_M16.append(EN[1,0]*1e-2)
                    sav_U_M16.append(A[nsta,2+n_mt*3]*1e-2)
                sav_E_sta.append(sav_E_M16)
                sav_N_sta.append(sav_N_M16)
                sav_U_sta.append(sav_U_M16)
        elif SRCtype==1 or SRCtype==2:
            #zrt is in forward model format
            for nsta,azrad in enumerate(azrads):
                Trans=np.matrix([[np.sin(azrad),np.cos(azrad)],[np.cos(azrad),np.sin(azrad)*-1.0]])
                AA=np.matrix(A[nsta,3:5])
                EN=Trans*AA.transpose()
                sav_E_sta.append([EN[0,0]*1e-2]) #cm to meter
                sav_N_sta.append([EN[1,0]*1e-2])
                sav_U_sta.append([A[nsta,2]*1e-2])
        E=np.matrix(sav_E_sta)
        N=np.matrix(sav_N_sta)
        U=np.matrix(sav_U_sta)
        #stack all the ENU into a matrix
        G=np.vstack([E,N,U])
        nSRC += 1
        if nSRC==1:
            #write README file so that you will know the columns later
            if SRCtype==3:
                RdMe=open(outdir+'README_Format.txt','w')
                RdMe.write('The file is sorted in ENU.\n 6 columns mean M1~M6, unit in meter.')
                RdMe.write('For example\n %s\n'%('-2.626887e-10 -9.267290e-10 -9.669541e-11 -1.501915e-10 -8.642910e-11 1.737191e-10'))
                RdMe.write('The 6 column means from the source #FILE_NAME to the n-idx station GF')
                RdMe.close()
            elif SRCtype==1 or SRCtype==2:
                RdMe=open(outdir+'README_Format.txt','w')
                RdMe.write('The file is sorted in ENU.\n only one column here, unit in meter.')
                RdMe.write('From the source #FILE_NAME to the n-idx station forward model')
                RdMe.close()
        if FMT==1:
            outname=SRC.split('/')[-1].split('.')[0]+'_ENU.npy'
            np.save(outdir+outname,G)
        if FMT==2:
            outname=SRC.split('/')[-1].split('.')[0]+'_ENU.txt'
            OUT1=open(outdir+outname,'w')
            for g in G:
                np.savetxt(OUT1,g,fmt='%e')
            OUT1.close()
        if FMT==3:
            outname=SRC.split('/')[-1].split('.')[0]+'_ENU.npy'
            np.save(outdir+outname,G)
            outname=SRC.split('/')[-1].split('.')[0]+'_ENU.txt'
            OUT1=open(outdir+outname,'w')
            for g in G:
                np.savetxt(OUT1,g,fmt='%e')
            OUT1.close()


# Test code
if __name__ == "__main__":
    print('Usage: set','zrt2enz.[inpdir/outdir/OutFormat]')
    print('Run model:','zrt2enz.run_convert()')
