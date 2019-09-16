#Parameter file for GPS CMT

#########################Parameters only need to be changed here#############################
#make forward model by fk
#---------------------parameters for forward model---------------------------
home='/Users/timlin/Documents/Project/GPSCMT/example/Nantou0602' #where's your working dir?
project_name='Test1' #create a folder under working dir
fk_home='/Users/timlin/Documents/Project/GPSInv/code/make_GFs/fk' #where is your fk dir
modelpath='/Users/timlin/Documents/Project/GPSCMT/example/Nantou0602/vel1D_CWB.mod' #full path of the velocity model

grdfile='Nantou0602.grd' #source location file should be under "home" directory
stafile='Nantou0602.sta' #recording stations
srcfile=None             #source file (doesn't matter in generating GFs step)

nprocess=8    #parallel processing
SRCtype=3 #type=1 (Mw/strike/dip/rake) or type=2 (m0/M1~M6) or type3 (generate Green's function for CMT)
GFs_out=1 # should be 1~3.  1.make .npy data, 2.make .txt data(for checking), 3.both

#-----------------------parameters for CMT inversion-------------------------
data_file='/Users/timlin/Documents/Project/GPSCMT/example/Nantou0602/Coseis_0602.gam' #data path 
name_col=9            #open your data file. at which column is your station name? 
LL_col=[0,1]          #at which column is the station lon. and lat. ?
ENU_col=[2,3,4]       #which column is E,N,U? 
sENU_col=[5,6,7]      #what are the corresponding std for ENU?
comp_INV=[0,1,2]      #component of inversion [0,1,2] for [E,N,U]. Similarly, [0,1] for just [E,N]
scale_of_obs=0.001    #data*scale_of_obs=meter (i.e. the unit of data are in mm so it needs to multiply by 0.001 to meter)
n_cores=2             #how many CPUs you are using? CAREFULLY USED when your GFs is HUGE! (10,000+ grid points, 500+ stations)

#-----------------------Switches, whether to run--------------------------------------
#Those True(or False) decide whether run the process
make_GreensFcn=0      #Generate GFs?
convert_GFs=0          #Convert ZRT to ENZ and save in npy format
run_CMTinv=1          #Run CMT inversion

##################################Parameters setting END######################################



if make_GreensFcn:
    from gpscmt import make_forward  ####If this fail, check if you add python path ->  export PYTHONPATH=$PYTHONPATH:YOUR_PATH_TO_GPS_CMT/src/python in environment
    make_forward.home=home
    make_forward.project_name=project_name
    make_forward.fk_home=fk_home
    make_forward.modelpath=modelpath
    make_forward.grdfile=grdfile
    make_forward.srcfile=srcfile
    make_forward.stafile=stafile
    make_forward.nprocess=nprocess
    make_forward.SRCtype=SRCtype
    make_forward.run_forward()

if convert_GFs:
    from gpscmt import zrt2enz
    inpdir=home+'/'+project_name+'/'
    outdir=home+'/'+project_name+'/GFs/'
    FMT=GFs_out #out
    zrt2enz.inpdir=inpdir
    zrt2enz.outdir=outdir
    zrt2enz.FMT=FMT
    zrt2enz.run_convert()

if run_CMTinv:
    try:
        from gpscmt import main_inv_mpi_detailed
        ENUdir=home+'/'+project_name+'/GFs' #GFs for all GRD with respect to STA(Check whether this is correct if you didn't run convert_GFs above)
        GRDfile=home+'/'+grdfile
        STAfile=home+'/'+stafile
        main_inv_mpi_detailed.ENUdir=ENUdir
        main_inv_mpi_detailed.GRDfile=GRDfile
        main_inv_mpi_detailed.STAfile=STAfile
        main_inv_mpi_detailed.data_file=data_file
        main_inv_mpi_detailed.name_col=name_col
        main_inv_mpi_detailed.LL_col=LL_col
        main_inv_mpi_detailed.ENU_col=ENU_col
        main_inv_mpi_detailed.sENU_col=sENU_col
        main_inv_mpi_detailed.comp_INV=comp_INV
        main_inv_mpi_detailed.scale_of_obs=scale_of_obs
        main_inv_mpi_detailed.n_cores=n_cores
        main_inv_mpi_detailed.Main_run()
    except:
        print('Single core not available in this version...')
        print('please >>pip install joblib and pip install multiprocessing')
    

