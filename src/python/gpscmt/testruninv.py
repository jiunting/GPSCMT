import main_inv

ENUdir='/Users/timlin/Documents/Project/GPSInv/GFs/TW_dense/ENU_out' #GFs for all GRD with respect to STA
GRDfile='/Users/timlin/Documents/Project/GPSInv/GFs/TW_dense/TW_dense.grid'
STAfile='/Users/timlin/Documents/Project/GPSInv/GFs/TW_dense/TPN_All_filt.sta'
#data_file='/Users/timlin/Documents/Project/GPSInv/Invercases/Chi_Chi/Wu_coseis.dat'
data_file='/Users/timlin/Documents/Project/GPSInv/code/Inversion/TEST2/SRC_sum.txt'

all_G=main_inv.load_G(ENUdir,GRDfile)

main_inv.run_inv(ENUdir,all_G,GRDfile,STAfile,data_file,'testQQQ.log')
#main_inv.run_inv(ENUdir,False,GRDfile,STAfile,data_file,'testQQQ.log')
