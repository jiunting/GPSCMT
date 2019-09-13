#plot result for GPSCMT

def res_map(res_file,model_file,data_file):
    import matplotlib.pyplot as plt
    import numpy as np
    res=np.genfromtxt(res_file)
    #residuals map
    plt.scatter(res[:,0],res[:,1],c=res[:,2],cmap=plt.cm.jet_r)
    plt.colorbar(label='residual')
    #vectors
    d=np.genfromtxt(data_file)
    dhat=np.genfromtxt(model_file)
    plt.quiver(d[:,0],d[:,1],d[:,2]*0.001,d[:,3]*0.001,scale=0.1,color=[0,0,0]) #data
    plt.quiver(dhat[:,0],dhat[:,1],dhat[:,2],dhat[:,3],scale=0.1,color=[1,0,0]) #model
    plt.xlabel('Lon.',fontsize=15)
    plt.ylabel('Lat.',fontsize=15)
    plt.show()
