# GPSCMT
![GitHub last commit](https://img.shields.io/github/last-commit/jiunting/GPSCMT?style=plastic) 
### The forward and inverse code for 1-D point source

* Lin, J. T., Chang, W. L., Melgar, D., Thomas, A., & Chiu, C. Y. (2019). Quick determination of earthquake source parameters from GPS measurements: a study of suitability for Taiwan. Geophysical Journal International, 219(2), 1148-1162.  
****
What it can/cannot do
```
-[x] Forward model: One source to single/many stations 
    e.g. You have a point source, what is the surface displacement looks like?
-[x] Forward model: Many sources to single/many stations
    e.g. You have many sources, what are the corresponding surface displacements?
-[x] Generate Green's function (Strike/Dip slip or moment tensors)
-[x] Multiple point source model
-[x] Quick centroid moment tensor inversion
    A nice parallelized CMT inversion based on your pre-generated Green's function. 
    The code allows GFs recycle and station inconsistency, which will search the available stations in the pre-built GFs for the inversion.
-[ ] Multiple point source inversion (In Prep)
```
****
## 1. Installation
#### cd to the place where you want to put the source code  
```console
cd Your_Local_Path  
git clone https://github.com/jiunting/GPSCMT.git
```
#### Add GPSCMT to PYTHONPATH

> Go to your environval variable file (.base_profile or .bashrc)  
```console
vi ~/.bashrc  
```
> or  
```console
vi ~/.bash_profile      
```
> and add the following line in the file

```bash
#set GPSCMT
export PYTHONPATH=$PYTHONPATH:YOUR_PATH_MARGE/GPSCMT/src/python
```    

#### Add fk to env
```bash
#fk package
export PATH=/usr/local/fk:$PATH
```   
> Note that the fk.pl calls fk defined by environtal variable, make sure fk work in any path


****
## 2. Forward & Inversion
#### Example code for forward calculation is provided in:  
```GPSCMT/example/Forward/Forwardtest.GPSCMT.py``` 
#### Example code for forward calculation is provided in:  
```GPSCMT/example/Forward/Nantou0602/Nantou0602.GPSCMT.py```

>Example work of the ```GPSCMT```  
![][fig1]
![][fig2]



You can also build finite fault inversion model by the GPSCMT, but I would recommand [Mudpy][Mudpy] which is a well-written package.

Add [Issues][Issue_lnk] if you have questions, ideas, or would like to contribute to the code via [Step-by-Step Fork tutorial][Fork_lnk].  
or simply email: ```jiunting AT uoregon DOT edu```

[Mudpy]:https://github.com/dmelgarm/MudPy "Diego's Mudpy link"
[Issue_lnk]:https://github.com/jiunting/GPSCMT/issues
[Fork_lnk]:https://guides.github.com/activities/forking/
[fig1]:./figs/GPSCMT_expfig1.png
[fig2]:./figs/GPSCMT_expfig2.png
