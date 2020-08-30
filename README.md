# GPSCMT
### The forward and inverse code for 1-D point source

* Lin, J. T., Chang, W. L., Melgar, D., Thomas, A., & Chiu, C. Y. (2019). Quick determination of earthquake source parameters from GPS measurements: a study of suitability for Taiwan. Geophysical Journal International, 219(2), 1148-1162.

What it can/cannot do
```
-[x] Forward model-One source-to-single/many stations 
    e.g. You have a point source, what are the surface displacements look like?
-[x] Many different sources-to-single/many stations
    e.g. You have many sources, what are the corresponding surface displacements?
-[x] Generate Green's function (Strike/Dip slip or moment tensors)
-[x] Multiple point source model
-[x] Quick centroid moment tensor inversion
    A nice parallelized CMT inversion based on your pre-generated Green's function. 
    The code allows GFs recycle and station inconsistency, which will search the available stations in the pre-built GFs for the inversion.
-[ ] Multiple point source inversion (In Prep)
```


You can also build finite fault inversion model by the GPSCMT, but I would recommand [Mudpy][Mudpy] which is a well-written package.

Add [Issues][Issue_lnk] if you have questions, ideas, or would like to contribute to the code via [Step-by-Step Fork tutorial][Fork_lnk].  
or simply email: ```jiunting AT uoregon DOT edu```

[Mudpy]:https://github.com/dmelgarm/MudPy "Diego's Mudpy link"
[Issue_lnk]:https://github.com/jiunting/GPSCMT/issues
[Fork_lnk]:https://guides.github.com/activities/forking/
