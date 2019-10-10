# GPSCMT

This is the forward and inverse code for 1-D point source used in the study: 
Lin, J. T., Chang, W. L., Melgar, D., Thomas, A., & Chiu, C. Y. (2019). Quick determination of earthquake source parameters from GPS measurements: a study of suitability for Taiwan. Geophysical Journal International, 219(2), 1148-1162.

What it can do:
    Forward model:
        One source-to-single/many stations
            e.g. You have a point source, what are the surface displacements look like?
        Many different sources-to-single/many stations
            e.g. You have many sources, what are the corresponding surface displacements?
        Multiple point source model
        Generate Green's function
    Inversion:
        Quick centroid moment tensor inversion
            A nice parallelized CMT inversion based on your pre-generated Green's function. The code allows GFs recycle and station inconsistency, which will search the available stations in the pre-built GFs for the inversion.
        Multiple point source inversion (working...)

    
You can also build finite fault inversion model by the GPSCMT, but I would recommand my advisor Diego Melgar's Mudpy https://github.com/dmelgarm/MudPy which is a very neat package.
