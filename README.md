# rhoLDFSSFoam_gpu

This is **gpu version** of the **rhoLDFSSFoam solver** from [rhoLDFSS](https://github.com/gauravkumar463/rhoLDFSSFoam.git). This solver can be compiled and run after installing  [RapidCFD-dev](https://github.com/SimFlowCFD/RapidCFD-dev.git). 

You can compile RapidCFD-dev in a container using cuda8.def and rapidcfd.def files from apptainer directory. First build cuda8.sif which contains essential packages needed by RapidCFD-dev (or OpenFOAM) such as nvidia-drivers, CUDA-8.1 etc. Then build RapidCFD-dev using cuda.sif. 

```
apptainer build cuda.sif cuda.def

apptainer build rapidcfd.sif rapidcfd.def
```

You'll need to download the **ThirdParty-dev** package from OpenFOAM's website for parallel operations using OPENMPI and SCOTCH.

Once RapidCFD-dev is built without error, rhoLDFSSFoam can be compile using the following command:

``apptainer run --nv ~/bin/rapidcfd.sif sh -c 'export WM_PRECISION_OPTION=SP; . /usr/local/RapidCFD/RapidCFD-dev/etc/bashrc -f; wmake'``

You should replace ``WM_PRECISION_OPTION=SP`` to ``WM_PRECISION_OPTION=DP`` if you want to do double precision calculations.

solver can be run using the following command:

``apptainer run --nv ~/bin/rapidcfd.sif sh -c 'export WM_PRECISION_OPTION=SP; . /usr/local/RapidCFD/RapidCFD-dev/etc/bashrc -f; rhoLDFSSFoam'``

Note: I had to change ``-arch=sm_30`` to ``-arch=sm_60`` in *RapidCFD-dev/wmake/rules/linux64Nvcc/c* and *RapidCFD-dev/wmake/rules/linux64Nvcc/c++* files for NVIDIA GeForce RTX 4090 GPU card.

Packages required to install cuda-8.1 in Ubuntu-16.04 container:
```
nvidia-diag-driver-local-repo-ubuntu1604_375.66-1_amd64.deb

cuda-repo-ubuntu1604-8-0-local-ga2_8.0.61-1_amd64.deb

cuda-repo-ubuntu1604-8-0-local-cublas-performance-update_8.0.61-1_amd64.deb 
```
