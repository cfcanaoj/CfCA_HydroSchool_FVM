# 3D visualization

[Go to top](../README.md)  

This document is used for CfCA hydro school 2022. 

## Make HDF files and xdmf text

To learn about VisIt and paraview, you need to make the data for it. First login the analysis server.

    ssh <your account>@an**.cfca.nao.ac.jp
    
To treat HDF5 format, you need library. See [the instrucion](./InstallHDF5.md).
If you just want to play with it, add the follwoing command in your`~/.bashrc`.
    
    module load intel
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/hydro16/hdf5/lib
    module load visit
    
Then copy the source code if you have not copy it yet.

    cd /cfca-work/<your account>
    cp -r /cfca-work/hydro16/HydroSchool2025/ .
To run the code, you need to compile `main.f90`.
    
    cd HydroSchool2022/vis3D
    make makedata.x
    
Then `makedata.x`is made in this directory.
    
    ./makedata.x
    
You can find the data in `hdfdata`.

## details of xdmf
See the links.
- https://qiita.com/hsimyu/items/6e163428477d19429576
- https://www.visitusers.org/index.php?title=Using_XDMF_to_read_HDF5

## Visualization

To use VisIt, go to the [instruction](https://moored-cave-326.notion.site/VisIt-535a672b1e024c018c04e38015d2a249).

## Preparation
Abobe, you used a library in `~/hydro17/hdf5`. After the school you cannot use it. You need to install `hdf5` in your home directory in analysis server. 

