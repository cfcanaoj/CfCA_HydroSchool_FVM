# Installing HDF5 with Fortran Support (User Installation)

This guide explains how to install **HDF5 1.14.3** in a **user directory (without root privileges)** on a Linux system.  
The installation enables the **Fortran interface**.

The library will be installed in: $HOME/hdf5


---

# 1. Create directory

Create directories for source code and installation.

```bash
mkdir -p $HOME/hdf5
```

# 2. Download HDF5 source code

Download the stable release of HDF5 1.14.3.
```bash
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-1.14.3/src/hdf5-1.14.3.tar.gz
```
Extract the archive
```bash
tar xvf hdf5-1.14.3.tar.gz
```

# 3. Configure the build
Here if you want to use different compiler, please change FC, e.g., ifort.
```bash
cd hdf5-1.14.3
FC=gfortran ./configure --prefix=$HOME/hdf5 --enable-fortran
```
The -j option allows parallel compilation.
```bash
make -j
make install
```
After installation, the directory structure will look like:
```bash
$HOME/hdf5
 ├── bin
 ├── include
 ├── lib
 └── share
```

#4 Update environment variables
```bash
export PATH=$HOME/hdf5/bin:$PATH
export LD_LIBRARY_PATH=$HOME/hdf5/lib:$LD_LIBRARY_PATH
export HDF5_DIR=$HOME/hdf5
```
Apply the chgenge

```bash
source ~/.bashrc
```

#5. Verify the installation

Check the HDF5 configuration:
```bash
h5cc -showconfig
```
If the output includes the following text, the Fortran interface has been successfully installed.
```bash
Fortran: yes
```

Check the Fortran wrapper:
```bash
h5fc -show
```

# 6. Compiling a Fortran program with HDF5

The easiest way is to use the HDF5 compiler wrapper:
```bash
h5fc program.f90
```
Alternatively, compile manually:
```bash
gfortran program.f90 \
  -I$HOME/local/include \
  -L$HOME/local/lib \
  -lhdf5_fortran -lhdf5
```
