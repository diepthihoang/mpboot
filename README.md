# mpboot
Asynchronous version for MPBoot

## Downloading source code
You can clone the source code from GitHub with:

`git clone https://github.com/diepthihoang/mpboot.git`  
Switch to branch which contains asynchronous version of MPBoot-MPI source code:
* `git checkout mpboot-mpi-async`

## Compiling under Linux
1. Create folder **build** outside folder **mpboot**.
2. Open a Terminal.
3. Change directory to **build**
4. Configure source code with CMake:  
`cmake ../source -DIQTREE_FLAGS=avx -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx`
> Replace avx by sse4 in above command if you decide to run MPBoot on SSE architecture

5. Run command `make -j4` to compile source code with 4 processes:  
> Option **j** specifies the number of processes used to compile source code with **make**.
  
The compiler will generate an executable file named **mpboot-avx**
> In case of running MPBoot on SSE architecture, the executable file is named **mpboot**.

6. To analyst file **example.phy** with 2 processes, run command:  
`mpirun -np 2 ./mpboot-avx -s example.phy`
> Option **np** specifies the number of processes used to run MPBoot-MPI

## Compiling under Mac OS X
1. Create folder **build** outside folder **mpboot**.
2. Open a Terminal.
3. Change directory to **build**
4. Configure source code with CMake:  
`cmake ../mpboot -DIQTREE_FLAGS=avx -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx`
> Replace avx by sse4 in above command if you decide to run MPBoot on SSE architecture

5. Run command `make -j4` to compile source code with 4 processes. Option **j** specifies the number of processes used to compile source code with **make**.
  
The compiler will generate an executable file named **mpboot-avx**
> In case of running MPBoot on SSE architecture, the executable file is named **mpboot**.

6. To analyst file **example.phy** with 2 processes, run command:  
`mpirun -np 2 ./mpboot-avx -s example.phy` 
> Option **np** specifies the number of processes used to run MPBoot-MPI

## Compiling under Windows
* Requirements:  
  * cmake version >= 3.21
  * TDM-GCC
  * MSMPI
1. Create folder **build** outside folder **mpboot**.
2. Open a Terminal.
3. Change directory to **build**
4. Configure source code with CMake:  
`cmake -G "MinGW Makefiles" -DIQTREE_FLAGS=mpiavx ../mpboot`
> Replace mpiavx by mpisse4 in above command if you decide to run MPBoot on SSE architecture.  
> Due to having conflicts with **Vectorization**, please not using **Clang** to configure source code.

5. Run command `mingw32-make -j4` to compile source code with 4 processes:  
> Option **j** specifies the number of processes used to compile source code with **make**.
  
The compiler will generate an executable file named **mpboot-avx** 
> In case of running MPBoot on SSE architecture, the executable file is named **mpboot**.
  

6. To analyst file **example.phy** with 2 processes, run command:  
`mpiexec -n 2 ./mpboot-avx -s example.phy`
> Option **n** specifies the number of processes used to run MPBoot-MPI
