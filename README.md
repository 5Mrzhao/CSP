# CSP
Comprehensively-Sparsified Preconditioner for Efficient Nonlinear Circuit Simulation

#### What is CSP(v1.0) ？
CSP transforms nonlinear components into symmetric Laplacian matrices, seamlessly blending nonlinear and linear elements into the sparsification process. By intersecting this refined sparsifier with the original Modified Nodal Analysis (MNA) matrix, CSP further trims down the sparsity. This clever approach not only sharpens the accuracy but also speeds up the factorization time of the preconditioner, making the overall process more efficient. 

#### Need to do !
Download gcc and g++:
```c++
sudo apt update
sudo apt install gcc
gcc --version
sudo apt-get install g++
```
Download CSP:
```c++
unzip CSP.zip
```
Begin testing:
```c++
mkdir saveMtx
sh run.sh
./main matrix_name
sh test.sh
```
Then Generate three types of output files:

(a) **name.csv**: Stored in the saveMtx directory, this file contains the name of the successfully solved matrix.

(b) **""/"_gpscp/"_csp"**: Also located in the saveMtx directory, this file "_CSP" includes the pre-conditioned matrix after CSP preprocessing.

(c) **getdata.csv**: Found in the root directory, this file records statistics on the CSP runtime.

#### What is CSP(v1.1) ？
CSP enhances the parallelization of the spectral sparsification strategy by incorporating block RMQ and point exclusivity algorithms, significantly accelerating preprocessing.

..........
