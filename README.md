

<!--
 * @Author: Rui Wang
 * @Date: 2020-12-10 11:06:29
 * @LastModifiedBy: Rui Wang
 * @LastEditTime: 2020-12-10 15:27:16
 * @Email: wangru25@msu.edu
 * @FilePath: /HERMERS/README.md
 * @Description: 
-->
# HERMERS

Highly Efficient Robust Multidimensional Evolutionary Spectra (HERMERS) is a library for evaluating both the harmonic and non-harmonic spectra of persistent Laplacian operators. In the present release, we consider an implementation in alpha complexes.  The theoretical literature can be found in [PSG](https://users.math.msu.edu/users/weig/paper/p243.pdf) and [HERMERS](xxx).

## Requirements
- cmake 3.1 or higher
- gcc 7.5.0
- GNU Make 4.1
- [CGAL](https://www.cgal.org/) version 4.14

## Install and Build
### You can intall HERMERS directly from the development repository:
```bash
git clone https://github.com/wangru25/HERMERS.git
```

### How to build the project:
```bash
mkdir build
cd build
cmake ..
make
```

### Notice:
Please make sure the MATLAB directory in the [CMakeLists.txt](https://github.com/wangru25/HERMERS/blob/main/CMakeLists.txt) (Line 15 and Line 20) matches with yours. 


## Examples
There are several examples projects in the [examples](https://github.com/wangru25/HERMERS/tree/main/examples).
### How to run
```bash
./HERMERS InputData Filtration Num P
```
- InputData: The point cloud data is allowed
- Filtration: The radius filtration parameters 
- Num: The number of eigenvalues that will be calculated
- P: The persistent value
### For example:
```bash
cd examples
./../build/HERMERS Test_C60.xyz filtration.txt 100 0.4
```
- The spectra of the *p*-persistent *0*-th combinatorial Laplacian matrix will be saved in [examples/snapshots_vertex.txt](https://github.com/wangru25/HERMERS/blob/main/examples/snapshots_vertex.txt)
- The spectra of the *p*-persistent *1*-th combinatorial Laplacian matrix will be saved in [examples/snapshots_edge.txt](https://github.com/wangru25/HERMERS/blob/main/examples/snapshots_edge.txt)
- The spectra of the *p*-persistent *2*-th combinatorial Laplacian matrix will be saved in [examples/snapshots_facet.txt](https://github.com/wangru25/HERMERS/blob/main/examples/snapshots_facet.txt)


## Documentation 

Documentation for HERMERS can be found [here](https://weilab.math.msu.edu/HERMERS/).

## Citing
You may use the following bibtex entry to cite HERMERS (will be updated with complete publication data):
```
@article{1908.02518,
	author = {xxxx},
	title = {HERMES: persistent spectral graph software},
	year = {2020},
	journal = {arXiv preprint arXiv:1804.06990}
}
```

## Contributors

HERMERS was created by [Rundong Zhao](https://github.com/rdzhao) and is maintained by [Yiying Tong](xxx), [Emily Ribando-Gros](xxx), [Jiahui Chen](https://github.com/Jiahuic), [Rui Wang](https://github.com/wangru25), and [Weilab at MSU Math](https://github.com/msuweilab).

