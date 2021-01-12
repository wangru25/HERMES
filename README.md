

<!--
 * @Author: Rui Wang
 * @Date: 2020-12-10 11:06:29
 * @LastModifiedBy: Rui Wang
 * @LastEditTime: 2021-01-12 17:02:18
 * @Email: wangru25@msu.edu
 * @FilePath: /HERMES/README.md
 * @Description: 
-->
# HERMES

HERMES is a software package for simultaneous topological data analysis (persistent Betti numbers) and geometric data analysis (persistent eigenvalues). It is realized through persistent spectral graph theory. In the present release, we consider an implementation in alpha complexes.

## Requirements
- cmake 3.1 or higher
- gcc 7.5.0
- GNU Make 4.1
- MATLAB
- [CGAL](https://www.cgal.org/) version 4.14

## Install and Build
### You can intall HERMES directly from the development repository:
```bash
git clone https://github.com/wangru25/HERMES.git
```

### How to build the project:
```bash
mkdir build
cd build
cmake ..
make
```

### Notice:
Please make sure the MATLAB directory in the [CMakeLists.txt](https://github.com/wangru25/HERMES/blob/main/CMakeLists.txt) (Line 15 and Line 20) matches with yours. 


## Examples
There are several examples projects in the [examples](https://github.com/wangru25/HERMES/tree/main/examples).
### How to run
```bash
./Snapshot InputData Filtration Num P
```
- InputData: The point cloud data is allowed
- Filtration: The radius filtration parameters 
- Num: The number of eigenvalues that will be calculated
- P: The persistent value
### For example:
```bash
cd examples
./../build/Snapshot Test_C60.xyz filtration.txt 100 0.4
```
- The spectra of the 0th-order 0.4-persistent Laplacian will be saved in [examples/snapshots_vertex.txt](https://github.com/wangru25/HERMES/blob/main/examples/snapshots_vertex.txt). Each line presents harmonic or non-harmonic eigenvalues at a specific filtration value.
- The spectra of the 1th-order 0.4-persistent will be saved in [examples/snapshots_edge.txt](https://github.com/wangru25/HERMES/blob/main/examples/snapshots_edge.txt). Each line presents harmonic or non-harmonic eigenvalues at a specific filtration value.
- The spectra of the 2th-order 0.4-persistent will be saved in [examples/snapshots_facet.txt](https://github.com/wangru25/HERMES/blob/main/examples/snapshots_facet.txt). Each line presents harmonic or non-harmonic eigenvalues at a specific filtration value.


## Documentation 

Documentation for HERMES can be found [here](https://weilab.math.msu.edu/HERMES/).

## Citing
You may use the following bibtex entry to cite HERMES:
```
@article{wang2020hermes,
  title={HERMES: Persistent spectral graph software},
  author={Wang, Rui and Zhao, Rundong and Ribando-Gros, Emily and Chen, Jiahui and Tong, Yiying and Wei, Guo-Wei},
  journal={arXiv preprint arXiv:2012.11065},
  year={2020}
}
```

## References
- R. Wang, R. Zhao, E. Ribando-Gros, J. Chen, Y. Tong, and G.-W. Wei. [Hermes: Persistent spectral graph software](https://arxiv.org/pdf/2012.11065.pdf), 2020.
- R. Wang, D. D. Nguyen, and G.-W. Wei. [Persistent spectral graph](https://users.math.msu.edu/users/weig/paper/p243.pdf). International Journal for Numerical Methods in Biomedical Engineering, page e3376, 2020.


## Contributors

HERMES was created by [Rundong Zhao](https://github.com/rdzhao) and is maintained by [Yiying Tong](xxx), [Emily Ribando-Gros](https://github.com/eribandogros), [Jiahui Chen](https://github.com/Jiahuic), [Rui Wang](https://github.com/wangru25), and [Weilab at MSU Math](https://github.com/msuweilab).

