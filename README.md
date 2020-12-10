

<!--
 * @Author: Rui Wang
 * @Date: 2020-12-10 11:06:29
 * @LastModifiedBy: Rui Wang
 * @LastEditTime: 2020-12-10 12:17:53
 * @Email: wangru25@msu.edu
 * @FilePath: /HERMERS/README.md
 * @Description: 
-->
# HERMERS

Highly Efficient Robust Multidimensional Evolutionary Spectra (HERMERS) is a library for evaluating both the harmonic and non-harmonic spectra of persistent Laplacian operators. In the present realse, we consider an implementation in alpha complexes.  

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
cd build
./HERMERS InputData Filtration Num
```
- InputData: The point cloud data is allowed
- Filtration: The radius filtration parameters 
- Num: The number of eigenvalues that will be calculated
### For example:
```bash
cd build
./HERMERS Test_C60.xyz filtration 100
```

## Documentation 

Documentation for HERMERS can be found [here](https://users.math.msu.edu/users/weig/HERMES).

## Contributors

HERMERS was created by [Rundong Zhao](https://github.com/rdzhao) and is maintained by [Yiying Tong](xxx), [Emily Ribando-Gros](xxx), [Jiahui Chen](https://github.com/Jiahuic), [Rui Wang](https://github.com/wangru25), and [Weilab at MSU Math](https://github.com/msuweilab).

