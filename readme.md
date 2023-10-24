## Closed shell Hartree-Fock and DFT (rhf/rks) calculation code

It is a mixture of python and C++ code. The purpose of this code is to gain a better understanding about the basic theory of quantum chemistry . 

### How To build

C++ code parts have to be compiled to be dynamic lib, so the python code can use functions in it. Of course, a pure python verison of funtions also prepared, which will be used when c++ libs are not presented. But the python verison functions are slower than c++ funtions.
The grid_int.cpp，analy_int.cpp and Lebedev-Laikov.c should be compiled like

```bash
g++ grid_int.cpp -fPIC -shared -O3 -o grid_int.so

# tgamma and gamma_p functions from boost is used in analy_int.cpp
g++ analy_int.cpp -fPIC -shared -O3 -o analy_int.so -IC:\pro_Program\boost_1_64_0 

# Lebedev-Laikov.c comes from https://github.com/Rufflewind/lebedev_laikov
gcc Lebedev-Laikov.c -fPIC -shared -O3 -o liblebedevlaikov.so 
```

And packages numpy and scipy are also needed

The pre-build file for windows is ready, which is [HF-DFT_win.zip](https://github.com/dcccc/Hartree_Fock/blob/master/HF-DFT_win.zip).

### Input File

The first line is the keyword line, only several keywords are recognized. They are 'dft', 'hf', 'mp2' for method, 'sto3g'，'3-21g' and '6-31g' to set the basis set, and 'contracted' to use the contracted basis set, or original basis set will be used. All other word will be omitted.

Only elements from H to Ne are supported. As calculation is not fast, so don't attempt large systems. And mp2 method always takes much more time than HF, so use it with caution.

A sample input file for 2 H2O is like 

```
hf  contracted  sto3g                                       
                                                           
H      -5.566915644    0.795031347    0.111482211          
O      -4.858802109    0.112000847    0.625433498 
H      -4.035840572    0.704027786    1.077471763 
H      -2.150596855    1.132059269   -0.159430515 
O      -1.442483321    0.449028769    0.354520772 
H      -0.619521783    1.041055708    0.806559037 
```

run 

```bash
python scf.py h2o.input
```

the output line will be printed on the screen



### Bugs


1. If the structure is far from stable conformer, the scf may will not converge

2. As the numerical integration in DFT is not accurate enough (functional is xα) , so the energy results will be a little different from the values by orca(using sto3g contracted basis set for CH4, the difference is less than 0.0005 Hartree). The difference will be large for big structure.

3. For hf or dft method, the final energy resuslts between using pure python functions and functions in c++ libs ame input context may be different, but usually the difference is small, and can assumed to be numerial error. 

3. mp2 calculation can only conducted when grid_int.cpp is compiled to be grid_int.so properly.

4. Other unknown bugs


---------------------------------------------------------------------


## 闭壳的Hatree-Fock和DFT(rhf/rks)小程序

使用python和C++混合写成的闭壳HF和DFT计算的代码，主要在于加深对量子化学理论的理解

### 使用方法


程序中的C++部分在运行前需要编译成动态链接库，以方便在后续运行中被python调用
如果没有编译c++部分，计算将使用内置的纯python函数。c++库版本的函数总是比纯python的函数快，特别是dft方法
其中grid_int.cpp，analy_int.cpp，Lebedev-Laikov.c编译方法如下


```bash
g++ grid_int.cpp -fPIC -shared -O3 -o grid_int.so

# analy_int.cpp 中使用了boost中的函数
g++ analy_int.cpp -fPIC -shared -O3 -o analy_int.so -IC:\pro_Program\boost_1_64_0 

# Lebedev-Laikov.c 引用于 https://github.com/Rufflewind/lebedev_laikov
gcc Lebedev-Laikov.c -fPIC -shared -O3 -o liblebedevlaikov.so 
```

此外运行还需numpy和scipy


之后便可运行

已经编译好可以在windows上运行的文件，如HF-DFT_win.zip




### 输入文件 


输入文件第一行为关键词行，仅支持几个关键词，分别为表示计算类型的“dft”，“hf”，“mp2”，表示基组类型的 “sto3g”，“3-21g”，“6-31g”，和表示使用收缩基计算的“contracted”. 其它词将被忽略

目前仅支持元素周期表前两周期的元素，并且由于计算速度较慢，请不要尝试大体系的计算. mp2的计算将花费比HF多好几倍的时间，请谨慎使用


输入文件见示例如下

```
hf  contrcted  sto3g                                       
                                                           
H      -5.566915644    0.795031347    0.111482211          
O      -4.858802109    0.112000847    0.625433498 
H      -4.035840572    0.704027786    1.077471763 
H      -2.150596855    1.132059269   -0.159430515 
O      -1.442483321    0.449028769    0.354520772 
H      -0.619521783    1.041055708    0.806559037 
```

直接运行

```
python scf.py h2o.input
```

计算输出将打印在屏幕上




### 存在的问题


1. 如果输入结构离平衡结构较远，可能scf计算中会发生震荡而一直不收敛

2. 由于格点积分精度不够高，DFT(泛函为最简单的xα)计算得到的能量与orca计算的能量之间存在一些差别(ch4使用sto3g收缩基与ORCA差别小于0.0005Hartree)，较大的体系可能差别较大而不可忽略

3. 对于hf和dft，c++库中的函数计算的能量结果和纯python函数计算的结果之间存在一定的差别，但一般而言差别非常小，可以认为只是数值计算误差

4. mp2的计算只能在grid_int.cpp被编译成grid_int.so库文件时才能计算，由于mp2计算中有一步非常耗时，所以没有纯python版本的函数可以计算

5. 其它未知的bug


### 主要参考资料


1. 由网友zyniso书写整理的关于RHF编程的资料，包含四个网页，分别为http://bbs.keinsci.com/thread-10317-1-1.html 、 http://bbs.keinsci.com/thread-10368-1-1.html 、 http://bbs.keinsci.com/thread-10427-1-1.html 、 http://bbs.keinsci.com/thread-10445-1-1.html
2. 网友sobereva书写整理的关于DFT方法中的积分计算内容资料，http://sobereva.com/69
3. 书籍《MOLECULAR ELECTRONIC-STRUCTURE THEORY》
4. 网友Shannon编写的HF代码，http://bbs.keinsci.com/thread-933-1-1.html
5. 网友978142355编写的使用HF，mp2方法计算Hn的能量的代码，http://bbs.keinsci.com/thread-6139-1-1.html
6. 关于DFT中格点积分方法介绍的文章《A program to generate a basis set adaptive radial quadrature grid for density functional theory》，以及选取积分点的文章《A standard grid for density functional calculations》
7. github用户Rufflewind编写的用于生成球面采样点的程序，https://github.com/Rufflewind/lebedev_laikov