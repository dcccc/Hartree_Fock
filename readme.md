## 闭壳的hatree fock和DFT(rhf/rks)小程序
使用python和C++混合写成的闭壳HF和DFT计算的代码，主要目的在于加深对量子化学理论的理解。

### 使用方法

程序中的C++部分在运行前需要编译成动态链接库，以方便在后续运行中被python调用
其中grid_int.cpp，analy_int.cpp，Lebedev-Laikov.c编译方法如下

'''
g++ grid_int.cpp -fPIC -shared -O3 -o grid_int.so

# analy_int.cpp 中使用了boost中的函数
g++ analy_int.cpp -fPIC -shared -O3 -o analy_int.so -IC:\pro_Program\boost_1_64_0 

# Lebedev-Laikov.c 引用于 https://github.com/Rufflewind/lebedev_laikov
gcc Lebedev-Laikov.c -fPIC -shared -O3 -o liblebedevlaikov.so 
'''

此外运行还需numpy和scipy

之后便可运行

这里已经整理好了一份win10上的，下载地址为 链接: https://pan.baidu.com/s/1fPsyLRKrgTTD9oQs1DJjDQ 提取码: xdpq

### 输入文件 

输入文件第一行为关键词行，仅支持几个个关键词，分别为表示计算类型的“dft”，“hf”，“mp2”，表示基组类型的 “sto3g”，“3-21g”，“6-31g”，和表示使用收缩基计算的“contrcted”，其它词将被忽略

目前仅支持元素周期表前两周期的元素，并且由于代码计算速度较慢，请不要尝试大体系的计算

输入文件见示例如下

'''
hf  contrcted  sto3g                                       #关键词行
                                                           #空行
H      -5.566915644    0.795031347    0.111482211          #原子坐标
O      -4.858802109    0.112000847    0.625433498 
H      -4.035840572    0.704027786    1.077471763 
H      -2.150596855    1.132059269   -0.159430515 
O      -1.442483321    0.449028769    0.354520772 
H      -0.619521783    1.041055708    0.806559037 
'''

直接运行
'''
python scf.py h2o.input
'''
计算输出将打印在屏幕屏幕上



### 存在的问题

1. 如果输入结构离平衡结构较远，可能scf计算中会发生震荡而一直不收敛
2. DFT(泛函为最简单的xα)计算得到的能量与orca计算的能量之间存在约万分之一的相对差别，计算能量绝对差别会随着体系的增大而等比例增大，该bug出现的原因暂时未知
3. 其它未知的bug
