# pic1d_zlin_myver
a version to push
基于pic1d_zlin.f90开发的系列程序
## 文件介绍
### documents 
说明文档
主要说明文档在 documents/PIC1D调试说明文档.md （以下简称《说明文档》）里。
其余文件为辅助文件、中间文件等。

### maxwell_particle_evolution 
Maxwell分布下电子相空间演化
fortran程序：pic1d_diy_2_full_more_explain.f90 
python绘制粒子云图程序：py_readparticle_many_picture.py
python动图绘制程序：py_imageio_combine2gif.py

### twostream_particle_cloud
双流初始条件下的演化，绘制云图。
fortran程序：pic1d_diy_2stream_fullf.f90
此fortran程序主要是按照《说明文档》中的方法对load模块进行修改，使电子初始速度分布为双流分布。涉及的新参数为nbeam，vbeam。
python绘制粒子云图程序：py_readparticle_many_picture.py
python动图绘制程序：py_imageio_combine2gif.py

### twostream_f_grid
最初的双流演化测试程序，采用对粒子加权到xv相空间网格的方法进行诊断。实质上是测试了对 movie.out 诊断数据进行绘图处理。
尽管平均到网格的做法过于粗糙，但可以直接绘制分布函数。

### twostream_growth
考察增长率。
fortran程序：pic1d_diy_2stream_fullf.f90
此fortran程序大致和twostream_particle_cloud 中的一样，但输出端做了更改，输出电场能量随时演化文件"Ef_energy.txt"。同时滤波器默认选取只通基频。

py_readEfieldEnergy_plot.py
读取Ef_energy.txt,绘制电场能量随时演化对数图，方便目测选取拟合区域。

py_fieldEnergy_fit.py
由py_readEfieldEnergy_plot.py更改而来，读取Ef_energy.txt，对数处理，选段线性拟合，输出拟合效果图和拟合参数（斜率、斜率的不确定度）。

数据保存：kv-gamma.xlsx
在拟合得到增长率后，将结果复制到excel里。

#### 理论解程序
pyplot_yx.py 绘制y-x曲线
pyplot_gamma.py 绘制gammaE-kv图
pyplot_energy_gamma.py 绘制能量增长率理论解和数值拟合结果，对比。 

### results 运行结果和中间文件
运行结果保存在results中，由于中间文件较多较大，在git中被忽略。

## 当前问题
### 文件管理：基于git解决
遇到比较大的问题是文件管理问题。有时候git在push时会有bug，通过换电脑解决。中间结果一大堆，不上传github，可用u盘保存results。
需要确定开发流程，有新文件后怎么处理，这种问题。
最理想情况：所有设备的git都没有问题，git也不用考虑文件大小。

## 下步工作
原始的控制方程
离散方程、离散方法及其修正方程
时间步和真实时间
各子程序计算流程图绘制

主要：增长率与解析解的对比
需要推导增长率解析式
1.基本：从国际（或高斯）单位制公式到模拟归一化单位、公式。参数选取。
2.