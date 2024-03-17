import numpy as np  
import matplotlib.pyplot as plt  
import os  

nparticle = 60000  
vbeam=1
ntime = 400  # 假设粒子位置数组的长度代表总的时间步数  
ndiag = 2  # 假设每ndiag个时间步进行一次诊断  
mystep=1 #每隔mystep个点选一个点

fileE = 'Ef_energy.txt' 
E_field_energy=np.loadtxt(fileE)
time_points=np.arange(1,ntime,ndiag)

# 使用matplotlib绘制图像  
plt.figure(figsize=(12, 8), dpi=100)  
plt.plot(time_points, E_field_energy, marker='o', markersize=1)  # marker='o' 表示在每个数据点上画一个圆圈  
plt.xlabel('Time(0.1/w_pe)')  # 设置横轴标签  
plt.ylabel('E_field_energy')  # 设置纵轴标签  
plt.title('E_field_energy over Time,vbeam={:}'.format(vbeam))  # 设置图像标题  
plt.grid(True)  # 显示网格  

# 设置纵轴为对数坐标  
plt.yscale('log') 

plt.show()  # 显示图像