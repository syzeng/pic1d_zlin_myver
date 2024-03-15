import numpy as np  
import matplotlib.pyplot as plt  
from scipy.optimize import curve_fit  

nparticle = 60000  
ntime = 400  # 假设粒子位置数组的长度代表总的时间步数  
ndiag = 5  # 假设每100个时间步进行一次诊断  
mystep=1 #每隔mystep个点选一个点

fileE = 'Ef_energy.txt' 
E_field_energy=np.loadtxt(fileE)
time_points=np.arange(1,ntime,ndiag)

tsteps=len(time_points)*4//6
E_field_energy=E_field_energy[2:tsteps]
time_points=time_points[2:tsteps]

# 对数据进行对数转换 
log_energy = np.log10(E_field_energy)  

# 选择你想要拟合的数据段，先看横坐标，再计算指标
start_t=50
end_t=180

start_idx = start_t//ndiag  
end_idx = end_t//ndiag  
time_fit = time_points[start_idx:end_idx]  
E_fit = log_energy[start_idx:end_idx]  

# 定义线性拟合函数  
def linear_fit(x, m, c):  
    return m * x + c  

# 使用curve_fit进行拟合  
popt, pcov = curve_fit(linear_fit, time_fit, E_fit)  
m_fit, c_fit = popt  

# 使用matplotlib绘制图像  
plt.figure(figsize=(12, 8), dpi=100)  

# 绘制原始数据  
plt.plot(time_points,log_energy, marker='o', markersize=1,label='Log-transformed Data')  # marker='o' 表示在每个数据点上画一个圆圈  

# 绘制拟合线，仅绘制你进行拟合的数据段范围内的线  
plt.plot(time_fit, linear_fit(time_fit, *popt), 'r', label='Fitted Line on Log-transformed Data')  

plt.xlabel('Time(0.1w_pe)')  # 设置横轴标签  
plt.ylabel('Log(E_field_energy)')  # 设置纵轴标签  
plt.title('Log-transformed E_field_energy over Time with Linear Fit')  # 设置图像标题  
plt.grid(True)  # 显示网格  


plt.show()  # 显示图像