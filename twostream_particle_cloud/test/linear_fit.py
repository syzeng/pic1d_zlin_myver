import numpy as np  
import matplotlib.pyplot as plt  
from scipy.optimize import curve_fit  
  
# 假设你已经有了一些数据点  
time_points = np.linspace(0, 10, 100)  
E_field_energy = 3 * time_points + 2 + np.random.normal(0, 1, 100)  # 线性关系加上一些噪声  
  
# 选择你想要拟合的数据段，例如从第20个点到第80个点  
start_idx = 20  
end_idx = 80  
time_fit = time_points[start_idx:end_idx]  
E_fit = E_field_energy[start_idx:end_idx]  
  
# 定义线性拟合函数  
def linear_fit(x, m, c):  
    return m * x + c  
  
# 使用curve_fit进行拟合  
popt, pcov = curve_fit(linear_fit, time_fit, E_fit)  
m_fit, c_fit = popt  
  
# 绘制原始数据  
plt.plot(time_points, E_field_energy, marker='o', label='Original Data')  
  
# 绘制拟合线，仅绘制你进行拟合的数据段范围内的线  
plt.plot(time_fit, linear_fit(time_fit, *popt), 'r', label='Fitted Line')  
  
# 设置轴标签和标题  
plt.xlabel('Time')  
plt.ylabel('E_field_energy')  
plt.title('E_field_energy over Time with Linear Fit')  
  
# 显示图例  
plt.legend()  
  
# 显示网格  
plt.grid(True)  
  
# 显示图像  
plt.show()