import numpy as np  
import matplotlib.pyplot as plt  
from scipy.optimize import curve_fit  
  
# 假设你已经有了一些数据点  
time_points = np.linspace(1, 100, 100)  
E_field_energy = time_points ** 2 + 50 * time_points + 100 + np.random.normal(0, 50, 100)  # 假设的非线性关系加上噪声  
  
# 对数据进行对数转换  
log_time = np.log10(time_points)  
log_energy = np.log10(E_field_energy)  
  
# 定义对数转换后的线性拟合函数  
def log_linear_fit(x, m, c):  
    return m * x + c  
  
# 使用curve_fit进行拟合  
popt, pcov = curve_fit(log_linear_fit, log_time, log_energy)  
m_fit, c_fit = popt  
  
# 绘制原始数据的对数转换结果和拟合线  
plt.plot(log_time, log_energy, marker='o', label='Log-transformed Data')  
plt.plot(log_time, log_linear_fit(log_time, *popt), 'r', label='Fitted Line on Log-transformed Data')  
  
# 设置轴标签和标题  
plt.xlabel('Log(Time)')  
plt.ylabel('Log(E_field_energy)')  
plt.title('Log-transformed E_field_energy over Time with Linear Fit')  
  
# 显示图例  
plt.legend()  
  
# 显示网格  
plt.grid(True)  
  
# 显示图像  
plt.show()