import numpy as np
import matplotlib.pyplot as plt

# 定义函数
def func_plus(x):
    return x**2 + 0.5 + 0.5 * np.sqrt(8*x**2 + 1)

def func_minus(x):
    return x**2 + 0.5 - 0.5 * np.sqrt(8*x**2 + 1)

# 创建x轴坐标范围
x = np.linspace(0, 1, 200)  # 根据需要调整范围和点数

# 计算对应的y值
y_plus = func_plus(x)
#y_minus = func_minus(x)
#y_minus = abs(func_minus(x))**0.5   #gamma
y_minus = 2*abs(func_minus(x))**0.5   #gamma energy

# 绘制图形
plt.figure(figsize=(12, 8))  # 设置图形大小
#plt.plot(x, y_plus, label=r'$y = x^2 + \frac{1}{2} + \frac{1}{2}\sqrt{8x^2+1}$', color='blue')
plt.plot(x*2.5, y_minus, label=r'$\gamma_{Eenergy} $', color='red')

plt.xlabel(r'$k_m V_b$')
plt.ylabel(r'$\gamma_{Eenergy} /w_{pe}$')
plt.title(r'$\gamma_{Eenergy} -k_m V_b$')
plt.legend()
plt.grid(True)

# # 指定文件名  
# file_mfit = 'mfit.txt'  
# file_v = 'v_list.txt'  
# # 读取数据文件  
# with open(file_mfit, 'r') as fm, open(file_v, 'r') as fv:  
#     arm = np.loadtxt(fm)  
#     arv = np.loadtxt(fv)  

# print(arm,arv)
arm=np.array([0.03375884, 0.05017035, 0.06108921, 0.07220455, 0.07428468, 0.07236756, 0.06139279, 0.05395995])
arv=np.array([0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.2])
plt.scatter(arv,arm*10,s=10, color='blue')

plt.show()