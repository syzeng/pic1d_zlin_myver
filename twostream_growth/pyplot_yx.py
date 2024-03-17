import numpy as np
import matplotlib.pyplot as plt

# 定义函数
def func_plus(x):
    return x**2 + 0.5 + 0.5 * np.sqrt(8*x**2 + 1)

def func_minus(x):
    return x**2 + 0.5 - 0.5 * np.sqrt(8*x**2 + 1)

# 创建x轴坐标范围
x = np.linspace(-2, 2, 400)  # 根据需要调整范围和点数

# 计算对应的y值
y_plus = func_plus(x)
y_minus = func_minus(x)
#y_minus = abs(func_minus(x))**0.5

# 绘制图形
plt.figure(figsize=(8, 6))  # 设置图形大小
#plt.plot(x, y_plus, label=r'$y = x^2 + \frac{1}{2} + \frac{1}{2}\sqrt{8x^2+1}$', color='blue')
plt.plot(x, y_minus, label=r'$y = x^2 + \frac{1}{2} - \frac{1}{2}\sqrt{8x^2+1}$', color='red')

plt.xlabel('x')
plt.ylabel('y')
plt.title('Plot of the given functions')
plt.legend()
plt.grid(True)
plt.show()