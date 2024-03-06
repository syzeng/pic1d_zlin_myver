import numpy as np  
import matplotlib.pyplot as plt  
  
# 假设每个数组有 m 个元素  
m = 500  # 替换为您的数组的实际大小  
n = 4   # 您想要生成的子图数量    
# 读取文件内容  
with open('phi_mode.out', 'r') as file:  
    data = np.array([float(line.strip()) for line in file]).reshape(-1, m)  
  
# 现在 data 是一个二维数组，其中每一行是一个数组  
# 计算数组中数组的数量  
num_arrays = data.shape[0]  
  
# 确保我们有足够的数据来填充四个子图  
if num_arrays >= n:  
    # 创建带有四个子图的图形窗口  
    fig, axs = plt.subplots(2, 2)  # 2行2列的子图  
  
    # 绘制每个子图的数据  
    for i in range(n):  
        ax = axs.flatten()[i]  # 获取子图的轴对象  
        ax.plot(data[i])  
        ax.set_title(f'phi_mode {i+1}')  
        ax.set_xlabel('Index')  
        ax.set_ylabel('Value')  
  
    # 调整子图之间的间距  
    plt.tight_layout()  
  
    # 显示图形窗口  
    plt.show()  
else:  
    print("Not enough arrays to fill 4 subplots.")


# # 绘制每个数组  
# for i in range(num_arrays):  
#     plt.figure(figsize=(10, 6))  # 创建新的图形窗口  
#     plt.plot(data[i])  # 绘制当前数组  
#     plt.title(f'Array {i+1}')  # 设置图表标题  
#     plt.xlabel('Index')  # 设置x轴标签  
#     plt.ylabel('Value')  # 设置y轴标签  
#     plt.show()  # 显示图表  
  
# 如果您希望将所有数组绘制在同一个图上，可以使用以下代码：  
plt.figure(figsize=(10, 6))  
for i in range(num_arrays):  
    plt.plot(data[i], label=f'Array {i+1}')  
plt.legend()  # 显示图例  
plt.xlabel('Index')  
plt.ylabel('Value')  
plt.title('All modes Together')  
plt.show()