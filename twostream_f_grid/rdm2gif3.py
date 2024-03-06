import numpy as np  #用于画三维动图
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from matplotlib import cm, rcParams
def read_and_store_64x64_chunks(filename, n_chunks=10):
    matrices_list = []

    with open(filename, 'r') as file:
        lines = file.readlines()
        
        chunk_size = 64 * 64
        total_chunks = len(lines) // chunk_size
        
        # 确保文件有足够的数据读取指定次数的64x64块
        assert total_chunks >= n_chunks, f"File does not contain enough data for {n_chunks} 64x64 chunks."
        
        for _ in range(n_chunks):
            start = (_ * chunk_size)
            end = start + chunk_size
            data_chunk = [float(line.strip()) for line in lines[start:end]]

            # 将一维数据转换为64x64二维数组并添加到列表中
            matrix = np.reshape(data_chunk, (64, 64))
            matrices_list.append(matrix)

    return matrices_list

# 使用示例
nc = 450
matrices = read_and_store_64x64_chunks('movie.out', n_chunks=nc)

fig, ax = plt.subplots()
cmap = 'viridis'

def init():
    ax.set_xlim(0, 63)
    ax.set_ylim(0, 63)  # 设置坐标轴范围以匹配64x64矩阵
    return []

def animate(i):
    contour_set = ax.contourf(np.arange(64), np.arange(64), matrices[i], cmap=cmap)
    return contour_set.collections  # 返回需要更新的对象

ani = animation.FuncAnimation(fig, animate, frames=np.arange(0, nc, 1), init_func=init, interval=50, blit=True)

# 保存为gif
Writer = animation.writers['pillow']
writer = Writer(fps=10, metadata=dict(artist='Me'), bitrate=1800)
ani.save('my_contour_animation.gif', writer=writer)

plt.show()
    
'''
plt.show()
# 然后你可以访问存储在列表中的任意矩阵
for i, matrix in enumerate(matrices):
    plt.subplot(2, 5, i+1)  # 假设你想在一个大图中显示所有矩阵
    plt.imshow(matrix, cmap='gray')
    plt.xticks([]), plt.yticks([])  # 移除坐标轴刻度
plt.show()

fig, ax = plt.subplots()
initial_image = ax.imshow(matrices[0], cmap='viridis')  # 初始化第一幅图像
def update_map(num):
    initial_image.set_data(matrices[num])  # 更新已有图像的数据
    # 可选：根据需求更新标题
    # ax.set_title(f't={num * dt}')  

ani = animation.FuncAnimation(
    fig=fig,
    func=update_map,
    frames=list(range(0, nc + 1, 5)),
    interval=100,
    repeat=True)

# 保存为gif
#ani.save('twostream_.gif', writer='pillow')

# 显示动画
plt.show()

'''
