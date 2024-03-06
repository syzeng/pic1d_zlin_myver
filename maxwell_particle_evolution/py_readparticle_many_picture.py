import numpy as np  
import matplotlib.pyplot as plt  
import os  
  
# 指定文件名  
filex = 'particle_x.txt'  
filev = 'particle_v.txt'  
  
# 读取数据文件  
with open(filex, 'r') as fx, open(filev, 'r') as fv:  
    arx = np.loadtxt(fx)  
    arv = np.loadtxt(fv)  
  
nparticle = 60000  
ntime = 5000  # 假设粒子位置数组的长度代表总的时间步数  
ndiag = 20  # 假设每100个时间步进行一次诊断  
mystep=1 #每隔mystep个点选一个点
  
# 创建一个子文件夹（如果不存在）  
subfolder = 'scatter_plots_nparticle_{:d}_ntime_{:d}_ndiag_{:d}'.format(nparticle,ntime,ndiag)  
  
# 拼接当前工作目录和子文件夹路径  
full_path_dir = os.path.join(os.getcwd(), subfolder)  
  
# 检查子文件夹是否存在  
if not os.path.exists(full_path_dir):  
    # 如果不存在，则创建子文件夹  
    os.makedirs(full_path_dir)
  
# 循环遍历每个诊断并生成散点图  
for idiag in range(1, ntime // ndiag + 1):  
    start_idx = (idiag - 1) * nparticle  
    end_idx = idiag * nparticle  
    #mystep=10 #每隔mystep个点选一个点
  
    plt.figure(figsize=(12, 8), dpi=100)  
    marker_size = 0.1  
    plt.scatter(arx[start_idx:end_idx:mystep], arv[start_idx:end_idx:mystep], s=marker_size, color='blue')  
  
    # 添加标题和坐标轴标签  
    title_str = 'time step={:d}'.format((idiag - 1) * ndiag + 1)  # 修正标题以显示正确的时间步  
    plt.title(title_str)  
    plt.xlabel('x')  
    plt.ylabel('v')  
  
    # 显示网格  
    plt.grid(True)  
  
    # 定义文件名，包含idiag的值  
    filename_base = 'scatter_plot_idiag_'  
    file_extension = '.png'  
    filename = f"{filename_base}{idiag}{file_extension}"  
  
    # 拼接完整的文件路径  
    full_path = os.path.join(os.getcwd(), subfolder, filename)  
  
    # 保存散点图到指定路径  
    plt.savefig(full_path)  
    plt.close()
    # 显示图形（如果需要）  
    # plt.show()  # 如果你想要在每次迭代中都显示图形，可以取消注释这一行  
  
# 注意：如果你不想在每次迭代中都显示图形，可以移除或注释掉plt.show()。  
# 当循环结束后，所有的图都会被保存到指定的子文件夹中。