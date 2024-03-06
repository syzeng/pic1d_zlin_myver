import numpy as np  
import matplotlib.pyplot as plt  
import os  
from matplotlib.font_manager import FontProperties 
# from pylab import mpl
# # 设置显示中文字体
# mpl.rcParams["font.sans-serif"] = ["SimHei"]
  
# 指定文件名  
filex = 'particle_x.txt'  
filev = 'particle_v.txt'  
  
# 读取数据文件  
with open(filex, 'r') as fx, open(filev, 'r') as fv:  
    arx = np.loadtxt(fx)  
    arv = np.loadtxt(fv)  
  
nparticle = 60000  
ntime = 5000  
ndiag = 100  
nspecies = 2  
  
# 计算第idiag次诊断的索引范围  
idiag = 3 
start_idx = (idiag - 1) * nparticle  
end_idx = idiag * nparticle  
mystep=10 #每隔mystep个点选一个点

plt.figure(figsize=(12, 8), dpi=100)  
# 创建散点云图  
marker_size = 0.1  

plt.scatter(arx[start_idx:end_idx:mystep], arv[start_idx:end_idx:mystep], s=marker_size, color='blue')  
  
# 添加标题和坐标轴标签  
title_str = 'time step={:d} '.format((idiag-1)*ndiag)  
#font = FontProperties(fname=r"c:\windows\ fonts\STFANGSO.TTF", size=14)  # 指定字体文件路径和大小
plt.title(title_str)  
plt.xlabel('x')  
plt.ylabel('v')  
  
# 显示网格  
plt.grid(True)  
  
# 定义子文件夹名称和文件名  
subfolder = 'scatter_plots_one_test_nparticle_{:d}'.format(nparticle)  
filename_base = 'scatter_plot_idiag_{:d}'.format(idiag)  
file_extension = '.png'  
  
# 拼接文件名，包含idiag的值  
filename = f"{filename_base}{idiag}{file_extension}"  
  
# 拼接完整的文件路径  
full_path = os.path.join(os.getcwd(), subfolder, filename)  
  
# 如果子文件夹不存在，则创建它  
if not os.path.exists(os.path.join(os.getcwd(), subfolder)):  
    os.makedirs(os.path.join(os.getcwd(), subfolder))  
  
# 保存散点图到指定路径  
plt.savefig(full_path)  
  
# 显示图形  
plt.show()