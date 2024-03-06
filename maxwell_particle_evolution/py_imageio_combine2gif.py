import imageio.v2 as imageio
import os  

nparticle = 60000  
ntime = 5000  # 假设粒子位置数组的长度代表总的时间步数  
ndiag = 20  # 假设每100个时间步进行一次诊断  
  
# 创建一个子文件夹（如果不存在）  
subfolder = 'scatter_plots_nparticle_{:d}_ntime_{:d}_ndiag_{:d}'.format(nparticle,ntime,ndiag)  
#subfolder = 'scatter_plots_py_many'   
# 设置GIF动画的帧率（例如，每秒10帧）  
fps = 10  
  
# 获取子文件夹中所有PNG文件的路径  
image_paths = [os.path.join(os.getcwd(), subfolder, f) for f in os.listdir(os.path.join(os.getcwd(), subfolder)) if f.endswith('.png')]  
  
# 排序文件路径，确保它们按数字顺序排列  
image_paths.sort(key=lambda x: int(os.path.splitext(x)[0].split('_')[-1]))  
  
# 使用imageio保存为GIF  
with imageio.get_writer(os.path.join(os.getcwd(), 'scatter_animation_nparticle_{:d}_ntime_{:d}_ndiag_{:d}.gif'.format(nparticle,ntime,ndiag) ), mode='I', fps=fps) as writer:  
    for path in image_paths:  
        image = imageio.imread(path)  
        writer.append_data(image)  
  
print("GIF animation created successfully!")