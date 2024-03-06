import numpy as np  
import matplotlib.pyplot as plt  
  
# 读取文件内容到numpy数组  
array = np.loadtxt('phi_mode.out', dtype=float)  
  
# 绘制折线图  
plt.plot(array)  
  
# 设置图表标题和坐标轴标签  
plt.title('Array Plot')  
plt.xlabel('Index')  
plt.ylabel('Value')  
  
# 显示图表  
plt.show()

# 应用快速傅里叶变换 (FFT)
fft_result = np.fft.fft(array)

# 获取频率轴上的值（通常是对称的，所以只需要展示一半）
n = len(array)
freqs = np.fft.fftfreq(n)  # 注意：这里返回的是归一化的频率
half_n = n // 2

# 绘制原始信号和频谱
plt.figure(figsize=(10, 4))

# 时间域信号
plt.subplot(121)
plt.plot(array)
plt.title('Original Signal')

# 频域信号
plt.subplot(122)
plt.plot(freqs[:half_n], np.abs(fft_result[:half_n]) / n)  # 取绝对值并除以样本数量以获取幅度谱
plt.xlabel('Frequency')
plt.ylabel('Amplitude')
plt.title('FFT Spectrum')

plt.tight_layout()
plt.show()