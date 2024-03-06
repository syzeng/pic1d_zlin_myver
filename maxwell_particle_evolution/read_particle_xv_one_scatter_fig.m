clear
% 指定文件名  
filex = 'particle_x.txt';  
filev='particle_v.txt';
  
% 打开文件用于读取  
fileIDx = fopen(filex, 'r');  
fileIDv = fopen(filev, 'r');  
  
% 使用textscan读取每一行的数据，%f表示读取浮点数  
arx = textscan(fileIDx, '%f', 'Delimiter', '\n');  
arv=textscan(fileIDv, '%f', 'Delimiter', '\n');  
  
% 关闭文件  
fclose(fileIDx);  
fclose(fileIDv);   
% textscan返回一个cell数组，我们需要将其转换为普通数组  
arx = arx{1};  
arv = arv{1}; 

nparticle=60000;
ntime=1000;
ndiag=100;
nspecies=2;
%诊断次数：ntime/ndiag;每次诊断nparticle个粒子的信息

% 创建散点云图  
% 设置散点的颜色为蓝色，标记为实心圆，并设定点的大小（例如：10个点的直径单位）
markerSize = 1;
idiag=3;    %第idiag次诊断
% 计算第idiag次诊断的索引范围  
startIdx = (idiag-1)*nparticle + 1;  
endIdx = idiag*nparticle;  
scatter(arx(startIdx:endIdx), arv(startIdx:endIdx), markerSize, 'b', 'filled');

% 添加标题和坐标轴标签  
% 创建格式化的标题字符串  
titleStr = sprintf('粒子散点云图（第 %d 次诊断）', idiag);  
% 设置标题  
title(titleStr);  
xlabel('x');  
ylabel('v');  
  
% 显示网格
grid on;

% 获取当前图形句柄  
xviFig = gcf;  
  
% 设置图形大小（以像素为单位）  
set(xviFig, 'Position', [100, 100, 800, 600]);  
  
% 保存散点图  
%title_save_fig=sprintf('粒子散点云图（第 %d 次诊断）.png', idiag);  
%saveas(xviFig, title_save_fig);

% ... (绘制散点图的代码)  
  
% 定义子文件夹名称和文件名  
subfolder = 'scatter_plots';  
filename_base = 'scatter_plot_idiag_';  
file_extension = '.png';  
  
% 拼接文件名，包含idiag的值  
filename = [filename_base, num2str(idiag), file_extension];  
  
% 拼接完整的文件路径  
full_path = fullfile(pwd, subfolder, filename);  
  
% 如果子文件夹不存在，则创建它  
if ~isfolder(fullfile(pwd, subfolder))  
    mkdir(fullfile(pwd, subfolder));  
end  
  
% 保存散点图到指定路径  
saveas(gcf, full_path);
