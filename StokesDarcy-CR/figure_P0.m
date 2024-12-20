% 创建一个简单的模型
model = createpde;

% 定义几何（如单位圆）
geometryFromEdges(model, @circleg);

% 生成三角形网格
mesh = generateMesh(model, 'Hmax', 0.2);

% 获取网格节点和单元信息
nodes = mesh.Nodes;        % 网格节点坐标
elements = mesh.Elements;  % 网格单元定义

% 假设有 P0 元数值解，每个单元对应一个解值
% 例如，给每个单元随机赋值
P0_solution = rand(size(elements, 2), 1);  % 每个单元一个常数解

% 提取三角形顶点坐标
x = nodes(1, :);
y = nodes(2, :);
tri = elements';

% trisurf 方法，利用三角形单元及每个单元的常数值绘制
trisurf(tri, x, y, P0_solution, 'EdgeColor', 'none'); 
colorbar;
title('P0 Solution Visualization');
xlabel('x');
ylabel('y');
view(2); % 设置为俯视视角