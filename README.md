## BoPaint 光栅图形学作业

#### 依赖

- OpenCV
- nlohmann/json: https://github.com/nlohmann/json

#### 基本使用

- 功能：
  - [x] Bresenham画线
  - [x] 中点法画圆弧
  - [x] 扫描线填充多边形
  - [x] 直线加权采样抗锯齿
  - [x] 圆弧加权采样抗锯齿
  - [x] 多边形边界加权采样抗锯齿
- 编译程序：``make``
- 运行程序：``bopaint [input json] [output filename]``，其中``input json``为图像的描述文件，``output filename``是输出的路径

#### 本地环境

- macOS 10.14.3
- OpenCV 3.4.3
- gcc 4.2.1 (Apple LLVM 10.0.0, clang-1000.11.45.5)

#### 输入文件格式

- 输入文件格式需为json，下面提供一个样例

  ```json
  {
      "width": 1024,
      "height": 1024,
      "waaa": true,
      "elaa": 1,
      "lines": [
          [640, 400, 960, 300]
      ],
      "arcs": [
          [410, 512, 256, -120, 180]
      ],
      "polygon": [
          {"fill": true, "points": [215, 440, 260, 430, 245, 475]}
      ]
  }
  ```

- 其中``width``、``height``参数分别为图片宽度和高度，均为整数类型；``waaa``为是否开启加权区域采样抗锯齿，为Boolean类型；``elaa``为分辨率放大倍数，为整数类型

- ``lines``为需要画出的直线，为数组类型，数组的每个元素为包含4个整数的子数组，4个整数分别表示起点的XY坐标和终点的XY坐标

- ``arcs``为需要画出的圆弧，为数组类型，数组的每个元素为包含5个整数的子数组，5个整数分别表示圆心的XY坐标、半径、起始角度和终止角度，其中角度的范围是$[-180^{\circ},180^{\circ}]$

- ``polygon``为需要画出的多边形，为数组类型，数组的每个元素为一个字典。字典的``fill``属性代表是否填充多边形，``points``属性是多边形的点，依次是每个点的XY坐标，这里暂时不支持自相交的多边形

- 如下为一组对应的输入输出

  ![](/Users/lyricz/Desktop/CG/BoPaint/example.jpeg)

#### 算法描述

- 上述画线、画圆弧、画多边形和抗锯齿的算法实现均在``graph.cpp/h``中
- 画线采用了Bresenham算法，有几点需要注意的事情
  - 需要注意斜率$k\ge 1$和$k<1$的两种情况，为了保证最好的效果需要对坐标进行翻转
  - 也需要注意$k<0$的情况，增量$e$等其他变量均为负数，也需要特殊判断
  - 对应的抗锯齿采用了加权区域采样，计算Bresenham需要画的像素点上下两个点的灰度并进行对应填充，其中线宽只取了1个像素，在高分辨率的时候效果并不是很好
- 画圆采用了中点画圆法
  - 利用了圆的对称性，只计算了八分之一的圆弧之后进行对称变换
  - 其中抗锯齿的灰度计算结果也同理对称，减少了大量的计算
  - 抗锯齿和直线唯一不同的是区域采样的判别不同，这里需要判断区域是否在线宽为1的两个圆之间
- 画多边形采用了教材上的扫描线法
  - 不同的是和扫描线相交的边通过``std:: set``来维护，提高了效率
  - 对于边界的抗锯齿，这里采用了再调用Bresenham算法对边界进行重新填充，对应的边界抗锯齿效果也会由补充上去