# Solve_for_Electrostatic_Field

#### 介绍
本程序是电子光学系统的计算机辅助设计程序(CAD)，用于求解多电极简单排布的空间静电场电位分布(一、二类像管)，可根据要求输出空间电位、轴上电位和扫描点电位及坐标，并绘制等位线。
#### 软件架构
软件架构说明
1.  main.cpp
    - Grid_Array.h
        - struct Grid_Array          创建结构电场网格
        - struct Iteration_Process   创建迭代链表
        - global_variables           声明全局变量
    - Grid_Initialize.h
        - void grid_initialize_mode_1  用于网格的初始化
        - void grid_initialize_mode_2
    - Iteration.h
        - double residual                  计算迭代残差
        - double convergence_criteria      计算迭代精度
        - void SOR                       一次SOR迭代
        - double select_accelerator_factor 选择加速因子
    - File_Operation.h
        - void readdata1  读取数据
        - void readdata2
        - void writedata    写数据
        - void scan_all_1   扫描电位及坐标输出
        - void scan_all_2
        - int rescan_1      扫描一次
        - int rescan_2    
    - Paint_Equipotential_Line.h
        - void paint_all_1  绘制等位线
        - void paint_all_2     
        - void scan_1       扫描等电位点
        - void scan_2     
        - void paint_1      绘制一条等位线
        - void paint_2
                                    
#### 安装教程

1.  打开vs就能跑

#### 使用说明

1.  一类像管如图 ![Image text](https://gitee.com/xx11nn6/solve_for_-electrostatic_-field/raw/master/%E4%B8%80%E7%B1%BB%E5%83%8F%E7%AE%A1.png)
2.  二类像管如图![Image text](https://gitee.com/xx11nn6/solve_for_-electrostatic_-field/raw/master/%E4%BA%8C%E7%B1%BB%E5%83%8F%E7%AE%A1.png)
3.  输入要求如图![Image text](https://gitee.com/xx11nn6/solve_for_-electrostatic_-field/raw/master/%E8%BE%93%E5%85%A5%E5%8F%98%E9%87%8F%E5%8F%8A%E6%A0%BC%E5%BC%8F%E8%A6%81%E6%B1%82.png)
4.  将待读取数据存入本程序根目录下
5.  在main.cpp中更改待读取文件和路径、输出文件名及路径
6.  运行程序
7.  空间电位等结果保存在指定路径下，等位线以窗口形式输出

#### 参与贡献

1.  魏千怀-负责 main.cpp , Grid_Array.h , Grid_Initialize.h , Iteration.h 和 File_Operation.h 的 writedata 函数
2.  李煜翔-负责 Paint_Equipotential_Line.h 和 File_Operation.h 的 scan_all, rescan 函数
3.  葛军韬-负责 File_Operation.h 函数
4.  xx11nn6-fork 负责可爱捏 ^^ 联系方式xx11nn6@126.com


#### 结语

1.  完结撒花！！
2.  当听完算法描述的时候，我们都小看了这次编程，都以为会很简单。然而实际上手起来才发现事情不对，实际工程量远大于想象。说起来，大学三年，我们编程只学过c语言，甚至没学过数据结构，却要完成工程量如此大的项目。在这短短一个月内，还有各种考试和大作业，实属是榨干了我们每一个人。大家断断续续努力着，所幸总算是在答辩之前把项目写完了。回顾整个工程，在开头的时候我还信誓旦旦想尝试一下面向对象，但是由于没有专门学过，写出来的代码总是奇奇怪怪，于是放弃了。现在回头看来，工程貌似庞大，但是里面完全没有复杂的代码，也许是因为我们太菜了写不出高级的代码，几乎都是if和for解决的。仔细想想我们代码似乎确实过于冗余，大家实现的方式都奇奇怪怪，但最终代码确实能跑，简直就是shit mountain。
3.  本次项目确实狠狠提升了我们的编程能力，尤其是对指针的运用。最重要的是，我们确实体会到了代码分块和问题分解的意义，每一个功能的实现似乎都很困难，但是将其分解成的多个简单的问题，似乎就迎刃而解了。总的来说，这次编程确实难，但带给我们收获也是巨大的。完成答辩的时候，大伙儿长叹了一口气，终于完了，我甚至有种空虚的感觉，感觉每天都在改代码，总算写完了却有些不知所措了。这次项目在最终的分数里只占40%，但我们依然用尽全力去完成了。在最后非常感谢我的队友李煜翔和葛军韬，大家有困难一起讨论，一起努力解决，这是我在大学组过最靠谱的小组了，谢谢你们！
