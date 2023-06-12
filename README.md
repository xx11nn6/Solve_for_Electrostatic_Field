# Solve_for_Electrostatic_Field

#### 介绍
本程序是电子光学系统的计算机辅助设计程序(CAD)，用于求解多电极简单排布的空间静电场电位分布(一、二类像管)，可根据要求输出空间电位、轴上电位和扫描点电位及坐标，并绘制等位线。
#### 软件架构11
软件架构说明
               定义网格，迭代链表和一些全局变量
main.cpp---┬---Grid_Array.h---------------┬------struct Grid_Array                  创建结构电场网格
           |                              |------struct Iteration_Process           创建迭代链表
           |                              └------global variables                   声明全局变量
           |
           |   初始化网格函数
           |---Grid_Initialize.h-----------------void grid_initialize_mode_1/2      用于网格的初始化
           |         
           |   迭代计算电位相关函数            
           |---Iteration.h----------------┬------residual                           计算迭代残差
           |                              |------convergence_criteria               计算迭代精度
           |                              |------SOR                                一次SOR迭代
           |                              └------select_accelerator_factor          选择加速因子
           | 
           |   文件操作相关函数
           |---File_Operation.h-----------┬------readdata1/2                        读取数据
           |                              |------writedata                          写数据
           |                              |------scan_all_1/2                       扫描电位及坐标输出
           |                              └------rescan_1/2                         扫描一次
           |
           |   等位线绘制相关函数
           └---Paint_Equipotential_Line.h-┬------paint_all_1/2                      绘制等位线
                                          |------scan_1/2                           扫描等电位点
                                          └------paint_1                            绘制一条等位线

                                          
#### 安装教程

1.  打开vs就能跑

#### 使用说明

1.  像管
1.  将待读取数据存入本程序根目录下
2.  

#### 参与贡献

1.  Fork 本仓库
2.  新建 Feat_xxx 分支
3.  提交代码
4.  新建 Pull Request


#### 特技

1.  使用 Readme\_XXX.md 来支持不同的语言，例如 Readme\_en.md, Readme\_zh.md
2.  Gitee 官方博客 [blog.gitee.com](https://blog.gitee.com)
3.  你可以 [https://gitee.com/explore](https://gitee.com/explore) 这个地址来了解 Gitee 上的优秀开源项目
4.  [GVP](https://gitee.com/gvp) 全称是 Gitee 最有价值开源项目，是综合评定出的优秀开源项目
5.  Gitee 官方提供的使用手册 [https://gitee.com/help](https://gitee.com/help)
6.  Gitee 封面人物是一档用来展示 Gitee 会员风采的栏目 [https://gitee.com/gitee-stars/](https://gitee.com/gitee-stars/)
