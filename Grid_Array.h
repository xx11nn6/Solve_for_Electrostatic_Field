#ifndef GRID_ARRAY_H  //定义网格，释放函数和一些全局变量
#define GRID_ARRAY_H

//声明全局变量M
extern int mode;			 //工作模式1为一类，2为二类
extern int M;			  	 //网格数MN
extern int N;
extern double delta;		 // 电极厚度
extern int n;				 // 电极数（包含荧光屏，不包含阴极）
extern double* dz;           // 每个电极的轴向间距，两类像管均用到此变量
extern int* _N;              // 相邻电极间划分步长(通用)
extern double* V;            // 各个电极电压
extern double* VI;           // 含鞍点电极电位(通用)
extern int M1, M2;           // 竖直方向格数划分的要求，第一类像管使用此变量
extern double r1,r2;         // r2为电极插入电场的深度，r1为电极底端到轴的距离(第一类像管使用)
extern int* _M;              // 电极之间径向所需要划分的网格数，第二类像管使用此变量
extern double* dr;           // 电极内孔半径(第二类)
extern double epsilon;       // 迭代精度(通用)
extern double omega;		 // 加速因子ω
extern int NST;              // 输出打印空间电位时网格点间隔数(通用)
extern int INS;              // 轴上电位做等距插值时步长数(通用)
extern int* V1;              // 要求扫描等电位线的电位间隔或者电位值(通用)

//文件声明
extern FILE* handle;
extern errno_t err;

struct Grid_Array  //定义电场网格结构
{
	int k;  //存储迭代次数k
	double voltage;  //存储各网格点迭代k次的电压
	double voltage_before;  //存储上一次迭代的电压
	bool is_margin;  //判断网格点是否是边界，1表示是边界
	bool is_inside;  //判断网格点是否在边界内部
	bool on_axis;  //判断是否在轴上；
	double h1, h2, h3, h4;  //存储网格点到周围的距离
	double r;  //存储网格点的径向距离
};

//使用delete来释放内存
void free_grid(Grid_Array** grid)
{
	for (int i = 0; i < M; i++) {
		delete[] grid[i];
	}
	delete[] grid;
}

#endif
