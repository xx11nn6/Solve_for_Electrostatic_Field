#include <iostream>
#include <malloc.h>
#include <numeric>
using namespace std;

struct Grid_Array  //定义电场网格结构
{
	float voltage;  //存储各网格点电压
	bool is_margin;  //判断网格点是否是边界，1表示是边界
	bool on_axis;  //判断是否在轴上；
	float h1, h2, h3, h4;  //存储网格点到周围的距离
	float r;  //存储网格点的径向距离
};


//定义全局变量
int M, N;//网格数M*N
int mode = 1;  //定义工作模式，1为第一类像管，2为第二类像管
//声明函数
void grid_initialize(Grid_Array** grid, float cathode_voltage, float screen_voltage, float* dz, int* _N, float r1, float r2, int M1, int M2, float delta);
void node_set(Grid_Array** grid, int m, int n, float voltage, bool is_margin);


int main()
{
	//初始化
	float V_c = 0, V_s = 100;  //定义光阴极和荧光屏电压分别为0和100V
	int n = 7;  //定义电极数（包含荧光屏，不包含阴极）
	float delta = 0.5;//定义电极宽度δ
	float z0 = 53.4, r0 = 32;  //定义径向和轴向的宽度
	float r2 = 12, r1 = r0 - r2;  //r2为电极插入电场的深度，r1为电极底端到轴的距离
	float dz[6];  //定义每个电极的间距
	int _N[6];  //定义每个电极之间所取的网格数（水平方向），加下划线为了与全局变量N区别
	int M1, M2;  //竖直方向网格数划分的要求
	dz[0] = 5.2; dz[1] = 8.6; dz[2] = 8.6; dz[3] = 8.6; dz[4] = 8.6; dz[5] = 5.2;
	_N[0] = 3; _N[1] = 5; _N[2] = 5; _N[3] = 5; _N[4] = 5; _N[5] = 2;
	M1 = 11; M2 = 7;
	//计算总网格数M与N
	M = M1 + M2;
	N = _N[0] + _N[1] + _N[2] + _N[3] + _N[4] + _N[5] + n - 1;
	float epsilon = 0.0005;  //迭代精度为0.005



	//c++中，new用于动态分配内存
	//例如：int* p = new int[10];
	//这句话表示分配10个int类型的空间，返回首地址存储在p中
	//因此可以用来创建一个可变大小的数组
	//要想创建二维数组，则必须使用双重指针
	Grid_Array** grid = new Grid_Array* [M];  //第一层指针指向行
	for (int i = 0; i < M; i++) 
	{
		grid[i] = new Grid_Array[N];  //第二层指向列
	}
	grid_initialize(grid,V_c,V_s,dz,_N,r1,r2,M1,M2,delta);
	node_set(grid,3,1,50,true);


	//c++中用cout来输出，使用cout输出矩阵
	cout << "grid:" << endl;
	for (int i = 0; i < M; i++) 
	{
		for (int j = 0; j < N; j++) 
		{
			cout << "(" << grid[i][j].r << ", " << grid[i][j].h3 << "," << grid[i][j].h4 << ") ";
		}
		cout << endl;
	}


	free(grid);//运行结束，释放内存

	return 0;
}




//用于初始化电场网格
void grid_initialize(Grid_Array** grid,float cathode_voltage,float screen_voltage,float* dz,int* _N,float r1,float r2,int M1,int M2,float delta)
					//输入电场网格，网格宽度、高度、阴极电压、荧光屏电压、网格间距（传入数组dz）、水平方向网格划分要求（传入数组_N）
					//电极底到轴距r1、电极深度r2，垂直方向网格划分要求M1,M2、电极宽度delta
{
	//先将电位置零，全部设置为非边界、非轴上点，以及竖直方向间距
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			grid[i][j].voltage = 0;
			grid[i][j].is_margin = false;
			grid[i][j].on_axis = false;
			//设置第一行
			if (i == 0)
			{
				grid[i][j].h4 = 0;
				grid[i][j].h3 = r2 / (M2 - 1);
				grid[i][j].r = r1 + r2;
			}
			//设置1~M2-1行
			else if (i < M2 - 1)
			{
				grid[i][j].h4 = r2 / (M2 - 1);
				grid[i][j].h3 = grid[i][j].h4;
				grid[i][j].r = r1 + r2 - ((r2 * i) / (M2 - 1));
			}
			//设置M2行
			else if (i == M2 - 1)
			{
				grid[i][j].h4 = r2 / (M2 - 1);
				grid[i][j].h3 = r1 / (M1 - 1);
				grid[i][j].r = r1;
			}
			//设置M2~M1-1行
			else if (i < (M2 + M1 - 1))
			{
				grid[i][j].h4 = r1 / (M1 - 1);
				grid[i][j].h3 = grid[i][j].h4;
				grid[i][j].r = r1 - ((r1 * (i - M2)) / (M1 - 1));
			}
			//设置轴上点
			else
			{
				grid[i][j].h3 = 0;
				grid[i][j].h4 = r1 / (M1 - 1);
				grid[i][j].r = 0;
				grid[i][j].on_axis = true;
				grid[i][j].is_margin = true;
			}
		}
	}
	//设置光阴极
	for (int i = 0; i < M; i++)
	{
		grid[i][0].voltage = cathode_voltage;
		grid[i][0].is_margin = true;
		grid[i][0].h1 = 0;  //光阴极为左边界，h1=0
		grid[i][0].h2 = dz[0] / (_N[0] - 1);  //光阴极右边到第一个电极距离为dz[0],网格数为_N[0]，则可算出h2
		
		
	}
	//设置荧光屏电位和边界情况
	for (int i = 0; i < M; i++)
	{
		grid[i][N-1].voltage = screen_voltage;
		grid[i][N-1].is_margin = true;
	}
}

//用于给指定坐标的网格节点设定电压值和边界状态，输入grid数组，坐标（m,n），电压，是否是边界
void node_set(Grid_Array** grid,int m,int n,float voltage,bool is_margin)
{
	grid[m][n].voltage = voltage;
	grid[m][n].is_margin = is_margin;
}

//边界封闭函数，用于使静电场
void boundary_close()
{

}

//c++中使用delete来释放内存
void free(Grid_Array** grid)
{
	for (int i = 0; i < M; i++) {
		delete[] grid[i];
	}
	delete[] grid;
}
