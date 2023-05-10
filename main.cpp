#include <iostream>
#include <malloc.h>
using namespace std;

struct Grid_Array  //定义电场网格结构
{
	float voltage;  //存储各网格点电压
	bool is_margin;  //判断是否是边界，1表示是边界
};

int M, N;//网格数M*N

//声明函数
void grid_initialize(Grid_Array** grid, float cathode_voltage, float screen_voltage);
void node_set(Grid_Array** grid, int m, int n, float voltage, bool is_margin);


int main()
{
	//初始化
	float V_c, V_s; //定义光阴极和荧光屏电压
	M = 20; N = 20;
	V_c = 0; V_s = 100; //定义光阴极和荧光屏电压分别为0和100V
	int mode = 1;  //定义工作模式，1为第一类像管，2为第二类像管
	float z0, r0;  //定义径向和轴向的宽度
	

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
	grid_initialize(grid,V_c,V_s);
	node_set(grid,3,1,50,true);


	//c++中用cout来输出，使用cout输出矩阵
	cout << "grid:" << endl;
	for (int i = 0; i < M; i++) 
	{
		for (int j = 0; j < N; j++) 
		{
			cout << "(" << grid[i][j].voltage << ", " << grid[i][j].is_margin << ") ";
		}
		cout << endl;
	}

	free(grid);//运行结束，释放内存

	return 0;
}




//用于初始化电场网格
void grid_initialize(Grid_Array** grid,float cathode_voltage,float screen_voltage)
					//输入电场网格，网格宽度、高度、阴极电压和荧光屏电压
{
	//先将电位置零，全部设置为非边界
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			grid[i][j].voltage = 0;
			grid[i][j].is_margin = false;
		}
	}
	//设置光阴极电位和边界情况
	for (int i = 0; i < M; i++)
	{
		grid[i][0].voltage = cathode_voltage;
		grid[i][0].is_margin = true;
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
