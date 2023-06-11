#include <iostream>
#include <malloc.h>
#include<stdlib.h>
#include <numeric>
#include <iomanip>
#include <graphics.h>		// 引用 EasyX 图形库
#include <conio.h>
#include "Grid_Array.h"			//定义网格，全局变量M、N和释放网格内存函数
#include "Grid_Initialize.h"	//用于初始化两类像管
#include "Iteration.h"			//包括SOR迭代，加速因子计算，残差计算和迭代精度判断函数
#include "File_Operation.h"		//包括文件的读写
using namespace std;

//定义全局变量
void potential_line(Grid_Array** grid, int n, int* _N, int M1, int M2);
void paint(int k, Grid_Array** grid);
int scan(double V, Grid_Array** grid);
void paint_all(int n, int M1, int M2, int M, int N, Grid_Array** grid, double z0, double r0, double* dz, double delta);
double* V_r, * V_z; //存扫描点坐标
double* zi, * ri;

int mode;			 //工作模式1为一类，2为二类
int M;			  	 //网格数MN
int N;
double delta;		 // 电极厚度
int n;				 // 电极数（包含荧光屏，不包含阴极）
double* dz;           // 每个电极的轴向间距，两类像管均用到此变量
int* _N;              // 相邻电极间划分步长(通用)
double* V;            // 各个电极电压
double* VI;           // 含鞍点电极电位(通用)
int M1, M2;           // 竖直方向格数划分的要求，第一类像管使用此变量
double r1, r2;         // r2为电极插入电场的深度，r1为电极底端到轴的距离(第一类像管使用)
int* _M;              // 电极之间径向所需要划分的网格数，第二类像管使用此变量
double* dr;           // 电极内孔半径(第二类)
double epsilon;       // 迭代精度(通用)
double omega;		 // 加速因子ω
int NST;              // 输出打印空间电位时网格点间隔数(通用)
int INS;              // 轴上电位做等距插值时步长数(通用)
int* V1;              // 要求扫描等电位线的电位间隔或者电位值(通用)
int i;                // 循环用参数
FILE* handle;
errno_t err;
// 一个计数器，判断扫描等位线的条数




int main()
{
#if 0
	//初始化
	mode = 2;
	double V_s;		//荧光屏电压
	//以下为第一类像管用到的变量
	//////////////////注意！！M*N为网数，而_N与_M是格数///////////////////////////////
	if (mode == 1)  //第一类像管的参数
	{
		//赋初值
		V[0] = 24; V[1] = 40; V[2] = 62; V[3] = 74; V[4] = 85; V[5] = 96;
		dz[0] = 5.2; dz[1] = 8.6; dz[2] = 8.6; dz[3] = 8.6; dz[4] = 8.6; dz[5] = 8.6; dz[6] = 5.2;
		_N[0] = 3; _N[1] = 5; _N[2] = 5; _N[3] = 5; _N[4] = 5; _N[5] = 5; _N[6] = 2;
		M1 = 11; M2 = 7;
		//计算总网数M与N
		M = M1 + M2 + 1;  //格数是M1+M2,网数还要加一
		N = _N[0] + _N[1] + _N[2] + _N[3] + _N[4] + _N[5] + _N[6] + n; //有_N+n-1个格数，网数需要加一
	}

	else if (mode == 2)  //第二类像管参数
	{
		//赋初值
		V[0] = 24; V[1] = 40; V[2] = 62; V[3] = 74; V[4] = 85; V[5] = 96;
		dr[0] = 5.2; dr[1] = 8.6; dr[2] = 8.6; dr[3] = 8.6; dr[4] = 8.6; dr[5] = 8.6;
		dz[0] = 5.2; dz[1] = 8.6; dz[2] = 8.6; dz[3] = 8.6; dz[4] = 8.6; dz[5] = 8.6; dz[6] = 8.6; dz[7] = 5.2;
		_N[0] = 3; _N[1] = 5; _N[2] = 5; _N[3] = 5; _N[4] = 5; _N[5] = 5; _N[6] = 5; _N[7] = 3;
		_M[0] = 3; _M[1] = 5; _M[2] = 5; _M[3] = 5; _M[4] = 5; _M[5] = 5;
		//计算总网数
		M = _M[0] + _M[1] + _M[2] + _M[3] + _M[4] + _M[5] + n;
		N = _N[0] + _N[1] + _N[2] + _N[3] + _N[4] + _N[5] + _N[6] +_N[7]+ 1;
	}


	//new用于动态分配内存
	//例如：int* p = new int[10];
	//这句话表示分配10个int类型的空间，返回首地址存储在p中
	//因此可以用来创建一个可变大小的数组
	//要想创建二维数组，则必须使用双重指针
	Grid_Array** grid = new Grid_Array * [M];  //第一层指针指向行
	for (int i = 0; i < M; i++)
	{
		grid[i] = new Grid_Array[N];  //第二层指向列
	}

	grid_initialize_mode_2(grid, V_s, n, dz, dr, _N, _M, V, delta);
	
	//grid_initialize_mode_1(grid, V_s, n, dz, _N, V, r1, r2, M1, M2, delta);  //初始化第一类像管
	///*
	omega = select_accelerator_factor(grid);  //加速因子omega选取
	do
	{
		SOR(grid, omega);  //进行迭代直至符合精度条件
	} while (convergence_criteria(grid) > epsilon);

	cout << "accelerator factor omega:" << endl;  //输出加速因子
	cout << omega << endl;
	cout << "iteration times:" << endl;  //输出迭代次数
	cout << grid[M-2][N-2].k << endl;
	//*/
	//输出网格电位
	cout << "grid:" << endl;
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout << " " << setfill(' ') << setw(4) << setprecision(6) << grid[i][j].voltage_before << " ";  //setfill等函数是iomanip库中的函数，用于控制输出格式
			//cout << "(" << setprecision(2) << grid[i][j].h3 << "," << setprecision(2) << grid[i][j].h4 << ") ";
		}
		cout << endl;
	}
	//cout << sum(grid)/(M*N) << endl;
	//绘图
	//paint_all(n, M1, M2, M, N, grid, z0, r0, dz, delta);
	system("pause");
	free_grid(grid);//运行结束，释放内存
#endif
	
	//读文件捏
	if ((err = fopen_s(&handle, "test2.txt", "rt")) != 0)
	{
		printf("找不到文件！错误码：%d\n", err);
		return -1;
	}


	fopen_s(&handle, "test2.txt", "rt");
	if (fscanf_s(handle, "%d\n", &mode) == 1 && mode == 2)
	{
		readdata2("test2.txt");
		fclose(handle);
	}
	else
	{
		readdata1("test2.txt");
		fclose(handle);
	}

	double V_s = V[n - 1];  //定义荧光屏电压
	//计算M和N
	if (mode == 1)
	{
		M = M1 + M2 + 1;
		N = 0;
		for (int j = 0; j < n; j++)
		{
			N = N + _N[j] + 1;
		}
	}
	else
	{
		M = 1;
		N = 1;
		for (int i = 0; i < n-1; i++)
		{
			M = M + _M[i] + 1;
		}
		for (int i = 0; i < n + 1; i++)
		{
			N = N + _N[i];
		}
	}

	//创建网格grid
	//new用于动态分配内存
	//例如：int* p = new int[10];
	//这句话表示分配10个int类型的空间，返回首地址存储在p中
	//因此可以用来创建一个可变大小的数组
	//要想创建二维数组，则必须使用双重指针
	Grid_Array** grid = new Grid_Array * [M];  //第一层指针指向行
	for (int i = 0; i < M; i++)
	{
		grid[i] = new Grid_Array[N];  //第二层指向列
	}

	//初始化像管
	if (mode == 1)
	{
		grid_initialize_mode_1(grid, V_s, n, dz, _N, V, r1, r2, M1, M2, delta);  //初始化第一类像管
	}
	else
	{
		grid_initialize_mode_2(grid, V_s, n, dz, dr, _N, _M, V, delta);
	}

	omega = select_accelerator_factor(grid);  //加速因子omega选取
	do
	{
		SOR(grid, omega);  //进行迭代直至符合精度条件
	} while (convergence_criteria(grid) > epsilon);

	cout << "accelerator factor omega:" << endl;  //输出加速因子
	cout << omega << endl;
	cout << "iteration times:" << endl;  //输出迭代次数
	cout << grid[M - 2][N - 2].k << endl;
	//*/
	//输出网格电位
	cout << "grid:" << endl;
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout << " " << setw(4)<< setfill(' ')  << setprecision(4) << grid[i][j].voltage_before << " ";  //setfill等函数是iomanip库中的函数，用于控制输出格式
			//cout << "(" << setprecision(2) << grid[i][j].h3 << "," << setprecision(2) << grid[i][j].h4 << ") ";
		}
		cout << endl;
	}

	free(dz);
	free(_N);
	free(V);
	free(VI);
	free(dr);
	free(_M);
	free(V1);

	return 0;
}




void paint_all(int n, int M1, int M2, int M, int N, Grid_Array** grid, double z0, double r0, double* dz, double delta)
//画电位线
{
	int i, j, k, num;
	double z_temp1, z_temp2, V_temp1, V_temp2, t;
	V_z = (double*)malloc(n * (M1 + M2) * sizeof(double));
	V_r = (double*)malloc(n * (M1 + M2) * sizeof(double));
	zi = (double*)malloc(N * sizeof(double));
	ri = (double*)malloc(M * sizeof(double));
	zi[0] = 0;
	for (i = 1; i < N; i++)
	{
		zi[i] = zi[i - 1] + grid[0][i - 1].h2;  //累加求和，计算每一点离原点的轴向距离，等效于matlab中cumsum
	}
	ri[0] = 0;
	for (i = 1; i < M; i++)
	{
		ri[i] = grid[i][0].r;
	}

	for (i = 0; i < n * (M1 + M2); i++)
	{
		V_z[i] = 0;
		V_r[i] = 0;
	}
	initgraph(10 * z0, 10 * r0);  //初始化绘图窗口，宽为z0*10，高为r0*10
	setorigin(0, 10 * r0);  //设置原点为（0，r0*10）
	setaspectratio(1, -1);  //设置当前缩放因子,翻转y轴，使y朝上为正

	t = dz[0];  //t表示轴向，用于定位电极位置
	for (int i = 1; i < n; i = i + 1)  //用蓝色线绘制电极左半部分
	{
		setcolor(RGB(0, 0, 255));
		line(t * 10, r0 * 10, t * 10, grid[M2][0].r * 10);  //line(x1,y1,x2,y2)函数，输入起点，终点
		t = t + dz[i] + delta;
	}

	//用蓝色线绘制电极右半部分
	t = dz[0] + delta;
	for (int i = 1; i < n; i = i + 1)
	{
		setcolor(RGB(0, 0, 255));
		line(t * 10, r0 * 10, t * 10, grid[M2][0].r * 10);
		t = t + dz[i] + delta;
	}

	//用蓝色线连接电极
	t = dz[0];
	for (int i = 1; i < n; i = i + 1)
	{
		setcolor(RGB(0, 0, 255));
		line(t * 10, grid[M2][0].r * 10, (t + delta) * 10, grid[M2][0].r * 10);
		t = t + dz[i] + delta;
	}
	k = scan(64, grid);
	paint(k, grid);
	system("pause");
	closegraph();
}

int scan(double V, Grid_Array** grid)
//扫描，输入待扫描电位，以及网格grid
{
	int i, j, k = 0;
	double r, z;
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N - 1; j++)
		{
			if ((grid[i][j].voltage <= V && grid[i][j + 1].voltage > V) || (grid[i][j].voltage >= V && grid[i][j + 1].voltage < V))  //判断待扫描电位的位置区间
			{
				z = zi[j] + grid[0][j].h2 * (V - grid[i][j].voltage) / (grid[i][j + 1].voltage - grid[i][j].voltage);  //插值法计算待扫描电位位置
				V_z[k] = 200 + 10 * z;
				V_r[k] = 310 + 10 * (ri[M - 1] - ri[i]);
				k++;
			}
		}
	}
	for (j = 0; j < N; j++)
	{
		for (i = 0; i < M - 1; i++)
		{
			if ((grid[i][j].voltage <= V && grid[i + 1][j].voltage > V) || (grid[i][j].voltage >= V && grid[i + 1][j].voltage < V))
			{
				r = ri[i] + grid[i][0].h4 * (V - grid[i][j].voltage) / (grid[i + 1][j].voltage - grid[i][j].voltage);
				V_z[k] = 200 + 10 * zi[j];
				V_r[k] = 310 + 10 * (ri[M - 1] - r);
				k++;
			}
		}
	}
	return k;
}

void paint(int k, Grid_Array** grid) //画线//
{
	int i, j;
	double temp, x, y, xy;
	x = grid[0][0].h2;
	for (i = 0; i < N - 1; i++)  //判断轴向间距最小值
	{
		if (x < grid[0][i].h2)
			x = grid[0][i].h2;
	}
	y = grid[0][0].h4;
	if (y < grid[M - 2][0].h4)
		y = grid[M - 2][0].h4;
	xy = sqrt(pow(x, 2) + pow(y, 2));
	for (i = 1; i < k; i++)
	{
		for (j = 1; j <= k - i; j++)  //冒泡排序
		{
			if (V_z[j - 1] >= V_z[j])
			{
				temp = V_z[j - 1];
				V_z[j - 1] = V_z[j];
				V_z[j] = temp;
				temp = V_r[j - 1];
				V_r[j - 1] = V_r[j];
				V_r[j] = temp;
			}
		}
	}
	for (i = 0; i < k; i++)
	{
		for (j = i + 1; j < k; j++)
		{
			if (sqrt(pow(V_r[j] - V_r[i], 2) + pow(V_z[j] - V_z[i], 2)) < 10.2 * xy &&
				fabs(V_r[j] - V_r[i]) < 10.2 * y && fabs(V_z[j] - V_z[i]) < 10.2 * x)
			{
				line(int(V_z[i]), int(V_r[i]), int(V_z[j]), int(V_r[j]));
				break;
			}
		}
	}
	for (i = 1; i < k; i++)
	{
		for (j = 1; j <= k - i; j++)
		{
			if (V_z[j - 1] <= V_z[j])
			{
				temp = V_z[j - 1];
				V_z[j - 1] = V_z[j];
				V_z[j] = temp;
				temp = V_r[j - 1];
				V_r[j - 1] = V_r[j];
				V_r[j] = temp;
			}
		}
	}
	for (i = 0; i < k; i++)
	{
		for (j = i + 1; j < k; j++)
		{
			if (sqrt(pow(V_r[j] - V_r[i], 2) + pow(V_z[j] - V_z[i], 2)) < 12 * xy &&
				fabs(V_r[j] - V_r[i]) < 12 * y && fabs(V_z[j] - V_z[i]) < 12 * x)
			{
				line(int(V_z[i]), int(V_r[i]), int(V_z[j]), int(V_r[j]));
				break;
			}
		}
	}
}