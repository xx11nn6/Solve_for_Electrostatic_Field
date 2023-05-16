#include <iostream>
#include <malloc.h>
#include <numeric>
#include <iomanip>
using namespace std;

struct Grid_Array  //定义电场网格结构
{
	double voltage;  //存储各网格点电压
	bool is_margin;  //判断网格点是否是边界，1表示是边界
	bool on_axis;  //判断是否在轴上；
	double h1, h2, h3, h4;  //存储网格点到周围的距离
	double r;  //存储网格点的径向距离
};


//定义全局变量
int M, N;//网格数M*N
int mode = 1;  //定义工作模式，1为第一类像管，2为第二类像管
//声明函数
void grid_initialize_mode_1(Grid_Array** grid, double cathode_voltage, double screen_voltage,int n, double* dz, int* _N, double* V, double r1, double r2, int M1, int M2, double delta);


int main()
{
	//初始化
	//////////////////注意！！M*N为网数，而_N与M1、M2是格数///////////////////////////////
	double V_c = 0, V_s = 100;		//定义光阴极和荧光屏电压分别为0和100V
	int n = 7;						//定义电极数（包含荧光屏，不包含阴极）
	double delta = 0.5;				//定义电极宽度δ
	double z0 = 53.4, r0 = 32;		//定义径向和轴向的宽度
	double r2 = 12, r1 = r0 - r2;   //r2为电极插入电场的深度，r1为电极底端到轴的距离
	double dz[7];				    //定义每个电极的间距
	int _N[7];						//定义每个电极之间所取的格数（水平方向），加下划线为了与全局变量N区别
	double V[6];					//定义各个电极电压
	int M1, M2;					    //竖直方向格数划分的要求

	V[0] = 24; V[1] = 40; V[2] = 62; V[3] = 74; V[4] = 85; V[5] = 96;
	dz[0] = 5.2; dz[1] = 8.6; dz[2] = 8.6; dz[3] = 8.6; dz[4] = 8.6; dz[5] = 8.6; dz[6] = 5.2;
	_N[0] = 3; _N[1] = 5; _N[2] = 5; _N[3] = 5; _N[4] = 5; _N[5] = 5; _N[6] = 2;
	M1 = 11; M2 = 7;
	//计算总网数M与N
	M = M1 + M2 + 1;  //格数是M1+M2,网数还要加一
	N = _N[0] + _N[1] + _N[2] + _N[3] + _N[4] + _N[5] +_N[6]+ n; //有_N+n-1个格数，网数需要加一
	double epsilon = 0.0005;  //迭代精度为0.005



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
	grid_initialize_mode_1(grid,V_c,V_s,n,dz,_N,V,r1,r2,M1,M2,delta);


	//c++中用cout来输出，使用cout输出矩阵
	cout << "grid:" << endl;
	for (int i = 0; i < M; i++) 
	{
		for (int j = 0; j < N; j++) 
		{
			cout << " " << setfill(' ') << setw(4) << setprecision(3) << grid[i][j].voltage << " ";  //setfill等函数是iomanip库中的函数，用于控制输出格式
			//cout << "(" << setprecision(2) << grid[i][j].h3 << "," << setprecision(2) << grid[i][j].h4 << ") ";
		}
		cout << endl;
	}

	free(grid);//运行结束，释放内存

	return 0;
}




//用于初始化第一类像管的电场网格
void grid_initialize_mode_1(Grid_Array** grid,double cathode_voltage,double screen_voltage,int n,double* dz,int* _N,double* V,double r1,double r2,int M1,int M2,double delta)
					//输入电场网格，网格宽度、高度、阴极电压、荧光屏电压、电极个数n、网格间距（传入数组dz）、水平方向网格划分要求（传入数组_N）、电极电压（数组V）
					//电极底到轴距r1、电极深度r2，垂直方向网格划分要求M1,M2、电极宽度delta
{
	//遍历网格进行初始化
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			grid[i][j].voltage = 0;  //初始电位全部赋零
			grid[i][j].is_margin = false;  //初始全部设置为非边界
			grid[i][j].on_axis = false;  //初始设置为非轴上点

			//设置垂直间距h3,h4
			//设置第一行
			if (i == 0)
			{
				grid[i][j].h4 = 0;
				grid[i][j].h3 = r2 / M2;
				grid[i][j].r = r1 + r2;
			}
			//设置1~M2行
			else if (i < M2)
			{
				grid[i][j].h4 = r2 / M2;
				grid[i][j].h3 = grid[i][j].h4;
				grid[i][j].r = r1 + r2 - ((r2 * i) / M2);
			}
			//设置M2行
			else if (i == M2)
			{
				grid[i][j].h4 = r2 / M2;
				grid[i][j].h3 = r1 / M1;
				grid[i][j].r = r1;
			}
			//设置M2~M1行
			else if (i < (M2 + M1))
			{
				grid[i][j].h4 = r1 / M1;
				grid[i][j].h3 = grid[i][j].h4;
				grid[i][j].r = r1 - ((r1 * (i - M2)) / M1);
			}
			//设置轴上点
			else
			{
				grid[i][j].h3 = 0;
				grid[i][j].h4 = r1 / M1;
				grid[i][j].r = 0;
				grid[i][j].on_axis = true;
				grid[i][j].is_margin = true;
			}

			//第一列（光阴极）
			if (j == 0)
			{
				grid[i][j].voltage = cathode_voltage;
				grid[i][j].h1 = 0;
				grid[i][j].h2 = dz[0] / (_N[0]);
				grid[i][j].is_margin = true;
			}
			else	
			{
				//后面几列需要用循环判断
				int sum = _N[0];
				int k = 0;
				while(k < n)
				{
					if (j < sum)
					{
						if (i == 0)  //如果i=0（第一行）直接做封闭边界处理
						{
							grid[i][j].is_margin = true;
							if (k == 0)  //若在阴极和第一个电极间，线性插值计算边界电压
							{
								grid[i][j].voltage = (V[0] - cathode_voltage) * j / _N[0];
							}
							else if (k != (n - 1))  //若不在最后一个电极到荧光屏区间
							{
								grid[i][j].voltage = V[k - 1] + (V[k] - V[k - 1]) * (_N[k] + j - sum) / _N[k];
							}
							else  //最后一个电极到荧光屏区间
							{
								grid[i][j].voltage = V[k - 1] + (screen_voltage - V[k - 1]) * (_N[k] + j - sum) / _N[k];
							}
						}
						grid[i][j].h1 = dz[k] / _N[k];
						grid[i][j].h2 = dz[k] / _N[k];
						break;
					}
					else if (j == sum && (k != (n - 1))) //在电极左侧且不在荧光屏
					{
						if (i <= M2)  //若高度在电极深度范围内
						{
							grid[i][j].is_margin = true;
							grid[i][j].voltage = V[k];
						}
						grid[i][j].h1 = dz[k] / _N[k];
						grid[i][j].h2 = delta;
						break;
					}
					else if (j == (sum + 1) && (k != (n - 1)))  //在电极右侧且不在荧光屏
					{
						if (i <= M2)  //若高度在电极深度范围内
						{
							grid[i][j].is_margin = true;
							grid[i][j].voltage = V[k];
						}
						grid[i][j].h1 = delta;
						grid[i][j].h2 = dz[k + 1] / _N[k + 1];
						break;
					}

					else if (k == (n-1) ) //最后一组
					{
						if (j == sum ) //在荧光屏上
						{
							grid[i][j].h1 = dz[k] / _N[k];
							grid[i][j].h2 = 0;
							grid[i][j].is_margin = true;
							break;
						}
					}
					else  //若越过第一个电极
					{
						k += 1;
						sum = sum + _N[k] + 1;
						continue;
					}
				}
			}
			
		}
	}
	
	//设置荧光屏电位和边界情况
	for (int i = 0; i < M; i++)
	{
		grid[i][N-1].voltage = screen_voltage;
		grid[i][N-1].is_margin = true;
	}
}



//c++中使用delete来释放内存
void free(Grid_Array** grid)
{
	for (int i = 0; i < M; i++) {
		delete[] grid[i];
	}
	delete[] grid;
}
