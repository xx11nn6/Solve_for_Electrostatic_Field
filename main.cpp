#include <iostream>
#include <malloc.h>
#include <stdlib.h>
#include <numeric>
#include <iomanip>
#include <graphics.h>		// 引用 EasyX 图形库
#include <conio.h>
#include "Grid_Array.h"			//定义网格，全局变量M、N和释放网格内存函数
#include "Grid_Initialize.h"	//用于初始化两类像管
#include "Iteration.h"			//包括SOR迭代，加速因子计算，残差计算和迭代精度判断函数
#include "File_Operation.h"		//包括文件的读写
#include "Paint_Equipotential_Line.h"
using namespace std;

//定义全局变量
//读入参数
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
double r1, r2;        // r2为电极插入电场的深度，r1为电极底端到轴的距离(第一类像管使用)
int* _M;              // 电极之间径向所需要划分的网格数，第二类像管使用此变量
double* dr;           // 电极内孔半径(第二类)
double epsilon;       // 迭代精度(通用)
int NST;              // 输出打印空间电位时网格点间隔数(通用)
int INS;              // 轴上电位做等距插值时步长数(通用)
int* V1;              // 要求扫描等电位线的电位间隔或者电位值(通用)

//迭代变量
double omega_1;		  // 加速因子ω，无鞍点
double omega_2;
int iteration_times_1;  // 迭代次数，无鞍点
int iteration_times_2;  // 迭代次数，有鞍点
bool round_p = false;   // 判断一轮迭代是否结束
int round_n;		    // 当前轮数
int times_n;			// 当前次数
Iteration_Process* head1;  //头结点，无鞍点
Iteration_Process* head2;  //头结点，有鞍点

//绘制等位线变量
double z0, r0;		   // 网格总长/宽
double* V_r, * V_z;    // 存扫描点坐标
double* zi, * ri;	   // 存轴向/径向距离
int count1;		       // 等位线计数器
int tmp;			   // 输出等位线要求

//读写文件变量
FILE* handle_read;	//读文件的句柄
FILE* handle_write; //写文件的句柄
errno_t err;


int main()
{
	//读文件捏
	count1 = 0;
	if ((err = fopen_s(&handle_read, "test6.txt", "rt")) != 0)
	{
		printf("找不到文件！错误码：%d\n", err);
		return -1;
	}


	fopen_s(&handle_read, "test6.txt", "rt");
	if (fscanf_s(handle_read, "%d\n", &mode) == 1 && mode == 2)
	{
		readdata2(handle_read);
		fclose(handle_read);
	}
	else
	{
		readdata1(handle_read);
		fclose(handle_read);
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

	//无鞍点网格
	Grid_Array** grid1 = new Grid_Array * [M];  //第一层指针指向行
	for (int i = 0; i < M; i++)
	{
		grid1[i] = new Grid_Array[N];  //第二层指向列
	}


	//有鞍点网格
	Grid_Array** grid2 = new Grid_Array * [M];  //第一层指针指向行
	for (int i = 0; i < M; i++)
	{
		grid2[i] = new Grid_Array[N];  //第二层指向列
	}

	cout << "\ncount:" << count1 << endl;

	//初始化像管
	if (mode == 1)
	{
		grid_initialize_mode_1(grid1, V_s, n, dz, _N, V, r1, r2, M1, M2, delta);  //初始化第一类像管(无鞍点)
		grid_initialize_mode_1(grid2, V_s, n, dz, _N, VI, r1, r2, M1, M2, delta);

		z0 = 0;
		for (int i = 0; i < n; i++) {
			z0 = dz[i] + z0;
		}
		z0 = z0 + (n - 1) * delta;
		r0 = r1 + r2;
	}
	else
	{
		grid_initialize_mode_2(grid1, V_s, n, dz, dr, _N, _M, V, delta);
		grid_initialize_mode_2(grid2, V_s, n, dz, dr, _N, _M, VI, delta);
	}

	//创建无鞍点网格迭代过程链表
	
	iteration_times_1 = 0;
	Iteration_Process* ite1 = new Iteration_Process();
	ite1->next = nullptr;
	head1 = ite1;
	head1->omega_r = omega_1;
	head1->avg_res = residual(grid1);
	head1->max_res = convergence_criteria(grid1);
	head1->round = 1;
	head1->times = 1;
	round_n = 1;
	times_n = 0;

	omega_1 = select_accelerator_factor(grid1, ite1,&iteration_times_1);  //加速因子omega选取
	do
	{
		SOR(grid1, omega_1, ite1, &iteration_times_1);  //进行迭代直至符合精度条件
	} 
	while (ite1->max_res >= epsilon);


	//有鞍点网格迭代过程链表
	iteration_times_2 = 0;
	Iteration_Process* ite2 = new Iteration_Process();
	ite2->next = nullptr;
	head2 = ite2;
	head2->omega_r = omega_2;
	head2->avg_res = residual(grid2);
	head2->max_res = convergence_criteria(grid2);
	head2->round = 1;
	head2->times = 1;
	round_n = 1;
	times_n = 0;

	omega_2 = select_accelerator_factor(grid2, ite2,&iteration_times_2);  //加速因子omega选取
	do
	{
		SOR(grid2, omega_2, ite2,&iteration_times_2);  //进行迭代直至符合精度条件
	} while (convergence_criteria(grid2) > epsilon);	


	//绘制等位线
	if (mode == 1)
	{
		paint_all_1(n, M1, M2, M, N, grid1, z0, r0, dz, delta, tmp, V1, count1);
		paint_all_1(n, M1, M2, M, N, grid2, z0, r0, dz, delta, tmp, V1, count1);
	}
	else
	{
		paint_all_2(n, M, N, grid1, dz, delta, _M, dr, _N, tmp, V1, count1);
		paint_all_2(n, M, N, grid2, dz, delta, _M, dr, _N, tmp, V1, count1);
	}
	
	

	cout << "accelerator factor omega:" << endl;  //输出加速因子
	cout << omega_1 << endl;
	cout << "iteration times:" << endl;  //输出迭代次数
	cout << iteration_times_1 << endl;
	//*/
	//输出网格电位
	cout << "grid:" << endl;
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			;
			//cout << " " << setw(4)<< setfill(' ')  << setprecision(4) << grid[i][j].r << " ";  //setfill等函数是iomanip库中的函数，用于控制输出格式
			//cout << "(" << setprecision(2) << grid[i][j].h3 << "," << setprecision(2) << grid[i][j].h4 << ") ";
		}
		//cout << endl;
	}


	//写文件，输出
	if ((err = fopen_s(&handle_write, "输出.res", "w")) != 0)
	{
		printf("找不到文件！错误码：%d\n", err);
		return -1;
	}

	writedata(handle_write, grid1, grid2, head1, head2);
	fclose(handle_write);



	free(dz);
	free(_N);
	free(V);
	free(VI);
	free(dr);
	free(_M);
	free(V1);
	free_grid(grid1);//运行结束，释放内存
	free_grid(grid2);//运行结束，释放内存
	return 0;
}
