#include <iostream>
#include <malloc.h>
#include <numeric>
#include <iomanip>
using namespace std;

struct Grid_Array  //定义电场网格结构
{
	int k;  //存储迭代次数k
	double voltage;  //存储各网格点迭代k次的电压
	double voltage_before;  //存储上一次迭代的电压
	bool is_margin;  //判断网格点是否是边界，1表示是边界
	bool on_axis;  //判断是否在轴上；
	double h1, h2, h3, h4;  //存储网格点到周围的距离
	double r;  //存储网格点的径向距离
};


//定义全局变量
int M, N;//网格数M*N
int mode = 1;  //定义工作模式，1为第一类像管，2为第二类像管
//声明函数
void grid_initialize_mode_1(Grid_Array** grid, double cathode_voltage, double screen_voltage, int n, double* dz, int* _N, double* V, double r1, double r2, int M1, int M2, double delta);
double residual(Grid_Array** grid);
double SOR(Grid_Array** grid, double omega);
double select_accelerator_factor(Grid_Array** grid);
double convergence_criteria(Grid_Array** grid);
void free(Grid_Array** grid);

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
	N = _N[0] + _N[1] + _N[2] + _N[3] + _N[4] + _N[5] + _N[6] + n; //有_N+n-1个格数，网数需要加一

	double omega;  //定义加速因子ω
	double epsilon = 0.005;  //迭代精度为0.005



	//c++中，new用于动态分配内存
	//例如：int* p = new int[10];
	//这句话表示分配10个int类型的空间，返回首地址存储在p中
	//因此可以用来创建一个可变大小的数组
	//要想创建二维数组，则必须使用双重指针
	Grid_Array** grid = new Grid_Array * [M];  //第一层指针指向行
	for (int i = 0; i < M; i++)
	{
		grid[i] = new Grid_Array[N];  //第二层指向列
	}
	grid_initialize_mode_1(grid, V_c, V_s, n, dz, _N, V, r1, r2, M1, M2, delta);
	omega = select_accelerator_factor(grid);
	do
	{
		SOR(grid, omega);
	} while (convergence_criteria(grid) > epsilon);
	cout << "iteration times:" << endl;
	cout << grid[5][5].k << endl;


	//c++中用cout来输出，使用cout输出矩阵
	cout << "grid:" << endl;
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout << " " << setfill(' ') << setw(4) << setprecision(6) << grid[i][j].voltage << " ";  //setfill等函数是iomanip库中的函数，用于控制输出格式
			//cout << "(" << setprecision(2) << grid[i][j].h3 << "," << setprecision(2) << grid[i][j].h4 << ") ";
		}
		cout << endl;
	}
	//cout << sum(grid)/(M*N) << endl;
	free(grid);//运行结束，释放内存

	return 0;
}




//用于初始化第一类像管的电场网格
void grid_initialize_mode_1(Grid_Array** grid, double cathode_voltage, double screen_voltage, int n, double* dz, int* _N, double* V, double r1, double r2, int M1, int M2, double delta)
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
			grid[i][j].k = 0;

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
				while (k < n)
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

					else if (k == (n - 1)) //最后一组
					{
						if (j == sum) //在荧光屏上
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


			//将初始迭代第1次电压设置与第0次相同
			for (int i = 0; i < M; i++)
			{
				for (int j = 0; j < N; j++)
				{
					grid[i][j].voltage_before = grid[i][j].voltage;
				}
			}
		}
	}

	//设置荧光屏电位和边界情况
	for (int i = 0; i < M; i++)
	{
		grid[i][N - 1].voltage = screen_voltage;
		grid[i][N - 1].is_margin = true;
	}
}

//用于计算残差的均值
double residual(Grid_Array** grid)
{
	double sum = 0;
	double diff;
	double res;
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			diff = grid[i][j].voltage - grid[i][j].voltage_before;
			sum += diff;
		}
	}
	res = sum / (M * N);
	return res;
}

//进行一次SOR迭代
double SOR(Grid_Array** grid, double omega)
{
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{

			if (grid[i][j].is_margin == false)//判断是否是边界，不是边界则参与迭代
			{
				double c1, c2, c3, c4, c0;
				//为了公式好看，先用变量暂存
				double phi = grid[i][j].voltage;  //迭代前的电压φ[k]
				double phi_after;  //迭代后的电位φ[k+1]
				double phi_avg;  //平均电压
				double r = grid[i][j].r;
				double h1 = grid[i][j].h1;
				double h2 = grid[i][j].h2;
				double h3 = grid[i][j].h3;
				double h4 = grid[i][j].h4;
				if (grid[i][j].on_axis == false)  //非轴上点的迭代
				{
					c1 = 2 / (h1 * (h1 + h2));
					c2 = 2 / (h2 * (h1 + h2));
					c3 = (2 * r - h4) / (r * h3 * (h3 + h4));
					c4 = (2 * r + h3) / (r * h4 * (h3 + h4));
					c0 = c1 + c2 + c3 + c4;
					//先进行一次普通（赛德尔-利伯曼）迭代
					phi_after = (c1 * grid[i][j - 1].voltage + c2 * grid[i][j + 1].voltage \
						+ c3 * grid[i + 1][j].voltage + c4 * grid[i - 1][j].voltage) / c0;
					//用SOR计算:φ[k+1]=φ[k]+ω(φ[k+1]_bar-φ[k])
					phi_avg = (phi + phi_after) / 2;
					phi_after = phi + omega * (phi_avg - phi);
					grid[i][j].voltage_before = phi;
					grid[i][j].voltage = phi_after;
					grid[i][j].k += 1;
				}
				else  //轴上点的迭代，其格式相同，只不过系数不同
				{
					c1 = 2 / (h1 * (h1 + h2));
					c2 = 2 / (h2 * (h1 + h2));
					c3 = 0;
					c4 = 4 / (h4 * h4);
					c0 = 2 * ((1 / (h1 * h2)) + (2 / (h4 * h4)));
					//迭代没有c3项
					phi_after = (c1 * grid[i][j - 1].voltage \
						+ c2 * grid[i][j + 1].voltage \
						+ c4 * grid[i - 1][j].voltage) / c0;
					phi_avg = (phi + phi_after) / 2;
					phi_after = phi + omega * (phi_avg - phi);
					grid[i][j].voltage_before = phi;
					grid[i][j].voltage = phi_after;
					grid[i][j].k += 1;
				}
			}
		}
	}
	return 0;
}


//最佳加速因子ω的选择
double select_accelerator_factor(Grid_Array** grid)
{
	double E_bar;
	double E_bar_after;
	double lambda = 0;  //先使λ=0
	double omega_l, mu_l, omega_m;  //ω_λ，μ_λ，ω_m
	double omega_m_bar;  //用于两轮迭代算出的加速因子ω_m平均值

	double omega = 1;  //取omega=1迭代一次
	SOR(grid, omega);
	omega = 1.375;
	do
	{
		lambda = 0;
		for (int i = 0; i < 12; i++)
		{
			E_bar = residual(grid);  //计算第一次残差
			SOR(grid, omega);  //迭代一次
			E_bar_after = residual(grid);  //再次计算残差
			if (i > 8)//最后三次迭代
			{
				lambda += (E_bar_after / E_bar);  //取最后三次λ相加
			}
		}
		lambda /= 3;  //计算最后三次迭代λ均值
		mu_l = (lambda + omega - 1) / (sqrt(lambda) * omega);
		omega_l = 2 / (1 + sqrt(1 - (pow(mu_l, 2))));
		omega_m = 1.25 * omega_l - 0.5;
		omega_m_bar = (omega + omega_m) / 2;
		omega = omega_m;
		//cout << omega << endl;
	} while (abs((omega_m_bar - omega_m) / (2 - omega_m)) >= 0.05);
	return omega;
}


//迭代精度判断
double convergence_criteria(Grid_Array** grid)
{
	double max_error, temp;
	max_error = 0;
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (grid[i][j].is_margin == false)//对非边界点判断
			{
				temp = abs(grid[i][j].voltage - grid[i][j].voltage_before);//用网格点与上一次迭代相比
				if (temp > max_error)
				{
					max_error = temp;  //取所有误差中的最大值
				}
			}
		}
	}
	return max_error;
}

//c++中使用delete来释放内存
void free(Grid_Array** grid)
{
	for (int i = 0; i < M; i++) {
		delete[] grid[i];
	}
	delete[] grid;
}
