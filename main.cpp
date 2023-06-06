#include <iostream>
#include <malloc.h>
#include <numeric>
#include <iomanip>
#include <graphics.h>		// 引用 EasyX 图形库
#include <conio.h>

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
int M, N;//网数M*N
int mode = 1;  //定义工作模式，1为第一类像管，2为第二类像管
//声明函数
void grid_initialize_mode_1(Grid_Array** grid, double screen_voltage, int n, double* dz, int* _N, double* V, double r1, double r2, int M1, int M2, double delta);
void grid_initialize_mode_2(Grid_Array** grid, double screen_voltage, int n, double* dz, double* dr, int* _N, int* _M, double* V, double delta);
double residual(Grid_Array** grid);
double SOR(Grid_Array** grid, double omega);
double select_accelerator_factor(Grid_Array** grid);
double convergence_criteria(Grid_Array** grid);
void free(Grid_Array** grid);
void potential_line(Grid_Array** grid, int n, int* _N, int M1, int M2);
void paint(int k, Grid_Array** grid);
int scan(double V, Grid_Array** grid);
void paint_all(int n, int M1, int M2, int M, int N, Grid_Array** grid, double z0, double r0, double* dz, double delta);
double* V_r, * V_z; //存扫描点坐标
double* z1, * r1;

int main()
{
	//初始化
	mode = 2;
	///以下为两类像管通用参数
	double V_s = 100;				//定义荧光屏电压为100V
	int n = 7;						//定义电极数（包含荧光屏，不包含阴极）
	double delta = 0.5;				//定义电极宽度δ
	double V[6];					//定义各个电极电压
	int _N[7];						//电极之间横向所需要划分的网格数，加下划线为了与全局变量N区别
	int _M[6];						//电极之间径向所需要划分的网格数，第二类像管使用此变量
	int M1, M2;					    //竖直方向格数划分的要求，第一类像管使用此变量
	double dz[7];				    //定义每个电极的轴向间距，两类像管均用到此变量
	double dr[6];					//定义电极径向间距，第二类像管用此变量
	//以下为第一类像管用到的变量
	double z0 = 56.4, r0 = 32;		//定义径向和轴向的宽度
	double r2 = 12, r1 = r0 - r2;   //r2为电极插入电场的深度，r1为电极底端到轴的距离
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
		dz[0] = 5.2; dz[1] = 8.6; dz[2] = 8.6; dz[3] = 8.6; dz[4] = 8.6; dz[5] = 8.6; dz[6] = 5.2;
		_N[0] = 3; _N[1] = 5; _N[2] = 5; _N[3] = 5; _N[4] = 5; _N[5] = 5; _N[6] = 3;
		_M[0] = 3; _M[1] = 5; _M[2] = 5; _M[3] = 5; _M[4] = 5; _M[5] = 5;
		//计算总网数
		M = _M[0] + _M[1] + _M[2] + _M[3] + _M[4] + _M[5] + n;
		N = _N[0] + _N[1] + _N[2] + _N[3] + _N[4] + _N[5] + _N[6] + 1;
	}

	double omega;  //定义加速因子ω
	double epsilon = 0.0005;  //迭代精度为0.0005



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
	/*grid_initialize_mode_1(grid, V_s, n, dz, _N, V, r1, r2, M1, M2, delta);  //初始化第一类像管
	omega = select_accelerator_factor(grid);  //加速因子omega选取
	do
	{
		SOR(grid, omega);  //进行迭代直至符合精度条件
	} while (convergence_criteria(grid) > epsilon);

	cout << "accelerator factor omega:" << endl;  //输出加速因子
	cout << omega << endl;
	cout << "iteration times:" << endl;  //输出迭代次数
	cout << grid[2][2].k << endl;
	*/
	//输出网格电位
	cout << "grid:" << endl;
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout << " " << setfill(' ') << setw(4) << setprecision(6) << grid[i][j].r << " ";  //setfill等函数是iomanip库中的函数，用于控制输出格式
			//cout << "(" << setprecision(2) << grid[i][j].h3 << "," << setprecision(2) << grid[i][j].h4 << ") ";
		}
		cout << endl;
	}
	//cout << sum(grid)/(M*N) << endl;
	//绘图
	//paint_all(n, M1, M2, M, N, grid, z0, r0, dz, delta);
	system("pause");
	free(grid);//运行结束，释放内存
	return 0;
}




//用于初始化第一类像管的电场网格
void grid_initialize_mode_1(Grid_Array** grid, double screen_voltage, int n, double* dz, int* _N, double* V, double r1, double r2, int M1, int M2, double delta)
//输入电场网格，荧光屏电压、电极个数n、网格间距（传入数组dz）、水平方向网格划分要求（传入数组_N）、电极电压（数组V）
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
				grid[i][j].voltage = 0;
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
								grid[i][j].voltage = V[0] * j / _N[0];
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
		grid[i][N - 1].voltage_before = screen_voltage;
	}
}



//用于初始化第二类像管的电场网格
void grid_initialize_mode_2(Grid_Array** grid, double screen_voltage, int n, double* dz, double* dr, int* _N, int* _M, double* V, double delta)
//输入电场网格，荧光屏电压、电极个数n、网格轴向间距（传入数组dz）、径向间距dr、水平方向网格划分要求（传入数组_N）、垂直方向网格划分要求_M、电极电压（数组V）、电极宽度delta
{
	//分块 必须有数组存储网格数划分要求的累加和，还有轴向距离
	int* Ms = new int[n-1];
	int* Ns = new int[n];
	double* rs = new double[n - 1];

	Ms[0] = _M[n - 2] + 1;  //本来是_M[0]的，但是没想到吧哈哈！要求网格划分和坐标系是反着的......
	Ns[0] = _N[0];
	rs[0] = dr[0];
	for (int i = 1; i < n; i++)
	{
		if (i != (n - 1))
		{
			Ms[i] = _M[n - i - 2] + Ms[i - 1] + 1;
			Ns[i] = _N[i] + Ns[i - 1];
			rs[i] = dr[i] + rs[i - 1];
		}
		else
		{
			Ns[i] = _N[i] + Ns[i - 1];
		}
	}
	//先遍历网格，把所有点置零
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			grid[i][j].voltage = 0;  //初始电位全部赋零
			grid[i][j].is_margin = false;  //初始全部设置为非边界
			grid[i][j].on_axis = false;  //初始设置为非轴上点
			grid[i][j].k = 0;
			grid[i][j].r = 0;
		}
	}

	//再次遍历，进行初始化和边界设置
	for (int k = 0; k < n; k++)
		//分块处理，每一块是两个电极间
	{
		if (k == 0)//从上往下 第一个电极到第二个之间
		{
			for (int i = 0; i < Ms[k]; i++)  
			{
				for (int j = 0; j < N; j++)
				{
					if (i == 0 || i == 1)  //如果是第一行或是第二行
					{
						if ((j >= Ns[n - k - 3]) && (j <= Ns[n - k - 2]))  //如果在电极区间
						{
							grid[i][j].voltage = V[n - k - 2];
							grid[i][j].is_margin = true;
						}
						else if (j > Ns[n - k - 2])  //如果在电极到荧光屏间
						{
							if (i == 0)  //设置第一行为边界插值
							{
								grid[i][j].voltage = V[n - k - 2] + double(j - Ns[n - k - 2]) / (Ns[n - k - 1] - Ns[n - k - 2]) * (screen_voltage - V[n - k - 2]);  //行间线性插值计算电压
								grid[i][j].is_margin = true;
							}
						}
					}

					else if ((i > 1) && (i < Ms[k]))  //如果是第三行到下一个电极间
					{
						if (j == Ns[n - k - 3])  //垂直补充边界
						{
							float r;  //网格点的轴向距离r
							r = rs[n - k - 2] - (double(i - 1) / _M[n - k - 2] * dr[n - k - 2] - delta);
							grid[i][j].voltage = V[n - k - 2] + (V[n - k - 3] - V[n - k - 2]) * \
								log(r / (rs[n - k - 2] - delta)) / log(rs[n - k - 3] / rs[n - k - 2]);  //对数插值，给我写晕了
							grid[i][j].is_margin = true;
							grid[i][j].r = r;
						}
					}
					else
					{
						break;
					}
				}
			}
		}
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
					phi_after = (1 - omega) * phi + omega * phi_after;
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
					phi_after = (1 - omega) * phi + omega * phi_after;
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

//使用delete来释放内存
void free(Grid_Array** grid)
{
	for (int i = 0; i < M; i++) {
		delete[] grid[i];
	}
	delete[] grid;
}

void paint_all(int n, int M1, int M2, int M, int N, Grid_Array** grid, double z0, double r0, double* dz, double delta)
//画电位线
{
	int i, j, k, num;
	double z_temp1, z_temp2, V_temp1, V_temp2, t;
	V_z = (double*)malloc(n * (M1 + M2) * sizeof(double));
	V_r = (double*)malloc(n * (M1 + M2) * sizeof(double));
	z1 = (double*)malloc(N * sizeof(double));
	r1 = (double*)malloc(M * sizeof(double));
	z1[0] = 0;
	for (i = 1; i < N; i++)
	{
		z1[i] = z1[i - 1] + grid[0][i - 1].h2;  //累加求和，计算每一点离原点的轴向距离，等效于matlab中cumsum
	}
	r1[0] = 0;
	for (i = 1; i < M; i++)
	{
		r1[i] = grid[i][0].r;
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
				z = z1[j] + grid[0][j].h2 * (V - grid[i][j].voltage) / (grid[i][j + 1].voltage - grid[i][j].voltage);  //插值法计算待扫描电位位置
				V_z[k] = 200 + 10 * z;
				V_r[k] = 310 + 10 * (r1[M - 1] - r1[i]);
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
				r = r1[i] + grid[i][0].h4 * (V - grid[i][j].voltage) / (grid[i + 1][j].voltage - grid[i][j].voltage);
				V_z[k] = 200 + 10 * z1[j];
				V_r[k] = 310 + 10 * (r1[M - 1] - r);
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