#define GRID_INITIALIZE_H  //用于初始化两类像管
#include "Grid_Array.h"


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
			grid[i][j].is_inside = true;  //第一类像管除边界，其余全部参与迭代
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

	//设置边界内部
	for (int i = 1; i < M; i++)
	{
		for (int j = 1; j < (N - 1); j++)
		{
			grid[i][j].is_inside = true;
		}
	}
}



//用于初始化第二类像管的电场网格
void grid_initialize_mode_2(Grid_Array** grid, double screen_voltage, int n, double* dz, double* dr, int* _N, int* _M, double* V, double delta)
//输入电场网格，荧光屏电压、电极个数n、网格轴向间距（传入数组dz）、径向间距dr、水平方向网格划分要求（传入数组_N）、垂直方向网格划分要求_M、电极电压（数组V）、电极宽度delta
{
	//分块 必须有数组存储网格数划分要求的累加和，还有轴向距离
	int* Ms = new int[n - 1];
	int* Ns = new int[n + 1];
	double* rs = new double[n - 1];

	Ms[0] = _M[n - 2] + 1;  //本来是_M[0]的，但是没想到吧哈哈！要求网格划分和坐标系是反着的......
	Ns[0] = _N[0];
	rs[0] = dr[0] + delta;
	for (int i = 1; i < (n + 1); i++)
	{
		if (i < (n - 1))
		{
			Ms[i] = _M[n - i - 2] + Ms[i - 1] + 1;
			Ns[i] = _N[i] + Ns[i - 1];
			rs[i] = dr[i] + rs[i - 1] + delta;
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
			grid[i][j].voltage = 0.0;  //初始电位全部赋零
			grid[i][j].is_margin = false;  //初始全部设置为非边界
			grid[i][j].on_axis = false;  //初始设置为非轴上点
			grid[i][j].is_inside = false;//初始设置为非内部
			grid[i][j].r = 0;
			grid[i][j].h1 = 0;
			grid[i][j].h2 = 0;
			grid[i][j].h3 = 0;
			grid[i][j].h4 = 0;
		}
	}

	//再次遍历，进行初始化和边界设置
	for (int k = 0; k < (n - 1); k++)
		//分块处理，每一块是两个电极间
	{
		//从上往下 第一个电极到第二个之间
		if (k == 0)
		{
			for (int i = 0; i < Ms[k]; i++)
			{
				for (int j = 0; j < N; j++)
				{
					if (i == 0 || i == 1)  //如果是第一行或是第二行
					{
						if ((j >= Ns[n - k - 2]) && (j <= Ns[n - k - 1]))  //如果在电极区间
						{
							grid[i][j].voltage = V[n - k - 2];	//设置电极电压
							grid[i][j].is_margin = true;
						}
						else if (j > Ns[n - k - 1] && j < (N - 1))  //如果在电极到荧光屏间
						{
							if (i == 0)  //设置第一行为边界插值
							{
								grid[i][j].voltage = V[n - k - 2] + double(j - Ns[n - k - 1]) / (Ns[n - k] - Ns[n - k - 1]) * (screen_voltage - V[n - k - 2]);  //行间线性插值计算电压
								grid[i][j].is_margin = true;
							}
							else if (i == 1) //设置第二行
							{
								grid[i][j].is_inside = true;
								grid[i][j].r = rs[n - k - 2] - delta;
								grid[i][j].h1 = dz[n - k] / _N[n - k];
								grid[i][j].h2 = grid[i][j].h1;
								grid[i][j].h3 = dr[n - k - 2] / _M[n - k - 2];
								grid[i][j].h4 = delta;
							}
						}
					}

					else if ((i > 1) && (i < Ms[k]))  //如果是第三行到下一个电极间
					{
						if (j == Ns[n - k - 2])  //如果在两电极间，垂直补充边界
						{
							float r;  //网格点的轴向距离r
							r = rs[n - k - 2] - (double(i - 1) / _M[n - k - 2] * dr[n - k - 2]) - delta;
							grid[i][j].voltage = V[n - k - 2] + (V[n - k - 3] - V[n - k - 2]) * \
								log(r / (rs[n - k - 2] - delta)) / log(rs[n - k - 3] / rs[n - k - 2]);  //对数插值，给我写晕了
							grid[i][j].is_margin = true;
							grid[i][j].r = r;
						}
						else if (j > Ns[n - k - 2] && j < (N - 1))  //设置边界内部
						{
							grid[i][j].is_inside = true;
							grid[i][j].h3 = dr[n - k - 2] / _M[n - k - 2];
							grid[i][j].h4 = grid[i][j].h3;
							grid[i][j].r = rs[n - k - 2] - (double(i - 1) / _M[n - k - 2] * dr[n - k - 2]) - delta;
							for (int l = 0; l < (k + 2); l++)  //轴向划分k+2个区域
							{
								if (j > Ns[n - k + l - 2] && j < Ns[n - k + l - 1])
								{
									grid[i][j].h1 = dz[n - k + l - 1] / _N[n - k + l - 1];
									grid[i][j].h2 = grid[i][j].h1;
								}
								else if (j == Ns[n - k + l - 1] && j != (N - 1))  //若在两个划分区间中，且不在荧光屏上
								{
									grid[i][j].h1 = dz[n - k + l - 1] / _N[n - k + l - 1];
									grid[i][j].h2 = dz[n - k + l] / _N[n - k + l];
								}
							}
						}
					}
					else
					{
						break;
					}
				}
			}
		}

		//第二个到第n-1个电极
		else  if (k != (n - 2))
		{
			for (int i = Ms[k - 1]; i < Ms[k]; i++)
			{
				for (int j = Ns[n - k - 2]; j < (N - 1); j++)
				{
					if (i == Ms[k - 1])  //如果是电极行
					{
						if ((j >= Ns[n - k - 2]) && j <= (Ns[n - k - 1]))  //如果在电极区间
						{
							grid[i][j].voltage = V[n - k - 2];	//设置电极电压
							grid[i][j].is_margin = true;
						}
						else
						{
							grid[i][j].is_inside = true;
							grid[i][j].h3 = delta;
							grid[i][j].h4 = dr[n - k - 1] / _M[n - k - 1];
							grid[i][j].r = rs[n - k - 2];
						}
					}
					else if (i == (Ms[k - 1] + 1))//如果是电极下一行
					{
						if ((j >= Ns[n - k - 2]) && (j <= Ns[n - k - 1]))  //如果在电极区间
						{
							grid[i][j].voltage = V[n - k - 2];	//设置电极电压
							grid[i][j].is_margin = true;
						}
						else
						{
							grid[i][j].is_inside = true;
							grid[i][j].h3 = dr[n - k - 2] / _M[n - k - 2];
							grid[i][j].h4 = delta;
							grid[i][j].r = rs[n - k - 2] - delta;
						}
					}
					else  //如果不是电极两行
					{
						if (j == Ns[n - k - 2])//如果在两电极列，垂直补充边界
						{
							{
								grid[i][j].voltage = V[n - k - 2];	//设置电极电压
								float r;  //网格点的轴向距离r
								r = rs[n - k - 2] - (double(i - Ms[k - 1] - 1) / _M[n - k - 2] * dr[n - k - 2]) - delta;
								grid[i][j].voltage = V[n - k - 2] + (V[n - k - 3] - V[n - k - 2]) * \
									log(r / (rs[n - k - 2] - delta)) / log(rs[n - k - 3] / rs[n - k - 2]);  //对数插值
								grid[i][j].is_margin = true;
								grid[i][j].r = r;
							}
						}
						else
						{
							grid[i][j].is_inside = true;
							grid[i][j].h3 = dr[n - k - 2] / _M[n - k - 2];
							grid[i][j].h4 = grid[i][j].h3;
							grid[i][j].r = rs[n - k - 2] - (double(i - Ms[k - 1] - 1) / _M[n - k - 2] * dr[n - k - 2]) - delta;
						}

					}

					for (int l = 0; l < (k + 2); l++)  //轴向划分k+2个区域
					{
						if (j > Ns[n - k + l - 2] && j < Ns[n - k + l - 1])
						{
							grid[i][j].h1 = dz[n - k + l - 1] / _N[n - k + l - 1];
							grid[i][j].h2 = grid[i][j].h1;
						}
						else if (j == Ns[n - k + l - 1] && j != (N - 1))  //若在两个划分区间中，且不在荧光屏上
						{
							grid[i][j].h1 = dz[n - k + l - 1] / _N[n - k + l - 1];
							grid[i][j].h2 = dz[n - k + l] / _N[n - k + l];
						}
					}
				}
			}
		}

		//如果在最后一个电极区间
		else if (k == (n - 2))
		{
			for (int i = Ms[k - 1]; i <= Ms[k]; i++)
			{
				for (int j = 0; j < (N - 1); j++)
				{
					if (i == Ms[k - 1])  //如果是电极行
					{
						if ((j >= Ns[n - k - 2]) && j <= (Ns[n - k - 1]))  //如果在电极区间
						{
							grid[i][j].voltage = V[n - k - 2];	//设置电极电压
							grid[i][j].is_margin = true;
						}
						else if (j > Ns[n - k - 1])
						{
							grid[i][j].is_inside = true;
							grid[i][j].h3 = delta;
							grid[i][j].h4 = dr[n - k - 1] / _M[n - k - 1];
							grid[i][j].r = rs[n - k - 2];
						}
					}
					else if (i == (Ms[k - 1] + 1))//如果是电极下一行
					{
						if (j > 0 && j < Ns[n - k - 2])	//线性插值计算电压
						{
							grid[i][j].voltage = double(j) / _N[0] * V[0];
							grid[i][j].is_margin = true;
						}
						else if ((j >= Ns[n - k - 2]) && (j <= Ns[n - k - 1]))  //如果在电极区间
						{
							grid[i][j].voltage = V[n - k - 2];	//设置电极电压
							grid[i][j].is_margin = true;
						}
						else if (j != 0)
						{
							grid[i][j].is_inside = true;
							grid[i][j].h3 = dr[n - k - 2] / _M[n - k - 2];
							grid[i][j].h4 = delta;
							grid[i][j].r = rs[n - k - 2] - delta;
						}
					}
					else  //如果不是电极两行
					{
						if (j != 0)
						{
							grid[i][j].is_inside = true;
							if (i != (M - 1))
							{
								grid[i][j].h3 = dr[n - k - 2] / _M[n - k - 2];
							}
							grid[i][j].h4 = dr[n - k - 2] / _M[n - k - 2];
							grid[i][j].r = rs[n - k - 2] - (double(i - Ms[k - 1] - 1) / _M[n - k - 2] * dr[n - k - 2]) - delta;
						}
					}

					for (int l = 0; l < (k + 2); l++)  //轴向划分k+3个区域
					{
						if (j > 0 && j < Ns[0])
						{
							grid[i][j].h1 = dz[0] / _N[0];
							grid[i][j].h2 = grid[i][j].h1;
						}
						else if (j == Ns[0])
						{
							grid[i][j].h1 = dz[0] / _N[0];
							grid[i][j].h2 = dz[1] / _N[1];
						}
						else if (j > Ns[n - k + l - 2] && j < Ns[n - k + l - 1])
						{
							grid[i][j].h1 = dz[n - k + l - 1] / _N[n - k + l - 1];
							grid[i][j].h2 = grid[i][j].h1;
						}
						else if (j == Ns[n - k + l - 1] && j != (N - 1))  //若在两个划分区间中，且不在荧光屏上
						{
							grid[i][j].h1 = dz[n - k + l - 1] / _N[n - k + l - 1];
							grid[i][j].h2 = dz[n - k + l] / _N[n - k + l];
						}
					}
				}
			}
		}

		//设置光阴极和荧光屏
		for (int i = 0; i < M; i++)
		{
			grid[i][0].voltage = 0;
			grid[i][0].is_margin = true;
			grid[i][N - 1].voltage = screen_voltage;
			grid[i][N - 1].is_margin = true;
		}

		//设置轴上点
		for (int j = 1; j < (N - 1); j++)
		{
			grid[M - 1][j].on_axis = true;
			grid[M - 1][j].is_inside = true;
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
