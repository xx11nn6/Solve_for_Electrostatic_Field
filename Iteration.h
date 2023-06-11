#define ITERATION_H  //包括SOR迭代，加速因子计算，残差计算和迭代精度判断函数
#include "Grid_Array.h"
#include <numeric>


//计算残差的均值
double residual(Grid_Array** grid)
{
	double diff;
	double res;
	double sum = 0;
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
void SOR(Grid_Array** grid, double omega)
{
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (grid[i][j].is_margin == false && grid[i][j].is_inside == true)//判断是否是边界及是否在边界内部，不是边界并且在内部则参与迭代
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
			else
			{
				grid[i][j].voltage = grid[i][j].voltage_before;
			}
		}
	}
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