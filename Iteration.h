#define ITERATION_H  //����SOR�������������Ӽ��㣬�в����͵��������жϺ���
#include "Grid_Array.h"
#include <numeric>


//����в�ľ�ֵ
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



//����һ��SOR����
void SOR(Grid_Array** grid, double omega)
{
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (grid[i][j].is_margin == false && grid[i][j].is_inside == true)//�ж��Ƿ��Ǳ߽缰�Ƿ��ڱ߽��ڲ������Ǳ߽粢�����ڲ���������
			{

				double c1, c2, c3, c4, c0;
				//Ϊ�˹�ʽ�ÿ������ñ����ݴ�
				double phi = grid[i][j].voltage;  //����ǰ�ĵ�ѹ��[k]
				double phi_after;  //������ĵ�λ��[k+1]
				double r = grid[i][j].r;
				double h1 = grid[i][j].h1;
				double h2 = grid[i][j].h2;
				double h3 = grid[i][j].h3;
				double h4 = grid[i][j].h4;
				if (grid[i][j].on_axis == false)  //�����ϵ�ĵ���
				{
					c1 = 2 / (h1 * (h1 + h2));
					c2 = 2 / (h2 * (h1 + h2));
					c3 = (2 * r - h4) / (r * h3 * (h3 + h4));
					c4 = (2 * r + h3) / (r * h4 * (h3 + h4));
					c0 = c1 + c2 + c3 + c4;
					//�Ƚ���һ����ͨ�����¶�-������������
					phi_after = (c1 * grid[i][j - 1].voltage + c2 * grid[i][j + 1].voltage \
						+ c3 * grid[i + 1][j].voltage + c4 * grid[i - 1][j].voltage) / c0;
					//��SOR����:��[k+1]=��[k]+��(��[k+1]_bar-��[k])
					phi_after = (1 - omega) * phi + omega * phi_after;
					grid[i][j].voltage_before = phi;
					grid[i][j].voltage = phi_after;
					grid[i][j].k += 1;
				}
				else  //���ϵ�ĵ��������ʽ��ͬ��ֻ����ϵ����ͬ
				{
					c1 = 2 / (h1 * (h1 + h2));
					c2 = 2 / (h2 * (h1 + h2));
					c3 = 0;
					c4 = 4 / (h4 * h4);
					c0 = 2 * ((1 / (h1 * h2)) + (2 / (h4 * h4)));
					//����û��c3��
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



//��Ѽ������Ӧص�ѡ��
double select_accelerator_factor(Grid_Array** grid)
{
	double E_bar;
	double E_bar_after;
	double lambda = 0;  //��ʹ��=0
	double omega_l, mu_l, omega_m;  //��_�ˣ���_�ˣ���_m
	double omega_m_bar;  //�������ֵ�������ļ������Ӧ�_mƽ��ֵ

	double omega = 1;  //ȡomega=1����һ��
	SOR(grid, omega);
	omega = 1.375;
	do
	{
		lambda = 0;
		for (int i = 0; i < 12; i++)
		{
			E_bar = residual(grid);  //�����һ�βв�
			SOR(grid, omega);  //����һ��
			E_bar_after = residual(grid);  //�ٴμ���в�
			if (i > 8)//������ε���
			{
				lambda += (E_bar_after / E_bar);  //ȡ������Φ����
			}
		}
		lambda /= 3;  //����������ε����˾�ֵ
		mu_l = (lambda + omega - 1) / (sqrt(lambda) * omega);
		omega_l = 2 / (1 + sqrt(1 - (pow(mu_l, 2))));
		omega_m = 1.25 * omega_l - 0.5;
		omega_m_bar = (omega + omega_m) / 2;
		omega = omega_m;
		//cout << omega << endl;
	} while (abs((omega_m_bar - omega_m) / (2 - omega_m)) >= 0.05);
	return omega;
}



//���������ж�
double convergence_criteria(Grid_Array** grid)
{
	double max_error, temp;
	max_error = 0;
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (grid[i][j].is_margin == false)//�ԷǱ߽���ж�
			{
				temp = abs(grid[i][j].voltage - grid[i][j].voltage_before);//�����������һ�ε������
				if (temp > max_error)
				{
					max_error = temp;  //ȡ��������е����ֵ
				}
			}
		}
	}
	return max_error;
}