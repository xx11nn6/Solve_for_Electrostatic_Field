#include <iostream>
#include <malloc.h>
#include <numeric>
#include <iomanip>
using namespace std;

struct Grid_Array  //����糡����ṹ
{
	int k;  //�洢��������k
	double voltage;  //�洢����������k�εĵ�ѹ
	double voltage_before;  //�洢��һ�ε����ĵ�ѹ
	bool is_margin;  //�ж�������Ƿ��Ǳ߽磬1��ʾ�Ǳ߽�
	bool on_axis;  //�ж��Ƿ������ϣ�
	double h1, h2, h3, h4;  //�洢����㵽��Χ�ľ���
	double r;  //�洢�����ľ������
};


//����ȫ�ֱ���
int M, N;//������M*N
int mode = 1;  //���幤��ģʽ��1Ϊ��һ����ܣ�2Ϊ�ڶ������
//��������
void grid_initialize_mode_1(Grid_Array** grid, double cathode_voltage, double screen_voltage, int n, double* dz, int* _N, double* V, double r1, double r2, int M1, int M2, double delta);
double residual(Grid_Array** grid);
double SOR(Grid_Array** grid, double omega);
double select_accelerator_factor(Grid_Array** grid);
double convergence_criteria(Grid_Array** grid);
void free(Grid_Array** grid);

int main()
{
	//��ʼ��
	//////////////////ע�⣡��M*NΪ��������_N��M1��M2�Ǹ���///////////////////////////////
	double V_c = 0, V_s = 100;		//�����������ӫ������ѹ�ֱ�Ϊ0��100V
	int n = 7;						//����缫��������ӫ������������������
	double delta = 0.5;				//����缫��Ȧ�
	double z0 = 53.4, r0 = 32;		//���徶�������Ŀ��
	double r2 = 12, r1 = r0 - r2;   //r2Ϊ�缫����糡����ȣ�r1Ϊ�缫�׶˵���ľ���
	double dz[7];				    //����ÿ���缫�ļ��
	int _N[7];						//����ÿ���缫֮����ȡ�ĸ�����ˮƽ���򣩣����»���Ϊ����ȫ�ֱ���N����
	double V[6];					//��������缫��ѹ
	int M1, M2;					    //��ֱ����������ֵ�Ҫ��

	V[0] = 24; V[1] = 40; V[2] = 62; V[3] = 74; V[4] = 85; V[5] = 96;
	dz[0] = 5.2; dz[1] = 8.6; dz[2] = 8.6; dz[3] = 8.6; dz[4] = 8.6; dz[5] = 8.6; dz[6] = 5.2;
	_N[0] = 3; _N[1] = 5; _N[2] = 5; _N[3] = 5; _N[4] = 5; _N[5] = 5; _N[6] = 2;
	M1 = 11; M2 = 7;
	//����������M��N
	M = M1 + M2 + 1;  //������M1+M2,������Ҫ��һ
	N = _N[0] + _N[1] + _N[2] + _N[3] + _N[4] + _N[5] + _N[6] + n; //��_N+n-1��������������Ҫ��һ

	double omega;  //����������Ӧ�
	double epsilon = 0.005;  //��������Ϊ0.005



	//c++�У�new���ڶ�̬�����ڴ�
	//���磺int* p = new int[10];
	//��仰��ʾ����10��int���͵Ŀռ䣬�����׵�ַ�洢��p��
	//��˿�����������һ���ɱ��С������
	//Ҫ�봴����ά���飬�����ʹ��˫��ָ��
	Grid_Array** grid = new Grid_Array * [M];  //��һ��ָ��ָ����
	for (int i = 0; i < M; i++)
	{
		grid[i] = new Grid_Array[N];  //�ڶ���ָ����
	}
	grid_initialize_mode_1(grid, V_c, V_s, n, dz, _N, V, r1, r2, M1, M2, delta);
	omega = select_accelerator_factor(grid);
	do
	{
		SOR(grid, omega);
	} while (convergence_criteria(grid) > epsilon);
	cout << "iteration times:" << endl;
	cout << grid[5][5].k << endl;


	//c++����cout�������ʹ��cout�������
	cout << "grid:" << endl;
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout << " " << setfill(' ') << setw(4) << setprecision(6) << grid[i][j].voltage << " ";  //setfill�Ⱥ�����iomanip���еĺ��������ڿ��������ʽ
			//cout << "(" << setprecision(2) << grid[i][j].h3 << "," << setprecision(2) << grid[i][j].h4 << ") ";
		}
		cout << endl;
	}
	//cout << sum(grid)/(M*N) << endl;
	free(grid);//���н������ͷ��ڴ�

	return 0;
}




//���ڳ�ʼ����һ����ܵĵ糡����
void grid_initialize_mode_1(Grid_Array** grid, double cathode_voltage, double screen_voltage, int n, double* dz, int* _N, double* V, double r1, double r2, int M1, int M2, double delta)
//����糡���������ȡ��߶ȡ�������ѹ��ӫ������ѹ���缫����n�������ࣨ��������dz����ˮƽ�������񻮷�Ҫ�󣨴�������_N�����缫��ѹ������V��
//�缫�׵����r1���缫���r2����ֱ�������񻮷�Ҫ��M1,M2���缫���delta
{
	//����������г�ʼ��
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			grid[i][j].voltage = 0;  //��ʼ��λȫ������
			grid[i][j].is_margin = false;  //��ʼȫ������Ϊ�Ǳ߽�
			grid[i][j].on_axis = false;  //��ʼ����Ϊ�����ϵ�
			grid[i][j].k = 0;

			//���ô�ֱ���h3,h4
			//���õ�һ��
			if (i == 0)
			{
				grid[i][j].h4 = 0;
				grid[i][j].h3 = r2 / M2;
				grid[i][j].r = r1 + r2;
			}
			//����1~M2��
			else if (i < M2)
			{
				grid[i][j].h4 = r2 / M2;
				grid[i][j].h3 = grid[i][j].h4;
				grid[i][j].r = r1 + r2 - ((r2 * i) / M2);
			}
			//����M2��
			else if (i == M2)
			{
				grid[i][j].h4 = r2 / M2;
				grid[i][j].h3 = r1 / M1;
				grid[i][j].r = r1;
			}
			//����M2~M1��
			else if (i < (M2 + M1))
			{
				grid[i][j].h4 = r1 / M1;
				grid[i][j].h3 = grid[i][j].h4;
				grid[i][j].r = r1 - ((r1 * (i - M2)) / M1);
			}
			//�������ϵ�
			else
			{
				grid[i][j].h3 = 0;
				grid[i][j].h4 = r1 / M1;
				grid[i][j].r = 0;
				grid[i][j].on_axis = true;
			}

			//��һ�У���������
			if (j == 0)
			{
				grid[i][j].voltage = cathode_voltage;
				grid[i][j].h1 = 0;
				grid[i][j].h2 = dz[0] / (_N[0]);
				grid[i][j].is_margin = true;
			}
			else
			{
				//���漸����Ҫ��ѭ���ж�
				int sum = _N[0];
				int k = 0;
				while (k < n)
				{
					if (j < sum)
					{
						if (i == 0)  //���i=0����һ�У�ֱ������ձ߽紦��
						{
							grid[i][j].is_margin = true;
							if (k == 0)  //���������͵�һ���缫�䣬���Բ�ֵ����߽��ѹ
							{
								grid[i][j].voltage = (V[0] - cathode_voltage) * j / _N[0];
							}
							else if (k != (n - 1))  //���������һ���缫��ӫ��������
							{
								grid[i][j].voltage = V[k - 1] + (V[k] - V[k - 1]) * (_N[k] + j - sum) / _N[k];
							}
							else  //���һ���缫��ӫ��������
							{
								grid[i][j].voltage = V[k - 1] + (screen_voltage - V[k - 1]) * (_N[k] + j - sum) / _N[k];
							}
						}
						grid[i][j].h1 = dz[k] / _N[k];
						grid[i][j].h2 = dz[k] / _N[k];
						break;
					}
					else if (j == sum && (k != (n - 1))) //�ڵ缫����Ҳ���ӫ����
					{
						if (i <= M2)  //���߶��ڵ缫��ȷ�Χ��
						{
							grid[i][j].is_margin = true;
							grid[i][j].voltage = V[k];
						}
						grid[i][j].h1 = dz[k] / _N[k];
						grid[i][j].h2 = delta;
						break;
					}
					else if (j == (sum + 1) && (k != (n - 1)))  //�ڵ缫�Ҳ��Ҳ���ӫ����
					{
						if (i <= M2)  //���߶��ڵ缫��ȷ�Χ��
						{
							grid[i][j].is_margin = true;
							grid[i][j].voltage = V[k];
						}
						grid[i][j].h1 = delta;
						grid[i][j].h2 = dz[k + 1] / _N[k + 1];
						break;
					}

					else if (k == (n - 1)) //���һ��
					{
						if (j == sum) //��ӫ������
						{
							grid[i][j].h1 = dz[k] / _N[k];
							grid[i][j].h2 = 0;
							grid[i][j].is_margin = true;
							break;
						}
					}
					else  //��Խ����һ���缫
					{
						k += 1;
						sum = sum + _N[k] + 1;
						continue;
					}
				}
			}


			//����ʼ������1�ε�ѹ�������0����ͬ
			for (int i = 0; i < M; i++)
			{
				for (int j = 0; j < N; j++)
				{
					grid[i][j].voltage_before = grid[i][j].voltage;
				}
			}
		}
	}

	//����ӫ������λ�ͱ߽����
	for (int i = 0; i < M; i++)
	{
		grid[i][N - 1].voltage = screen_voltage;
		grid[i][N - 1].is_margin = true;
	}
}

//���ڼ���в�ľ�ֵ
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

//����һ��SOR����
double SOR(Grid_Array** grid, double omega)
{
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{

			if (grid[i][j].is_margin == false)//�ж��Ƿ��Ǳ߽磬���Ǳ߽���������
			{
				double c1, c2, c3, c4, c0;
				//Ϊ�˹�ʽ�ÿ������ñ����ݴ�
				double phi = grid[i][j].voltage;  //����ǰ�ĵ�ѹ��[k]
				double phi_after;  //������ĵ�λ��[k+1]
				double phi_avg;  //ƽ����ѹ
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
					phi_avg = (phi + phi_after) / 2;
					phi_after = phi + omega * (phi_avg - phi);
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

//c++��ʹ��delete���ͷ��ڴ�
void free(Grid_Array** grid)
{
	for (int i = 0; i < M; i++) {
		delete[] grid[i];
	}
	delete[] grid;
}
