#include <iostream>
#include <malloc.h>
#include <numeric>
#include <iomanip>
#include <graphics.h>		// ���� EasyX ͼ�ο�
#include <conio.h>

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
void potential_line(Grid_Array** grid, int n, int* _N, int M1, int M2);
void paint(int k, Grid_Array** grid);
int scan(double V, Grid_Array** grid);
void paint_all(int n, int M1, int M2, int M, int N, Grid_Array** grid, double z0, double r0, double* dz, double delta);
double* V_r, * V_z; //��ɨ�������
double* z1, * r1;

int main()
{
	//��ʼ��
	//////////////////ע�⣡��M*NΪ��������_N��M1��M2�Ǹ���///////////////////////////////
	double V_c = 0, V_s = 100;		//�����������ӫ������ѹ�ֱ�Ϊ0��100V
	int n = 7;						//����缫��������ӫ������������������
	double delta = 0.5;				//����缫��Ȧ�
	double z0 = 56.4, r0 = 32;		//���徶�������Ŀ��
	double r2 = 12, r1 = r0 - r2;   //r2Ϊ�缫����糡����ȣ�r1Ϊ�缫�׶˵���ľ���
	double dz[7];				    //����ÿ���缫�ļ��
	int _N[7];						//����ÿ���缫֮����ȡ�ĸ�����ˮƽ���򣩣����»���Ϊ����ȫ�ֱ���N����
	double V[6];					//��������缫��ѹ
	int M1, M2;					    //��ֱ����������ֵ�Ҫ��
	double delta_V;                 //������λ���
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
	paint_all(n, M1, M2, M, N, grid, z0, r0, dz, delta);
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

void paint_all(int n, int M1, int M2, int M, int N, Grid_Array** grid, double z0, double r0, double* dz, double delta)
//����λ��
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
		z1[i] = z1[i - 1] + grid[0][i - 1].h2;  //�ۼ���ͣ�����ÿһ����ԭ���������룬��Ч��matlab��cumsum
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
	initgraph(10 * z0, 10 * r0);  //��ʼ����ͼ���ڣ���Ϊz0*10����Ϊr0*10
	setorigin(0, 10 * r0);  //����ԭ��Ϊ��0��r0*10��
	setaspectratio(1, -1);  //���õ�ǰ��������,��תy�ᣬʹy����Ϊ��

	t = dz[0];  //t��ʾ�������ڶ�λ�缫λ��
	for (int i = 1; i < n; i = i + 1)  //����ɫ�߻��Ƶ缫��벿��
	{
		setcolor(RGB(0, 0, 255));
		line(t * 10, r0 * 10, t * 10, grid[M2][0].r * 10);  //line(x1,y1,x2,y2)������������㣬�յ�
		t = t + dz[i] + delta;
	}

	//����ɫ�߻��Ƶ缫�Ұ벿��
	t = dz[0] + delta;
	for (int i = 1; i < n; i = i + 1)
	{
		setcolor(RGB(0, 0, 255));
		line(t * 10, r0 * 10, t * 10, grid[M2][0].r * 10);
		t = t + dz[i] + delta;
	}

	//����ɫ�����ӵ缫
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
//ɨ�裬�����ɨ���λ���Լ�����grid
{
	int i, j, k = 0;
	double r, z;
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N - 1; j++)
		{
			if ((grid[i][j].voltage <= V && grid[i][j + 1].voltage > V) || (grid[i][j].voltage >= V && grid[i][j + 1].voltage < V))  //�жϴ�ɨ���λ��λ������
			{
				z = z1[j] + grid[0][j].h2 * (V - grid[i][j].voltage) / (grid[i][j + 1].voltage - grid[i][j].voltage);  //��ֵ�������ɨ���λλ��
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

void paint(int k, Grid_Array** grid) //����//
{
	int i, j;
	double temp, x, y, xy;
	x = grid[0][0].h2;
	for (i = 0; i < N - 1; i++)  //�ж���������Сֵ
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
		for (j = 1; j <= k - i; j++)  //ð������
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