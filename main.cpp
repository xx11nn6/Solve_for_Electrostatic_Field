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
	bool is_inside;  //�ж�������Ƿ��ڱ߽��ڲ�
	bool on_axis;  //�ж��Ƿ������ϣ�
	double h1, h2, h3, h4;  //�洢����㵽��Χ�ľ���
	double r;  //�洢�����ľ������
};


//����ȫ�ֱ���
int M, N;//����M*N
int mode = 1;  //���幤��ģʽ��1Ϊ��һ����ܣ�2Ϊ�ڶ������
//��������
void grid_initialize_mode_1(Grid_Array** grid, double screen_voltage, int n, double* dz, int* _N, double* V, double r1, double r2, int M1, int M2, double delta);
void grid_initialize_mode_2(Grid_Array** grid, double screen_voltage, int n, double* dz, double* dr, int* _N, int* _M, double* V, double delta);
double residual(Grid_Array** grid);
double SOR(Grid_Array** grid, double omega);
double select_accelerator_factor(Grid_Array** grid);
double convergence_criteria(Grid_Array** grid);
void free_grid(Grid_Array** grid);
void potential_line(Grid_Array** grid, int n, int* _N, int M1, int M2);
void paint(int k, Grid_Array** grid);
int scan(double V, Grid_Array** grid);
void paint_all(int n, int M1, int M2, int M, int N, Grid_Array** grid, double z0, double r0, double* dz, double delta);
double* V_r, * V_z; //��ɨ�������
double* z1, * r1;

int main()
{
	//��ʼ��
	mode = 2;
	///����Ϊ�������ͨ�ò���
	double V_s = 100;				//����ӫ������ѹΪ100V
	int n = 7;						//����缫��������ӫ������������������
	double delta = 0.5;				//����缫��Ȧ�
	double V[6];					//��������缫��ѹ
	int* _N = new int[n+mode-1];			//�缫֮���������Ҫ���ֵ������������»���Ϊ����ȫ�ֱ���N����
	double* dz = new double[n+mode-1];		//����ÿ���缫�������࣬������ܾ��õ��˱���						
	int _M[6];						//�缫֮�侶������Ҫ���ֵ����������ڶ������ʹ�ô˱���
	int M1, M2;					    //��ֱ����������ֵ�Ҫ�󣬵�һ�����ʹ�ô˱���				    
	double dr[6];					//����缫�����࣬�ڶ�������ô˱���
	//����Ϊ��һ������õ��ı���
	double z0 = 56.4, r0 = 32;		//���徶�������Ŀ��
	double r2 = 12, r1 = r0 - r2;   //r2Ϊ�缫����糡����ȣ�r1Ϊ�缫�׶˵���ľ���
	//////////////////ע�⣡��M*NΪ��������_N��_M�Ǹ���///////////////////////////////
	if (mode == 1)  //��һ����ܵĲ���
	{
		//����ֵ
		V[0] = 24; V[1] = 40; V[2] = 62; V[3] = 74; V[4] = 85; V[5] = 96;
		dz[0] = 5.2; dz[1] = 8.6; dz[2] = 8.6; dz[3] = 8.6; dz[4] = 8.6; dz[5] = 8.6; dz[6] = 5.2;
		_N[0] = 3; _N[1] = 5; _N[2] = 5; _N[3] = 5; _N[4] = 5; _N[5] = 5; _N[6] = 2;
		M1 = 11; M2 = 7;
		//����������M��N
		M = M1 + M2 + 1;  //������M1+M2,������Ҫ��һ
		N = _N[0] + _N[1] + _N[2] + _N[3] + _N[4] + _N[5] + _N[6] + n; //��_N+n-1��������������Ҫ��һ
		//grid_initialize_mode_1(grid, V_s, n, dz, _N, V, r1, r2, M1, M2, delta);
	}

	else if (mode == 2)  //�ڶ�����ܲ���
	{
		//����ֵ
		V[0] = 24; V[1] = 40; V[2] = 62; V[3] = 74; V[4] = 85; V[5] = 96;
		dr[0] = 5.2; dr[1] = 8.6; dr[2] = 8.6; dr[3] = 8.6; dr[4] = 8.6; dr[5] = 8.6;
		dz[0] = 5.2; dz[1] = 8.6; dz[2] = 8.6; dz[3] = 8.6; dz[4] = 8.6; dz[5] = 8.6; dz[6] = 8.6; dz[7] = 5.2;
		_N[0] = 3; _N[1] = 5; _N[2] = 5; _N[3] = 5; _N[4] = 5; _N[5] = 5; _N[6] = 5; _N[7] = 3;
		_M[0] = 3; _M[1] = 5; _M[2] = 5; _M[3] = 5; _M[4] = 5; _M[5] = 5;
		//����������
		M = _M[0] + _M[1] + _M[2] + _M[3] + _M[4] + _M[5] + n;
		N = _N[0] + _N[1] + _N[2] + _N[3] + _N[4] + _N[5] + _N[6] +_N[7]+ 1;
		//grid_initialize_mode_2(grid, V_s, n, dz, dr, _N, _M, V, delta);
	}

	double omega;  //����������Ӧ�
	double epsilon = 0.0005;  //��������Ϊ0.0005



	//new���ڶ�̬�����ڴ�
	//���磺int* p = new int[10];
	//��仰��ʾ����10��int���͵Ŀռ䣬�����׵�ַ�洢��p��
	//��˿�����������һ���ɱ��С������
	//Ҫ�봴����ά���飬�����ʹ��˫��ָ��
	Grid_Array** grid = new Grid_Array * [M];  //��һ��ָ��ָ����
	for (int i = 0; i < M; i++)
	{
		grid[i] = new Grid_Array[N];  //�ڶ���ָ����
	}

	grid_initialize_mode_2(grid, V_s, n, dz, dr, _N, _M, V, delta);
	
	//grid_initialize_mode_1(grid, V_s, n, dz, _N, V, r1, r2, M1, M2, delta);  //��ʼ����һ�����
	///*
	omega = select_accelerator_factor(grid);  //��������omegaѡȡ
	do
	{
		SOR(grid, omega);  //���е���ֱ�����Ͼ�������
	} while (convergence_criteria(grid) > epsilon);

	cout << "accelerator factor omega:" << endl;  //�����������
	cout << omega << endl;
	cout << "iteration times:" << endl;  //�����������
	cout << grid[M-2][N-2].k << endl;
	//*/
	//��������λ
	cout << "grid:" << endl;
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout << " " << setfill(' ') << setw(4) << setprecision(6) << grid[i][j].voltage_before << " ";  //setfill�Ⱥ�����iomanip���еĺ��������ڿ��������ʽ
			//cout << "(" << setprecision(2) << grid[i][j].h3 << "," << setprecision(2) << grid[i][j].h4 << ") ";
		}
		cout << endl;
	}
	//cout << sum(grid)/(M*N) << endl;
	//��ͼ
	//paint_all(n, M1, M2, M, N, grid, z0, r0, dz, delta);
	system("pause");
	free_grid(grid);//���н������ͷ��ڴ�
	return 0;
}




//���ڳ�ʼ����һ����ܵĵ糡����
void grid_initialize_mode_1(Grid_Array** grid, double screen_voltage, int n, double* dz, int* _N, double* V, double r1, double r2, int M1, int M2, double delta)
//����糡����ӫ������ѹ���缫����n�������ࣨ��������dz����ˮƽ�������񻮷�Ҫ�󣨴�������_N�����缫��ѹ������V��
//�缫�׵����r1���缫���r2����ֱ�������񻮷�Ҫ��M1,M2���缫���delta
{
	//����������г�ʼ��
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			grid[i][j].voltage = 0;  //��ʼ��λȫ������
			grid[i][j].is_margin = false;  //��ʼȫ������Ϊ�Ǳ߽�
			grid[i][j].is_inside = true;  //��һ����ܳ��߽磬����ȫ���������
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
				grid[i][j].voltage = 0;
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
								grid[i][j].voltage = V[0] * j / _N[0];
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
		grid[i][N - 1].voltage_before = screen_voltage;
	}

	//���ñ߽��ڲ�
	for (int i = 1; i < M ; i++)
	{
		for (int j = 1; j < (N - 1); j++)
		{
			grid[i][j].is_inside = true;
		}
	}
}



//���ڳ�ʼ���ڶ�����ܵĵ糡����
void grid_initialize_mode_2(Grid_Array** grid, double screen_voltage, int n, double* dz, double* dr, int* _N, int* _M, double* V, double delta)
//����糡����ӫ������ѹ���缫����n�����������ࣨ��������dz����������dr��ˮƽ�������񻮷�Ҫ�󣨴�������_N������ֱ�������񻮷�Ҫ��_M���缫��ѹ������V�����缫���delta
{
	//�ֿ� ����������洢����������Ҫ����ۼӺͣ������������
	int* Ms = new int[n-1];
	int* Ns = new int[n+1];
	double* rs = new double[n - 1];

	Ms[0] = _M[n - 2] + 1;  //������_M[0]�ģ�����û�뵽�ɹ�����Ҫ�����񻮷ֺ�����ϵ�Ƿ��ŵ�......
	Ns[0] = _N[0];
	rs[0] = dr[0] + delta;
	for (int i = 1; i < (n+1); i++)
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
	//�ȱ������񣬰����е�����
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			grid[i][j].voltage = 0.0;  //��ʼ��λȫ������
			grid[i][j].is_margin = false;  //��ʼȫ������Ϊ�Ǳ߽�
			grid[i][j].on_axis = false;  //��ʼ����Ϊ�����ϵ�
			grid[i][j].is_inside = false;//��ʼ����Ϊ���ڲ�
			grid[i][j].k = 0;
			grid[i][j].r = 0;
			grid[i][j].h1 = 0;
			grid[i][j].h2 = 0;
			grid[i][j].h3 = 0;
			grid[i][j].h4 = 0;
		}
	}

	//�ٴα��������г�ʼ���ͱ߽�����
	for (int k = 0; k < (n-1); k++)
		//�ֿ鴦��ÿһ���������缫��
	{
		//�������� ��һ���缫���ڶ���֮��
		if (k == 0)
		{
			for (int i = 0; i < Ms[k]; i++)
			{
				for (int j = 0; j < N; j++)
				{
					if (i == 0 || i == 1)  //����ǵ�һ�л��ǵڶ���
					{
						if ((j >= Ns[n - k - 2]) && (j <= Ns[n - k - 1]))  //����ڵ缫����
						{
							grid[i][j].voltage = V[n - k - 2];	//���õ缫��ѹ
							grid[i][j].is_margin = true;
						}
						else if (j > Ns[n - k - 1] && j < (N - 1))  //����ڵ缫��ӫ������
						{
							if (i == 0)  //���õ�һ��Ϊ�߽��ֵ
							{
								grid[i][j].voltage = V[n - k - 2] + double(j - Ns[n - k - 1]) / (Ns[n - k] - Ns[n - k - 1]) * (screen_voltage - V[n - k - 2]);  //�м����Բ�ֵ�����ѹ
								grid[i][j].is_margin = true;
							}
							else if (i == 1) //���õڶ���
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

					else if ((i > 1) && (i < Ms[k]))  //����ǵ����е���һ���缫��
					{
						if (j == Ns[n - k - 2])  //��������缫�䣬��ֱ����߽�
						{
							float r;  //�������������r
							r = rs[n - k - 2] - (double(i - 1) / _M[n - k - 2] * dr[n - k - 2]) - delta;
							grid[i][j].voltage = V[n - k - 2] + (V[n - k - 3] - V[n - k - 2]) * \
								log(r / (rs[n - k - 2] - delta)) / log(rs[n - k - 3] / rs[n - k - 2]);  //������ֵ������д����
							grid[i][j].is_margin = true;
							grid[i][j].r = r;
						}
						else if (j > Ns[n - k - 2] && j < (N - 1))  //���ñ߽��ڲ�
						{
							grid[i][j].is_inside = true;
							grid[i][j].h3 = dr[n - k - 2] / _M[n - k - 2];
							grid[i][j].h4 = grid[i][j].h3;
							grid[i][j].r = rs[n - k - 2] - (double(i - 1) / _M[n - k - 2] * dr[n - k - 2]) - delta;
							for (int l = 0; l < (k + 2); l++)  //���򻮷�k+2������
							{
								if (j > Ns[n - k + l - 2] && j < Ns[n - k + l - 1])
								{
									grid[i][j].h1 = dz[n - k + l - 1] / _N[n - k + l - 1];
									grid[i][j].h2 = grid[i][j].h1;
								}
								else if (j == Ns[n - k + l - 1] && j != (N - 1))  //�����������������У��Ҳ���ӫ������
								{
									grid[i][j].h1 = dz[n - k + l - 1 ] / _N[n - k + l - 1];
									grid[i][j].h2 = dz[n - k + l ] / _N[n - k + l ];
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

		//�ڶ�������n-1���缫
		else  if (k != (n - 2))
		{
			for (int i = Ms[k - 1]; i < Ms[k]; i++)
			{
				for (int j = Ns[n - k - 2]; j < (N - 1); j++)
				{
					if (i == Ms[k - 1])  //����ǵ缫��
					{
						if ((j >= Ns[n - k - 2]) && j <= (Ns[n - k - 1]))  //����ڵ缫����
						{
							grid[i][j].voltage = V[n - k - 2];	//���õ缫��ѹ
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
					else if (i == (Ms[k - 1] + 1))//����ǵ缫��һ��
					{
						if ((j >= Ns[n - k - 2]) && (j <= Ns[n - k - 1]))  //����ڵ缫����
						{
							grid[i][j].voltage = V[n - k - 2];	//���õ缫��ѹ
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
					else  //������ǵ缫����
					{
						if (j == Ns[n - k - 2])//��������缫�У���ֱ����߽�
						{
							{
								grid[i][j].voltage = V[n - k - 2];	//���õ缫��ѹ
								float r;  //�������������r
								r = rs[n - k - 2] - (double(i - Ms[k-1] - 1) / _M[n - k - 2] * dr[n - k - 2]) - delta;
								grid[i][j].voltage = V[n - k - 2] + (V[n - k - 3] - V[n - k - 2]) * \
									log(r / (rs[n - k - 2] - delta)) / log(rs[n - k - 3] / rs[n - k - 2]);  //������ֵ
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

					for (int l = 0; l < (k + 2); l++)  //���򻮷�k+2������
					{
						if (j > Ns[n - k + l - 2] && j < Ns[n - k + l - 1])
						{
							grid[i][j].h1 = dz[n - k + l - 1] / _N[n - k + l - 1];
							grid[i][j].h2 = grid[i][j].h1;
						}
						else if (j == Ns[n - k + l - 1] && j != (N - 1))  //�����������������У��Ҳ���ӫ������
						{
							grid[i][j].h1 = dz[n - k + l - 1] / _N[n - k + l - 1];
							grid[i][j].h2 = dz[n - k + l ] / _N[n - k + l ];
						}
					}
				}
			}
		}

		//��������һ���缫����
		else if (k == (n - 2))
		{
			for (int i = Ms[k - 1]; i <= Ms[k]; i++)
			{
				for (int j = 0; j < (N - 1); j++)
				{
					if (i == Ms[k - 1])  //����ǵ缫��
					{
						if ((j >= Ns[n - k - 2]) && j <= (Ns[n - k - 1]))  //����ڵ缫����
						{
							grid[i][j].voltage = V[n - k - 2];	//���õ缫��ѹ
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
					else if (i == (Ms[k - 1] + 1))//����ǵ缫��һ��
					{
						if (j > 0 && j < Ns[n - k - 2])	//���Բ�ֵ�����ѹ
						{
							grid[i][j].voltage = double(j) / _N[0] * V[0];
							grid[i][j].is_margin = true;
						}
						else if ((j >= Ns[n - k - 2]) && (j <= Ns[n - k - 1]))  //����ڵ缫����
						{
							grid[i][j].voltage = V[n - k - 2];	//���õ缫��ѹ
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
					else  //������ǵ缫����
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

					for (int l = 0; l < (k + 2); l++)  //���򻮷�k+3������
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
						else if (j == Ns[n - k + l - 1] && j != (N - 1))  //�����������������У��Ҳ���ӫ������
						{
							grid[i][j].h1 = dz[n - k + l - 1] / _N[n - k + l - 1];
							grid[i][j].h2 = dz[n - k + l] / _N[n - k + l];
						}
					}
				}
			}
		}

		//���ù�������ӫ����
		for (int i = 0; i < M; i++)  
		{
			grid[i][0].voltage = 0;
			grid[i][0].is_margin = true;
			grid[i][N-1].voltage = screen_voltage;
			grid[i][N-1].is_margin = true;
		}

		//�������ϵ�
		for (int j = 1; j < (N - 1); j++)
		{
			grid[M - 1][j].on_axis = true;
			grid[M - 1][j].is_inside = true;
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
			diff = abs(grid[i][j].voltage - grid[i][j].voltage_before);
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

//ʹ��delete���ͷ��ڴ�
void free_grid(Grid_Array** grid)
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