#include <iostream>
#include <malloc.h>
#include<stdlib.h>
#include <numeric>
#include <iomanip>
#include <graphics.h>		// ���� EasyX ͼ�ο�
#include <conio.h>
#include "Grid_Array.h"			//��������ȫ�ֱ���M��N���ͷ������ڴ溯��
#include "Grid_Initialize.h"	//���ڳ�ʼ���������
#include "Iteration.h"			//����SOR�������������Ӽ��㣬�в����͵��������жϺ���
#include "File_Operation.h"		//�����ļ��Ķ�д
using namespace std;

//����ȫ�ֱ���
void potential_line(Grid_Array** grid, int n, int* _N, int M1, int M2);
void paint(int k, Grid_Array** grid);
int scan(double V, Grid_Array** grid);
void paint_all(int n, int M1, int M2, int M, int N, Grid_Array** grid, double z0, double r0, double* dz, double delta);
double* V_r, * V_z; //��ɨ�������
double* zi, * ri;

int mode;			 //����ģʽ1Ϊһ�࣬2Ϊ����
int M;			  	 //������MN
int N;
double delta;		 // �缫���
int n;				 // �缫��������ӫ������������������
double* dz;           // ÿ���缫�������࣬������ܾ��õ��˱���
int* _N;              // ���ڵ缫�仮�ֲ���(ͨ��)
double* V;            // �����缫��ѹ
double* VI;           // ������缫��λ(ͨ��)
int M1, M2;           // ��ֱ����������ֵ�Ҫ�󣬵�һ�����ʹ�ô˱���
double r1, r2;         // r2Ϊ�缫����糡����ȣ�r1Ϊ�缫�׶˵���ľ���(��һ�����ʹ��)
int* _M;              // �缫֮�侶������Ҫ���ֵ����������ڶ������ʹ�ô˱���
double* dr;           // �缫�ڿװ뾶(�ڶ���)
double epsilon;       // ��������(ͨ��)
double omega;		 // �������Ӧ�
int NST;              // �����ӡ�ռ��λʱ���������(ͨ��)
int INS;              // ���ϵ�λ���Ⱦ��ֵʱ������(ͨ��)
int* V1;              // Ҫ��ɨ��ȵ�λ�ߵĵ�λ������ߵ�λֵ(ͨ��)
int i;                // ѭ���ò���
FILE* handle;
errno_t err;
// һ�����������ж�ɨ���λ�ߵ�����




int main()
{
#if 0
	//��ʼ��
	mode = 2;
	double V_s;		//ӫ������ѹ
	//����Ϊ��һ������õ��ı���
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
	}


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
#endif
	
	//���ļ���
	if ((err = fopen_s(&handle, "test2.txt", "rt")) != 0)
	{
		printf("�Ҳ����ļ��������룺%d\n", err);
		return -1;
	}


	fopen_s(&handle, "test2.txt", "rt");
	if (fscanf_s(handle, "%d\n", &mode) == 1 && mode == 2)
	{
		readdata2("test2.txt");
		fclose(handle);
	}
	else
	{
		readdata1("test2.txt");
		fclose(handle);
	}

	double V_s = V[n - 1];  //����ӫ������ѹ
	//����M��N
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

	//��������grid
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

	//��ʼ�����
	if (mode == 1)
	{
		grid_initialize_mode_1(grid, V_s, n, dz, _N, V, r1, r2, M1, M2, delta);  //��ʼ����һ�����
	}
	else
	{
		grid_initialize_mode_2(grid, V_s, n, dz, dr, _N, _M, V, delta);
	}

	omega = select_accelerator_factor(grid);  //��������omegaѡȡ
	do
	{
		SOR(grid, omega);  //���е���ֱ�����Ͼ�������
	} while (convergence_criteria(grid) > epsilon);

	cout << "accelerator factor omega:" << endl;  //�����������
	cout << omega << endl;
	cout << "iteration times:" << endl;  //�����������
	cout << grid[M - 2][N - 2].k << endl;
	//*/
	//��������λ
	cout << "grid:" << endl;
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout << " " << setw(4)<< setfill(' ')  << setprecision(4) << grid[i][j].voltage_before << " ";  //setfill�Ⱥ�����iomanip���еĺ��������ڿ��������ʽ
			//cout << "(" << setprecision(2) << grid[i][j].h3 << "," << setprecision(2) << grid[i][j].h4 << ") ";
		}
		cout << endl;
	}

	free(dz);
	free(_N);
	free(V);
	free(VI);
	free(dr);
	free(_M);
	free(V1);

	return 0;
}




void paint_all(int n, int M1, int M2, int M, int N, Grid_Array** grid, double z0, double r0, double* dz, double delta)
//����λ��
{
	int i, j, k, num;
	double z_temp1, z_temp2, V_temp1, V_temp2, t;
	V_z = (double*)malloc(n * (M1 + M2) * sizeof(double));
	V_r = (double*)malloc(n * (M1 + M2) * sizeof(double));
	zi = (double*)malloc(N * sizeof(double));
	ri = (double*)malloc(M * sizeof(double));
	zi[0] = 0;
	for (i = 1; i < N; i++)
	{
		zi[i] = zi[i - 1] + grid[0][i - 1].h2;  //�ۼ���ͣ�����ÿһ����ԭ���������룬��Ч��matlab��cumsum
	}
	ri[0] = 0;
	for (i = 1; i < M; i++)
	{
		ri[i] = grid[i][0].r;
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
				z = zi[j] + grid[0][j].h2 * (V - grid[i][j].voltage) / (grid[i][j + 1].voltage - grid[i][j].voltage);  //��ֵ�������ɨ���λλ��
				V_z[k] = 200 + 10 * z;
				V_r[k] = 310 + 10 * (ri[M - 1] - ri[i]);
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
				r = ri[i] + grid[i][0].h4 * (V - grid[i][j].voltage) / (grid[i + 1][j].voltage - grid[i][j].voltage);
				V_z[k] = 200 + 10 * zi[j];
				V_r[k] = 310 + 10 * (ri[M - 1] - r);
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