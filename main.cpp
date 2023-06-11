#include <iostream>
#include <malloc.h>
#include <stdlib.h>
#include <numeric>
#include <iomanip>
#include <graphics.h>		// ���� EasyX ͼ�ο�
#include <conio.h>
#include "Grid_Array.h"			//��������ȫ�ֱ���M��N���ͷ������ڴ溯��
#include "Grid_Initialize.h"	//���ڳ�ʼ���������
#include "Iteration.h"			//����SOR�������������Ӽ��㣬�в����͵��������жϺ���
#include "File_Operation.h"		//�����ļ��Ķ�д
#include "Paint_Equipotential_Line.h"
using namespace std;

//����ȫ�ֱ���
//�������
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
double r1, r2;        // r2Ϊ�缫����糡����ȣ�r1Ϊ�缫�׶˵���ľ���(��һ�����ʹ��)
int* _M;              // �缫֮�侶������Ҫ���ֵ����������ڶ������ʹ�ô˱���
double* dr;           // �缫�ڿװ뾶(�ڶ���)
double epsilon;       // ��������(ͨ��)
int NST;              // �����ӡ�ռ��λʱ���������(ͨ��)
int INS;              // ���ϵ�λ���Ⱦ��ֵʱ������(ͨ��)
int* V1;              // Ҫ��ɨ��ȵ�λ�ߵĵ�λ������ߵ�λֵ(ͨ��)

//��������
double omega_1;		  // �������Ӧأ��ް���
double omega_2;
int iteration_times_1;  // �����������ް���
int iteration_times_2;  // �����������а���
bool round_p = false;   // �ж�һ�ֵ����Ƿ����
int round_n;		    // ��ǰ����
int times_n;			// ��ǰ����
Iteration_Process* head1;  //ͷ��㣬�ް���
Iteration_Process* head2;  //ͷ��㣬�а���

//���Ƶ�λ�߱���
double z0, r0;		   // �����ܳ�/��
double* V_r, * V_z;    // ��ɨ�������
double* zi, * ri;	   // ������/�������
int count1;		       // ��λ�߼�����
int tmp;			   // �����λ��Ҫ��

//��д�ļ�����
FILE* handle_read;	//���ļ��ľ��
FILE* handle_write; //д�ļ��ľ��
errno_t err;


int main()
{
	//���ļ���
	count1 = 0;
	if ((err = fopen_s(&handle_read, "test6.txt", "rt")) != 0)
	{
		printf("�Ҳ����ļ��������룺%d\n", err);
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

	//�ް�������
	Grid_Array** grid1 = new Grid_Array * [M];  //��һ��ָ��ָ����
	for (int i = 0; i < M; i++)
	{
		grid1[i] = new Grid_Array[N];  //�ڶ���ָ����
	}


	//�а�������
	Grid_Array** grid2 = new Grid_Array * [M];  //��һ��ָ��ָ����
	for (int i = 0; i < M; i++)
	{
		grid2[i] = new Grid_Array[N];  //�ڶ���ָ����
	}

	cout << "\ncount:" << count1 << endl;

	//��ʼ�����
	if (mode == 1)
	{
		grid_initialize_mode_1(grid1, V_s, n, dz, _N, V, r1, r2, M1, M2, delta);  //��ʼ����һ�����(�ް���)
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

	//�����ް������������������
	
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

	omega_1 = select_accelerator_factor(grid1, ite1,&iteration_times_1);  //��������omegaѡȡ
	do
	{
		SOR(grid1, omega_1, ite1, &iteration_times_1);  //���е���ֱ�����Ͼ�������
	} 
	while (ite1->max_res >= epsilon);


	//�а������������������
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

	omega_2 = select_accelerator_factor(grid2, ite2,&iteration_times_2);  //��������omegaѡȡ
	do
	{
		SOR(grid2, omega_2, ite2,&iteration_times_2);  //���е���ֱ�����Ͼ�������
	} while (convergence_criteria(grid2) > epsilon);	


	//���Ƶ�λ��
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
	
	

	cout << "accelerator factor omega:" << endl;  //�����������
	cout << omega_1 << endl;
	cout << "iteration times:" << endl;  //�����������
	cout << iteration_times_1 << endl;
	//*/
	//��������λ
	cout << "grid:" << endl;
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			;
			//cout << " " << setw(4)<< setfill(' ')  << setprecision(4) << grid[i][j].r << " ";  //setfill�Ⱥ�����iomanip���еĺ��������ڿ��������ʽ
			//cout << "(" << setprecision(2) << grid[i][j].h3 << "," << setprecision(2) << grid[i][j].h4 << ") ";
		}
		//cout << endl;
	}


	//д�ļ������
	if ((err = fopen_s(&handle_write, "���.res", "w")) != 0)
	{
		printf("�Ҳ����ļ��������룺%d\n", err);
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
	free_grid(grid1);//���н������ͷ��ڴ�
	free_grid(grid2);//���н������ͷ��ڴ�
	return 0;
}
