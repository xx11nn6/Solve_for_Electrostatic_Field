#define PAINT_EQUIPOTENTIAL_LINE_H
//���ڻ��Ƶ�λ��
#include "Grid_Array.h"
#include <graphics.h>
#include <iostream>
using namespace std;


void potential_line(Grid_Array** grid, int n, int* _N, int M1, int M2);
void paint_1(int k, Grid_Array** grid, double r0);
int scan_1(double V, Grid_Array** grid);
void paint_all_1(int n, int M1, int M2, int M, int N, Grid_Array** grid, double z0, double r0, double* dz, double delta, int tmp, int* V1, int count1);
void paint_2(int k, Grid_Array** grid, int* _N, int* _M, double* dz, double* dr);
int scan_2(int V, Grid_Array** grid, int n, double* dz, double* dr, double delta);
void paint_all_2(int n, int M, int N, Grid_Array** grid, double* dz, double delta, int* _M, double* dr, int* _N, int tmp, int* V1, int count1);
void axis_voltage(double z_sum, double* z, int INS);

void paint_all_1(int n, int M1, int M2, int M, int N, Grid_Array** grid, double z0, double r0, double* dz, double delta, int tmp, int* V1, int count1)//����λ��,
{
	int i, j, k, num;
	double z_1, v_1, z_0, v_0, V_temp1, V_temp2, t;
	V_z = (double*)malloc(n * (M1 + M2) * sizeof(double));//��ɨ���Z����
	V_r = (double*)malloc(n * (M1 + M2) * sizeof(double));//��ɨ���R����
	zi = (double*)malloc(N * sizeof(double));//��ÿ�е�r�����
	ri = (double*)malloc(M * sizeof(double));//��ÿ�о������
	zi[0] = 0;
	for (i = 1; i < N; i++)//����z1
	{
		zi[i] = zi[i - 1] + grid[0][i - 1].h2;
	}
	for (i = 0; i < M; i++)//����
	{
		ri[i] = grid[i][0].r;
	}

	for (i = 0; i < n * (M1 + M2); i++)//������ɨ���������ʼ��Ϊ0
	{
		V_z[i] = 0;
		V_r[i] = 0;
	}
	initgraph(10 * z0, 10 * r0);
	setorigin(0, 10 * r0);
	setaspectratio(1, -1);

	t = dz[0];
	for (int i = 1; i < n; i = i + 1)//���缫��߽�
	{
		setcolor(RGB(0, 0, 255));
		line(t * 10, r0 * 10, t * 10, grid[M2][0].r * 10);
		t = t + dz[i] + delta;
	}
	t = dz[0] + delta;
	for (int i = 1; i < n; i = i + 1)//���缫���ұ߽�
	{
		setcolor(RGB(0, 0, 255));
		line(t * 10, r0 * 10, t * 10, grid[M2][0].r * 10);
		t = t + dz[i] + delta;
	}
	t = dz[0];
	for (int i = 1; i < n; i = i + 1)//���缫���±߽�
	{
		setcolor(RGB(0, 0, 255));
		line(t * 10, grid[M2][0].r * 10, (t + delta) * 10, grid[M2][0].r * 10);
		t = t + dz[i] + delta;
	}

	if (tmp == 3) {
		for (i = 0; i < 100; i = i + V1[0]) {
			setcolor(RGB(rand() * 255, rand() * 255, rand() * 255));
			k = scan_1(i, grid);//����Ҫɨ��ĵ�λ������ָ����λ��ɨ�裬����kΪɨ�������
			paint_1(k, grid, r0);
		}
	}
	else {

		for (i = 0; i < count1; i++) {
			setcolor(RGB(rand() * 255, rand() * 255, rand() * 255));
			k = scan_1(V1[i], grid);//����Ҫɨ��ĵ�λ������ָ����λ��ɨ�裬����kΪɨ�������
			paint_1(k, grid, r0);
		}



	}
	system("pause");
	closegraph();
	//for (j = 0; j < k; j++) { //��ͼ�ϵ��ɨ���
		//putpixel(V_z[j], V_r[j], 0xffffff);
	initgraph(10 * z0, 10 * r0);
	setorigin(0, 10 * r0);
	setaspectratio(1, -1);

	double z_pace = z0 / INS;
	z_0 = 0;
	v_0 = 0;

	//cout << z_0 << "��ѹ" << v_0 << endl;
	line(100, 100, 100 + z0 * 3, 100);
	line(100, 100, 100, 100 + 300);

	for (i = 0; i < INS - 1; i++) {
		z_1 = z_0 + z_pace;


		for (j = 0; j < (N - 1); j++) {
			if ((zi[j] <= z_1) && (z_1 < zi[j + 1])) {
				v_1 = (z_1 - zi[j]) * (z_1 - zi[j + 1]) / (zi[j - 1] - zi[j]) / (zi[j - 1] - zi[j + 1]) * grid[M - 1][j - 1].voltage\
					+ (z_1 - zi[j - 1]) * (z_1 - zi[j + 1]) / (zi[j] - zi[j - 1]) / (zi[j] - zi[j + 1]) * grid[M - 1][j].voltage\
					+ (z_1 - zi[j - 1]) * (z_1 - zi[j]) / (zi[j + 1] - zi[j - 1]) / (zi[j + 1] - zi[j]) * grid[M - 1][j + 1].voltage;
				//cout << z_1 << "��ѹ" << v_1 << endl;
				line(100 + 3 * z_0, 100 + 3 * v_0, 100 + 3 * z_1, 100 + 3 * v_1);
				z_0 = z_1;
				v_0 = v_1;
			}
		}

	}
	//cout << z0 << "��ѹ" << 100 << endl;


	system("pause");
	closegraph();
}

int scan_1(double V, Grid_Array** grid)//ɨ��//
{
	int i, j, k = 0;
	double r, z;
	//cout << V << "V��λɨ�������(�а���)��" << endl;
	for (i = 0; i < M; i++) //��ɨ��
	{
		for (j = 0; j < N - 1; j++)
		{

			if (grid[i][j].voltage == V) {//�����λ��ɨ���λ��ȣ����õ�����걣�棬10Ϊͼ��ķŴ���
				V_z[k] = 10 * zi[j];
				V_r[k] = 10 * ri[i];
				k++;
				//printf(" ( %4.2f , %4.2f )", zi[i], ri[i]);
				
			}

			if ((grid[i][j].voltage <= V && grid[i][j + 1].voltage > V) || (grid[i][j].voltage >= V && grid[i][j + 1].voltage < V) || (grid[i][j].voltage > V && grid[i][j + 1].voltage <= V)\
				|| (grid[i][j].voltage < V && grid[i][j + 1].voltage >= V))//��ɨ������һ��ĵ�λ����ȣ�����в�ֵ
			{
				z = (V - grid[i][j].voltage) / (grid[i][j + 1].voltage - grid[i][j].voltage) * (zi[j + 1] - zi[j]) + zi[j];//���Բ�ֵ
				V_z[k] = 10 * z;
				V_r[k] = 10 * ri[i];
				k++;
				//printf(" ( %4.2f , %4.2f )", z, ri[i]);
				//cout << " " << setw(4) << setfill(' ') << setprecision(4) << "(" << z << "," << ri[i] << ")" << " ";

			}
		}
	}
	for (j = 0; j < N; j++)//��ɨ�裬���岽��ͬ��
	{
		for (i = 0; i < M - 1; i++)
		{


			if (grid[i][j].voltage == V) {//�����λ��ɨ���λ��ȣ����õ�����걣�棬10Ϊͼ��ķŴ���
				V_z[k] = 10 * zi[j];
				V_r[k] = 10 * ri[i];
				k++;
				//printf(" ( %4.2f , %4.2f )", zi[i], ri[i]);
				//cout << " " << setw(4) << setfill(' ') << setprecision(4) << "(" << zi[i] << "," << ri[i] << ")" << " ";
			}
			if ((grid[i][j].voltage <= V && grid[i + 1][j].voltage > V) || (grid[i][j].voltage >= V && grid[i + 1][j].voltage < V) || (grid[i][j].voltage > V && grid[i + 1][j].voltage <= V)\
				|| (grid[i][j].voltage < V && grid[i + 1][j].voltage >= V))//��ɨ������һ��ĵ�λ����ȣ�����в�ֵ
			{


				r = (V - grid[i][j].voltage) / (grid[i + 1][j].voltage - grid[i][j].voltage) * (ri[i + 1] - ri[i]) + ri[i];//���Բ�ֵ
				V_z[k] = 10 * zi[j];
				V_r[k] = 10 * r;
				k++;
				//printf(" ( %4.2f , %4.2f )", zi[i], r);
				//cout << " " << setw(4) << setfill(' ') << setprecision(4) << "(" << zi[i] << "," << r << ")" << " ";

			}
		}
	}
	cout << endl;
	return k;//������ɨ�赽���������
}

void paint_1(int k, Grid_Array** grid, double r0) //����//
{
	int i, j;
	double temp, x, y, xy, minlengthx_1, minlengthy_1, minlengthx_2, minlengthy_2, minlength;
	x = grid[0][0].h2;
	for (i = 0; i < N - 1; i++)//�ҳ�����z���
	{
		if (x < grid[0][i].h2)
			x = grid[0][i].h2;
	}
	y = grid[0][0].h4;
	for (i = 0; i < M - 1; i++)//�ҳ�����R���
	{
		if (y < grid[i][0].h4)
			y = grid[i][0].h4;
	}

	xy = sqrt(pow(x, 2) + pow(y, 2));//����֮����ܵ�������
	for (i = 1; i < k; i++)//��ɨ��ĵ�����������z������������
	{
		for (j = 1; j <= k - i; j++)
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
	for (i = 0; i < k; i++)  //�ӵ�һ���㿪ʼ��ÿ�����ҳ�����������������������֮��֮�ڣ�z��r���������z��r������δ���ӹ��ĵ�������ӣ�����ԭ�����Ǻ������
	{
		for (j = i + 1; j < k; j++)
		{
			if (sqrt(pow(V_r[j] - V_r[i], 2) + pow(V_z[j] - V_z[i], 2)) <= 15 * xy &&
				fabs(V_r[j] - V_r[i]) <= 15 * y && fabs(V_z[j] - V_z[i]) <= 15 * x)
			{
				line(int(V_z[i]), int(V_r[i]), int(V_z[j]), int(V_r[j]));
				break;
			}
		}
	}
	for (i = 1; i < k; i++)//��ɨ��ĵ�����������z���꽵������
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
	for (i = 0; i < k; i++)//�ӵ�һ���㿪ʼ��ÿ�����ҳ�����������������������֮��֮�ڣ�z��r���������z��r������δ���ӹ��ĵ�������ӣ�����ԭ�����Ǻ������
	{
		for (j = i + 1; j < k; j++)
		{
			if (sqrt(pow(V_r[j] - V_r[i], 2) + pow(V_z[j] - V_z[i], 2)) <= 15 * xy &&
				fabs(V_r[j] - V_r[i]) <= 15 * y && fabs(V_z[j] - V_z[i]) <= 15 * x)
			{
				line(int(V_z[i]), int(V_r[i]), int(V_z[j]), int(V_r[j]));
				break;
			}
		}
	}
}

void paint_all_2(int n, int M, int N, Grid_Array** grid, double* dz, double delta, int* _M, double* dr, int* _N, int tmp, int* V1, int count1)//����λ��,
{
	int i, j, k, num;
	FILE* file;
	double z_temp1, z_temp2, V_temp1, V_temp2, r_0, z_0, t, s, v_2, z_2, z_1, v_1;
	V_z = (double*)malloc(n * (M) * sizeof(double));//��ɨ���Z����
	V_r = (double*)malloc(n * (M) * sizeof(double));//��ɨ���R����
	zi = (double*)malloc(N * sizeof(double));//��ÿ�е�r�����
	ri = (double*)malloc(M * sizeof(double));//��ÿ�о������
	zi[0] = 0;
	grid[M - 1][0].h2 = dz[0] / _N[0];
	for (i = 1; i < N; i++)//����z1
	{
		zi[i] = zi[i - 1] + grid[M - 1][i - 1].h2;
	}
	for (i = 0; i < M; i++)//����r1
	{
		ri[i] = grid[i][N - 2].r;
	}
	for (i = 0; i < n * (M); i++)//������ɨ���������ʼ��Ϊ0
	{
		V_z[i] = 0;
		V_r[i] = 0;
	}
	r_0 = 0;
	for (i = 0; i < n - 1; i++) {
		r_0 = r_0 + dr[i] + delta;
	}
	z_0 = 0;
	for (i = 0; i < n + 1; i++) {
		z_0 = z_0 + dz[i];
	}
	initgraph(8 * z_0, 8 * r_0);
	setorigin(0, 8 * r_0);
	setaspectratio(1, -1);
	t = dz[0];
	s = dr[0];
	for (i = 1; i < n; i = i + 1)//���缫�±߽�
	{
		setcolor(RGB(0, 0, 255));
		line(t * 8, s * 8, (t + dz[i]) * 8, s * 8);
		t = t + dz[i];
		s = s + dr[i] + delta;
	}
	t = dz[0];
	s = dr[0] + delta;
	for (i = 1; i < n; i = i + 1)//���缫�ϱ߽�
	{
		setcolor(RGB(0, 0, 255));
		line(t * 8, s * 8, (t + dz[i]) * 8, s * 8);
		t = t + dz[i];
		s = s + dr[i] + delta;
	}
	t = dz[0];
	s = dr[0];
	for (i = 1; i < n; i = i + 1)//���缫��߽�
	{
		setcolor(RGB(0, 0, 255));
		line(t * 8, s * 8, t * 8, (s + delta) * 8);
		t = t + dz[i];
		s = s + dr[i] + delta;
	}
	t = dz[0];
	s = dr[0];
	for (i = 1; i < n; i = i + 1)//���缫�ұ߽�
	{
		setcolor(RGB(0, 0, 255));
		line((t + dz[i]) * 8, s * 8, (t + dz[i]) * 8, (s + delta) * 8);
		t = t + dz[i];
		s = s + dr[i] + delta;
	}
	if (tmp == 3) {
		for (i = V1[0]; i < 100; i = i + V1[0]) {
			setcolor(RGB(rand() * 255, rand() * 255, rand() * 255));
			k = scan_2(i, grid, n, dz, dr, delta);//����Ҫɨ��ĵ�λ������ָ����λ��ɨ�裬����kΪɨ�������
			paint_2(k, grid, _N, _M, dz, dr);

		}
	}
	else {

		for (i = 0; i < count1; i++) {
			setcolor(RGB(rand() * 255, rand() * 255, rand() * 255));
			k = scan_2(V1[i], grid, n, dz, dr, delta);//����Ҫɨ��ĵ�λ������ָ����λ��ɨ�裬����kΪɨ�������
			paint_2(k, grid, _N, _M, dz, dr);
		}
	}
	system("pause");
	closegraph();

	initgraph(8 * z_0, 8 * r_0);
	setorigin(0, 8 * r_0);
	setaspectratio(1, -1);

	double z_pace = z_0 / INS;
	z_2 = 0;
	v_2 = 0;

	cout << z_2 << "��ѹ" << v_2 << endl;
	line(10, 10, 10 + z_0 * 2, 10);
	line(10, 10, 10, 10 + 200);

	fopen_s(&file, "D:\\�����ļ�.res", "w");
	for (i = 0; i < INS - 1; i++) {
		z_1 = z_2 + z_pace;


		for (j = 0; j < (N - 1); j++) {
			if ((zi[j] <= z_1) && (z_1 < zi[j + 1])) {
				v_1 = (z_1 - zi[j]) * (z_1 - zi[j + 1]) / (zi[j - 1] - zi[j]) / (zi[j - 1] - zi[j + 1]) * grid[M - 1][j - 1].voltage\
					+ (z_1 - zi[j - 1]) * (z_1 - zi[j + 1]) / (zi[j] - zi[j - 1]) / (zi[j] - zi[j + 1]) * grid[M - 1][j].voltage\
					+ (z_1 - zi[j - 1]) * (z_1 - zi[j]) / (zi[j + 1] - zi[j - 1]) / (zi[j + 1] - zi[j]) * grid[M - 1][j + 1].voltage;
				//cout << "���꣺" << z_1 << "��ѹ��" << v_1 << endl;
				//fprintf(file, "����:%lf ��ѹ��%lf", z_1, v_1);
				line(10 + 2 * z_2, 10 + 2 * v_2, 10 + 2 * z_1, 10 + 2 * v_1);
				z_2 = z_1;
				v_2 = v_1;
			}
		}

	}
	//cout << z_0 << "��ѹ" << 100 << endl;
	system("pause");
	closegraph();

}

int scan_2(int V, Grid_Array** grid, int n, double* dz, double* dr, double delta)//ɨ��//
{
	int i, j, k = 0;
	int p, q;
	double r, z, t, s;
	//cout << V << "V��λɨ������꣨�ް��㣩��" << endl;
	for (i = 0; i < M; i++) //��ɨ��
	{
		for (j = 0; j < N - 1; j++)
		{
			if (grid[i][j].voltage == V) {//�����λ��ɨ���λ��ȣ����õ�����걣�棬10Ϊͼ��ķŴ���
				V_z[k] = zi[j];
				V_r[k] = ri[i];
				k++;

			}
			if ((grid[i][j].voltage <= V && grid[i][j + 1].voltage > V) || (grid[i][j].voltage >= V && grid[i][j + 1].voltage < V) || (grid[i][j].voltage > V && grid[i][j + 1].voltage <= V)\
				|| (grid[i][j].voltage < V && grid[i][j + 1].voltage >= V))//��ɨ������һ��ĵ�λ����ȣ�����в�ֵ
			{
				z = (V - grid[i][j].voltage) / (grid[i][j + 1].voltage - grid[i][j].voltage) * (zi[j + 1] - zi[j]) + zi[j];//���Բ�ֵ
				V_z[k] = z;
				V_r[k] = ri[i];
				k++;

			}
		}
	}
	for (j = 0; j < N; j++)//��ɨ�裬���岽��ͬ��
	{
		for (i = 0; i < M - 1; i++)
		{
			if (grid[i][j].voltage == V) {//�����λ��ɨ���λ��ȣ����õ�����걣�棬10Ϊͼ��ķŴ���
				V_z[k] = zi[j];
				V_r[k] = ri[i];
				k++;
			}
			if ((grid[i][j].voltage <= V && grid[i + 1][j].voltage > V) || (grid[i][j].voltage >= V && grid[i + 1][j].voltage < V) || (grid[i][j].voltage > V && grid[i + 1][j].voltage <= V)\
				|| (grid[i][j].voltage < V && grid[i + 1][j].voltage >= V))//��ɨ������һ��ĵ�λ����ȣ�����в�ֵ
			{
				r = (V - grid[i][j].voltage) / (grid[i + 1][j].voltage - grid[i][j].voltage) * (ri[i + 1] - ri[i]) + ri[i];//���Բ�ֵ
				V_z[k] = zi[j];
				V_r[k] = r;

				k++;
			}
		}
	}
	q = 0;
	for (i = 0; i < k; i++) {

		t = 0;
		s = 0;
		for (p = 0; p < n + 1; p++) {

			if (p == 0) {
				if ((V_z[i] <= dz[0]) && (V_r[i] <= dr[0]))
				{
					V_z[q] = 8 * V_z[i];
					V_r[q] = 8 * V_r[i];
					q++;
					//printf(" ( %4.2f , %4.2f )", V_z[i], V_r[i]);
					//cout << " " << setw(4) << setfill(' ') << setprecision(4) << "(" << V_z[i] << "," << V_r[i] << ")" << " ";

				}
				t = t + dz[0];
			}
			else if (p == 1) {
				if ((V_z[i] <= (t + dz[p])) && (V_r[i] <= dr[0]) && (V_z[i] > t))
				{
					V_z[q] = 8 * V_z[i];
					V_r[q] = 8 * V_r[i];
					//printf(" ( %4.2f , %4.2f )", V_z[i], V_r[i]);
					//cout << " " << setw(4) << setfill(' ') << setprecision(4) << "(" << V_z[i] << "," << V_r[i] << ")" << " ";
					q++;
				}
				t = t + dz[p];
				s = s + dr[p - 1] + delta;
			}

			else if (p == n) {
				if ((V_z[i] <= (t + dz[p] + 1)) && (V_r[i] <= s) && (V_z[i] > t))
				{
					V_z[q] = 8 * V_z[i];
					V_r[q] = 8 * V_r[i];
					//printf(" ( %4.2f , %4.2f )", V_z[i], V_r[i]);
					//cout << " " << setw(4) << setfill(' ') << setprecision(4) << "(" << V_z[i] << "," << V_r[i] << ")" << " ";
					q++;
				}
			}

			else {
				if ((V_z[i] <= ((t + dz[p]))) && (V_r[i] <= ((s + dr[p - 1]))) && (V_z[i] > (t)))
				{
					V_z[q] = 8 * V_z[i];
					V_r[q] = 8 * V_r[i];
					//printf(" ( %4.2f , %4.2f )", V_z[i], V_r[i]);
					//cout << " " << setw(4) << setfill(' ') << setprecision(4) << "(" << V_z[i] << "," << V_r[i] << ")" << " ";
					q++;

				}
				t = t + dz[p];
				s = s + dr[p - 1] + delta;
			}
		}
	}
	cout << endl;
	k = q;
	return k;//������ɨ�赽���������
}
void paint_2(int k, Grid_Array** grid, int* _N, int* _M, double* dz, double* dr) //����//
{
	int i, j;
	double temp, x, y, xy, minlengthx_1, minlengthy_1, minlengthx_2, minlengthy_2, minlength;
	x = dz[0] / _N[0];
	for (i = 0; i < N - 1; i++)//�ҳ�����z���
	{
		if (x < grid[0][i].h2)
			x = grid[0][i].h2;
	}
	y = dr[0] / _M[0];
	for (i = 0; i < M - 1; i++)//�ҳ�����R���
	{
		if (y < grid[i][0].h4)
			y = grid[i][0].h4;
	}
	xy = sqrt(pow(x, 2) + pow(y, 2));//����֮����ܵ�������
	for (i = 1; i < k; i++)//��ɨ��ĵ�����������z������������
	{
		for (j = 1; j <= k - i; j++)
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
	for (i = 0; i < k; i++)  //�ӵ�һ���㿪ʼ��ÿ�����ҳ�����������������������֮��֮�ڣ�z��r���������z��r������δ���ӹ��ĵ�������ӣ�����ԭ�����Ǻ������
	{
		for (j = i + 1; j < k; j++)
		{
			if (sqrt(pow(V_r[j] - V_r[i], 2) + pow(V_z[j] - V_z[i], 2)) <= 10 * xy &&
				fabs(V_r[j] - V_r[i]) <= 10 * y && fabs(V_z[j] - V_z[i]) <= 10 * x)
			{
				line(int(V_z[i]), int(V_r[i]), int(V_z[j]), int(V_r[j]));
				break;
			}
		}
	}
	for (i = 1; i < k; i++)//��ɨ��ĵ�����������z���꽵������
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
	for (i = 0; i < k; i++)//�ӵ�һ���㿪ʼ��ÿ�����ҳ�����������������������֮��֮�ڣ�z��r���������z��r������δ���ӹ��ĵ�������ӣ�����ԭ�����Ǻ������
	{
		for (j = i + 1; j < k; j++)
		{
			if (sqrt(pow(V_r[j] - V_r[i], 2) + pow(V_z[j] - V_z[i], 2)) <= 10 * xy &&
				fabs(V_r[j] - V_r[i]) <= 10 * y && fabs(V_z[j] - V_z[i]) <= 10 * x)
			{
				line(int(V_z[i]), int(V_r[i]), int(V_z[j]), int(V_r[j]));
				break;
			}
		}
	}
}