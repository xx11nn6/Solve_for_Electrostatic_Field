#include <iostream>
#include <malloc.h>
#include <numeric>
using namespace std;

struct Grid_Array  //����糡����ṹ
{
	float voltage;  //�洢��������ѹ
	bool is_margin;  //�ж�������Ƿ��Ǳ߽磬1��ʾ�Ǳ߽�
	bool on_axis;  //�ж��Ƿ������ϣ�
	float h1, h2, h3, h4;  //�洢����㵽��Χ�ľ���
	float r;  //�洢�����ľ������
};


//����ȫ�ֱ���
int M, N;//������M*N
int mode = 1;  //���幤��ģʽ��1Ϊ��һ����ܣ�2Ϊ�ڶ������
//��������
void grid_initialize(Grid_Array** grid, float cathode_voltage, float screen_voltage, float* dz, int* _N, float r1, float r2, int M1, int M2, float delta);
void node_set(Grid_Array** grid, int m, int n, float voltage, bool is_margin);


int main()
{
	//��ʼ��
	float V_c = 0, V_s = 100;  //�����������ӫ������ѹ�ֱ�Ϊ0��100V
	int n = 7;  //����缫��������ӫ������������������
	float delta = 0.5;//����缫��Ȧ�
	float z0 = 53.4, r0 = 32;  //���徶�������Ŀ��
	float r2 = 12, r1 = r0 - r2;  //r2Ϊ�缫����糡����ȣ�r1Ϊ�缫�׶˵���ľ���
	float dz[6];  //����ÿ���缫�ļ��
	int _N[6];  //����ÿ���缫֮����ȡ����������ˮƽ���򣩣����»���Ϊ����ȫ�ֱ���N����
	int M1, M2;  //��ֱ�������������ֵ�Ҫ��
	dz[0] = 5.2; dz[1] = 8.6; dz[2] = 8.6; dz[3] = 8.6; dz[4] = 8.6; dz[5] = 5.2;
	_N[0] = 3; _N[1] = 5; _N[2] = 5; _N[3] = 5; _N[4] = 5; _N[5] = 2;
	M1 = 11; M2 = 7;
	//������������M��N
	M = M1 + M2;
	N = _N[0] + _N[1] + _N[2] + _N[3] + _N[4] + _N[5] + n - 1;
	float epsilon = 0.0005;  //��������Ϊ0.005



	//c++�У�new���ڶ�̬�����ڴ�
	//���磺int* p = new int[10];
	//��仰��ʾ����10��int���͵Ŀռ䣬�����׵�ַ�洢��p��
	//��˿�����������һ���ɱ��С������
	//Ҫ�봴����ά���飬�����ʹ��˫��ָ��
	Grid_Array** grid = new Grid_Array* [M];  //��һ��ָ��ָ����
	for (int i = 0; i < M; i++) 
	{
		grid[i] = new Grid_Array[N];  //�ڶ���ָ����
	}
	grid_initialize(grid,V_c,V_s,dz,_N,r1,r2,M1,M2,delta);
	node_set(grid,3,1,50,true);


	//c++����cout�������ʹ��cout�������
	cout << "grid:" << endl;
	for (int i = 0; i < M; i++) 
	{
		for (int j = 0; j < N; j++) 
		{
			cout << "(" << grid[i][j].r << ", " << grid[i][j].h3 << "," << grid[i][j].h4 << ") ";
		}
		cout << endl;
	}


	free(grid);//���н������ͷ��ڴ�

	return 0;
}




//���ڳ�ʼ���糡����
void grid_initialize(Grid_Array** grid,float cathode_voltage,float screen_voltage,float* dz,int* _N,float r1,float r2,int M1,int M2,float delta)
					//����糡���������ȡ��߶ȡ�������ѹ��ӫ������ѹ�������ࣨ��������dz����ˮƽ�������񻮷�Ҫ�󣨴�������_N��
					//�缫�׵����r1���缫���r2����ֱ�������񻮷�Ҫ��M1,M2���缫���delta
{
	//�Ƚ���λ���㣬ȫ������Ϊ�Ǳ߽硢�����ϵ㣬�Լ���ֱ������
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			grid[i][j].voltage = 0;
			grid[i][j].is_margin = false;
			grid[i][j].on_axis = false;
			//���õ�һ��
			if (i == 0)
			{
				grid[i][j].h4 = 0;
				grid[i][j].h3 = r2 / (M2 - 1);
				grid[i][j].r = r1 + r2;
			}
			//����1~M2-1��
			else if (i < M2 - 1)
			{
				grid[i][j].h4 = r2 / (M2 - 1);
				grid[i][j].h3 = grid[i][j].h4;
				grid[i][j].r = r1 + r2 - ((r2 * i) / (M2 - 1));
			}
			//����M2��
			else if (i == M2 - 1)
			{
				grid[i][j].h4 = r2 / (M2 - 1);
				grid[i][j].h3 = r1 / (M1 - 1);
				grid[i][j].r = r1;
			}
			//����M2~M1-1��
			else if (i < (M2 + M1 - 1))
			{
				grid[i][j].h4 = r1 / (M1 - 1);
				grid[i][j].h3 = grid[i][j].h4;
				grid[i][j].r = r1 - ((r1 * (i - M2)) / (M1 - 1));
			}
			//�������ϵ�
			else
			{
				grid[i][j].h3 = 0;
				grid[i][j].h4 = r1 / (M1 - 1);
				grid[i][j].r = 0;
				grid[i][j].on_axis = true;
				grid[i][j].is_margin = true;
			}
		}
	}
	//���ù�����
	for (int i = 0; i < M; i++)
	{
		grid[i][0].voltage = cathode_voltage;
		grid[i][0].is_margin = true;
		grid[i][0].h1 = 0;  //������Ϊ��߽磬h1=0
		grid[i][0].h2 = dz[0] / (_N[0] - 1);  //�������ұߵ���һ���缫����Ϊdz[0],������Ϊ_N[0]��������h2
		
		
	}
	//����ӫ������λ�ͱ߽����
	for (int i = 0; i < M; i++)
	{
		grid[i][N-1].voltage = screen_voltage;
		grid[i][N-1].is_margin = true;
	}
}

//���ڸ�ָ�����������ڵ��趨��ѹֵ�ͱ߽�״̬������grid���飬���꣨m,n������ѹ���Ƿ��Ǳ߽�
void node_set(Grid_Array** grid,int m,int n,float voltage,bool is_margin)
{
	grid[m][n].voltage = voltage;
	grid[m][n].is_margin = is_margin;
}

//�߽��պ���������ʹ���糡
void boundary_close()
{

}

//c++��ʹ��delete���ͷ��ڴ�
void free(Grid_Array** grid)
{
	for (int i = 0; i < M; i++) {
		delete[] grid[i];
	}
	delete[] grid;
}
