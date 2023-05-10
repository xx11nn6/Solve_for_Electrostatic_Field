#include <iostream>
#include <malloc.h>
using namespace std;

struct Grid_Array  //����糡����ṹ
{
	float voltage;  //�洢��������ѹ
	bool is_margin;  //�ж��Ƿ��Ǳ߽磬1��ʾ�Ǳ߽�
};

int M, N;//������M*N

//��������
void grid_initialize(Grid_Array** grid, float cathode_voltage, float screen_voltage);
void node_set(Grid_Array** grid, int m, int n, float voltage, bool is_margin);


int main()
{
	//��ʼ��
	float V_c, V_s; //�����������ӫ������ѹ
	M = 20; N = 20;
	V_c = 0; V_s = 100; //�����������ӫ������ѹ�ֱ�Ϊ0��100V
	int mode = 1;  //���幤��ģʽ��1Ϊ��һ����ܣ�2Ϊ�ڶ������
	float z0, r0;  //���徶�������Ŀ��
	

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
	grid_initialize(grid,V_c,V_s);
	node_set(grid,3,1,50,true);


	//c++����cout�������ʹ��cout�������
	cout << "grid:" << endl;
	for (int i = 0; i < M; i++) 
	{
		for (int j = 0; j < N; j++) 
		{
			cout << "(" << grid[i][j].voltage << ", " << grid[i][j].is_margin << ") ";
		}
		cout << endl;
	}

	free(grid);//���н������ͷ��ڴ�

	return 0;
}




//���ڳ�ʼ���糡����
void grid_initialize(Grid_Array** grid,float cathode_voltage,float screen_voltage)
					//����糡���������ȡ��߶ȡ�������ѹ��ӫ������ѹ
{
	//�Ƚ���λ���㣬ȫ������Ϊ�Ǳ߽�
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			grid[i][j].voltage = 0;
			grid[i][j].is_margin = false;
		}
	}
	//���ù�������λ�ͱ߽����
	for (int i = 0; i < M; i++)
	{
		grid[i][0].voltage = cathode_voltage;
		grid[i][0].is_margin = true;
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
