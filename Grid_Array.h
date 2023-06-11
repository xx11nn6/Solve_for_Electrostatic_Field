#ifndef GRID_ARRAY_H  //���������ͷź�����һЩȫ�ֱ���
#define GRID_ARRAY_H

//����ȫ�ֱ���M
extern int mode;			 //����ģʽ1Ϊһ�࣬2Ϊ����
extern int M;			  	 //������MN
extern int N;
extern double delta;		 // �缫���
extern int n;				 // �缫��������ӫ������������������
extern double* dz;           // ÿ���缫�������࣬������ܾ��õ��˱���
extern int* _N;              // ���ڵ缫�仮�ֲ���(ͨ��)
extern double* V;            // �����缫��ѹ
extern double* VI;           // ������缫��λ(ͨ��)
extern int M1, M2;           // ��ֱ����������ֵ�Ҫ�󣬵�һ�����ʹ�ô˱���
extern double r1,r2;         // r2Ϊ�缫����糡����ȣ�r1Ϊ�缫�׶˵���ľ���(��һ�����ʹ��)
extern int* _M;              // �缫֮�侶������Ҫ���ֵ����������ڶ������ʹ�ô˱���
extern double* dr;           // �缫�ڿװ뾶(�ڶ���)
extern double epsilon;       // ��������(ͨ��)
extern double omega;		 // �������Ӧ�
extern int NST;              // �����ӡ�ռ��λʱ���������(ͨ��)
extern int INS;              // ���ϵ�λ���Ⱦ��ֵʱ������(ͨ��)
extern int* V1;              // Ҫ��ɨ��ȵ�λ�ߵĵ�λ������ߵ�λֵ(ͨ��)

//�ļ�����
extern FILE* handle;
extern errno_t err;

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

//ʹ��delete���ͷ��ڴ�
void free_grid(Grid_Array** grid)
{
	for (int i = 0; i < M; i++) {
		delete[] grid[i];
	}
	delete[] grid;
}

#endif
