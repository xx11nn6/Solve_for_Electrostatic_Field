#define GRID_INITIALIZE_H  //���ڳ�ʼ���������
#include "Grid_Array.h"


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
	for (int i = 1; i < M; i++)
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
	int* Ms = new int[n - 1];
	int* Ns = new int[n + 1];
	double* rs = new double[n - 1];

	Ms[0] = _M[n - 2] + 1;  //������_M[0]�ģ�����û�뵽�ɹ�����Ҫ�����񻮷ֺ�����ϵ�Ƿ��ŵ�......
	Ns[0] = _N[0];
	rs[0] = dr[0] + delta;
	for (int i = 1; i < (n + 1); i++)
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
			grid[i][j].r = 0;
			grid[i][j].h1 = 0;
			grid[i][j].h2 = 0;
			grid[i][j].h3 = 0;
			grid[i][j].h4 = 0;
		}
	}

	//�ٴα��������г�ʼ���ͱ߽�����
	for (int k = 0; k < (n - 1); k++)
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
									grid[i][j].h1 = dz[n - k + l - 1] / _N[n - k + l - 1];
									grid[i][j].h2 = dz[n - k + l] / _N[n - k + l];
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
								r = rs[n - k - 2] - (double(i - Ms[k - 1] - 1) / _M[n - k - 2] * dr[n - k - 2]) - delta;
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
							grid[i][j].h2 = dz[n - k + l] / _N[n - k + l];
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
			grid[i][N - 1].voltage = screen_voltage;
			grid[i][N - 1].is_margin = true;
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
