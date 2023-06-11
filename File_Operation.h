#define FILE_OPERATION_H  //�����ļ��Ķ�д
#include "Grid_Array.h"
#include <stdio.h>
#include <stdlib.h>




int readdata1(FILE* file)
{
    tmp = 0;
    //int count = 0;
    int i;
    // ��ȡ����
    //fopen_s(&file, filepath, "rt");

    fscanf_s(file, "delta = %lf mm; n = %d;\n", &delta, &n); // �缫���

    dz = (double*)malloc(n * sizeof(double)); // ���ڵ缫�����
    fscanf_s(file, "Zi =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(file, " %lf", &dz[i]);
    }
    fscanf_s(file, "mm;\n");

    _N = (int*)malloc(n * sizeof(int)); // ���ڵ缫�䲽��
    fscanf_s(file, "Ni =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(file, " %d", &_N[i]);
    }
    fscanf_s(file, ";\n");

    V = (double*)malloc(n * sizeof(double)); // �缫��λ
    fscanf_s(file, "Vi =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(file, " %lf", &V[i]);
    }
    fscanf_s(file, "V;\n");

    VI = (double*)malloc(n * sizeof(double)); // ������缫��λ
    fscanf_s(file, "Vi =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(file, " %lf", &VI[i]);
    }
    fscanf_s(file, "V;\n");

    fscanf_s(file, "r1 = %lf mm; M1 = %d; r2 = %lf mm; M2 = %d;\n", &r1, &M1, &r2, &M2);

    fscanf_s(file, "epsilon =%lf V;NST = %d;INS = %d.\n", &epsilon, &NST, &INS);

    if (fscanf_s(file, "%d\n", &tmp) == 1 && tmp == 3)
    {
        V1 = new int();
        fscanf_s(file, "dengweixian:%d V\n", &V1[0]);
    } 
    else
    {
        V1 = (int*)malloc((n + 10) * sizeof(int)); // ��λ��
        fscanf_s(file, "dengweixian:");
        for (i = 0; i < n + 10; i++)
        {
            fscanf_s(file, "%d", &V1[i]);
            count1 = count1++;
            if (feof(file))
            {
                break;
            }
        }
    }

    fclose(file);
    return 0;
}

int readdata2(FILE* file)
{

#if 0
    double a;        // �缫���
    int n;           // �缫����
    double* z;       // �缫���
    int* N;          // ���ڵ缫�仮�ֲ���
    double* V;       // �缫��λ
    double* VI;      //������缫��λ
    double* r;       // �缫�ڿװ뾶
    //double r2;       // �ӵ缫�ڿ׾����ص���ձ߽紦�ľ������
    int* M;      // r1��Χ�ڵȲ������ֵ�������,�ӵ缫�ڿ׾����ص���ձ߽紦�ľ������
    int NST;         // �����ӡ�ռ��λʱ���������
    int INS;         // ���ϵ�λ���Ⱦ��ֵʱ������
    double e;        // ��������
    int m_V;         // ɨ���λ����
    //double* V_scan;  // �ƶ�ɨ���λʱ�ݴ�ɨ���λ
    //int I,m = 0;    // mΪ���������IΪ�������
    double V1;       // Ҫ��ɨ��ȵ�λ�ߵĵ�λ���

    int i;
    m_V = 0;
#endif
    tmp = 0;
    //int count = 0;
    int i;
    //fopen_s(&file, filepath, "rt");

    fscanf_s(file, "delta = %lf mm; n = %d;\n", &delta, &n); // �缫���

    dz = (double*)malloc((n+1) * sizeof(double)); // ���ڵ缫�����
    fscanf_s(file, "Zi =");
    for (i = 0; i < n+1; i++)
    {
        fscanf_s(file, " %lf", &dz[i]);
    }
    fscanf_s(file, "mm;\n");

    _N = (int*)malloc((n + 1) * sizeof(int)); // ���ڵ缫�䲽��
    fscanf_s(file, "Ni =");
    for (i = 0; i < n + 1; i++)
    {
        fscanf_s(file, " %d", &_N[i]);
    }
    fscanf_s(file, ";\n");

    V = (double*)malloc(n * sizeof(double)); // �缫��λ
    fscanf_s(file, "Vi =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(file, " %lf", &V[i]);
    }
    fscanf_s(file, "V;\n");

    VI = (double*)malloc(n * sizeof(double)); // ������缫��λ
    fscanf_s(file, "Vi =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(file, " %lf", &VI[i]);
    }
    fscanf_s(file, "V;\n");


    dr = (double*)malloc((n) * sizeof(double)); // dr
    fscanf_s(file, "r =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(file, " %lf", &dr[i]);
    }
    fscanf_s(file, "mm;\n");


    _M = (int*)malloc((n) * sizeof(int)); // dr
    fscanf_s(file, "M =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(file, " %d", &_M[i]);
    }
    fscanf_s(file, ";\n");

    fscanf_s(file, "epsilon =%lf V;NST = %d;INS = %d.\n", &epsilon, &NST, &INS);

    if (fscanf_s(file, "%d\n", &tmp) == 1 && tmp == 3)
    {
        V1 = new int();
        fscanf_s(file, "dengweixian:%d V\n", &V1[0]);
    }
    else
    {
        V1 = (int*)malloc((n + 10) * sizeof(int)); // ��λ��
        fscanf_s(file, "dengweixian:");
        for (i = 0; i < n + 10; i++)
        {
            fscanf_s(file, "%d", &V1[i]);
            count1 = count1++;
            if (feof(file))
            {
                break;
            }
        }
    }

    fclose(file);
    return 0;
}



void writedata(FILE* file, Grid_Array** grid1, Grid_Array** grid2, Iteration_Process* head1,Iteration_Process* head2)
{
    //���������Ϣ
    int i;
    fprintf(file, "�������: %d����� \n", mode);
    fprintf(file, "�缫�ߴ�� = %.1lf mm\n", delta);
    fprintf(file, "�缫����n = %d\n", n);
    fprintf(file, "�缫���Zi = ");
    if (mode == 1)
    {
        for (i = 0; i < n; i++)
        {
            fprintf(file, "%.1lf ", dz[i]);
        }
    }
    else
    {
        for (i = 0; i < n + 1; i++)
        {
            fprintf(file, "%.1lf ", dz[i]);
        }
    }
    fprintf(file, "mm\n");
    fprintf(file, "���򻮷�Ҫ��Ni = ");
    if (mode == 1)
    {
        for (i = 0; i < n; i++)
        {
            fprintf(file, "%d ", _N[i]);
        }
    }
    else
    {
        for (i = 0; i < n + 1; i++)
        {
            fprintf(file, "%d ", _N[i]);
        }
    }
    
    fprintf(file, "\n");

    fprintf(file, "�缫��ѹVi = ");
    for (i = 0; i < n; i++)
    {
        fprintf(file, "%.1lf ", V[i]);
    }
    fprintf(file, "V\n");

    fprintf(file,"�������ѹVi = ");
    for (i = 0; i < n; i++)
    {
        fprintf(file, "%.1lf ", VI[i]);
    }
    fprintf(file, "V\n");

    if (mode == 2)
    {
        fprintf(file, "��ֱ���r = ");
        for (i = 0; i < n - 1; i++)
        {
            fprintf(file, " %.1lf", dr[i]);
        }
        fprintf(file, "mm;\n");
    }
    else
    {
        fprintf(file, "��ֱ���r1 = %lf mm, r2 = %lf mm\n", r1, r2);
    }
    

    if (mode == 2)
    {
        fprintf(file, "���򻮷�Ҫ��M = ");
        for (i = 0; i < n - 1; i++)
        {
            fprintf(file, " %d", _M[i]);
        }
        fprintf(file, ";\n");
    }
    else
    {
        fprintf(file, "���򻮷�Ҫ��M1 = %d, M2 = %d\n", M1, M2);
    }
    


    fprintf(file, "�������Ȧ� = %.7lf V\n", epsilon);

    fprintf(file, "�����λ�����NST = %d\n", NST);

    fprintf(file, "���ϵ�λ��ֵ������INS = %d\n", INS);

    if (tmp == 3)
    {
        fprintf(file, "�ȿ�ɨ���λ�߼����%d V\n", V1[0]);
    }
    else
    {
        fprintf(file, "����ɨ���λ�߼��:");
        for (i = 0; i < count1; i++)
        {
            fprintf(file, " %d", V1[i]);
        }
    }
    fprintf(file, "\n");
    fprintf(file, "��������M = %d, N = %d\n\n", M, N);

    //������������
    //������������z
    double* z = new double[N];
    z[0] = 0;
    for (i = 1; i < N; i++)
    {
        if (z[1] == 0)
        {
            z[1] = grid1[M - 1][2].h1;
        }
        z[i] = z[i - 1] + grid1[M-1][i - 1].h2;
    }

    fprintf(file, "��������꣺\n         �������� ");
    for (i = 0; i < N; i++)
    {
        fprintf(file, "%6.2lf  ", z[i]);
    }
    fprintf(file, "\n��������\n");
    for (i = 0; i < M; i++)
    {
        for (int j = 0; j < N + 1; j++)
        {
            {
                if (j == 0)
                {
                    if (grid1[0][M - 1].r == 0)
                    {
                        if (grid1[1][M - 1].r == 0)
                        {
                            grid1[1][M - 1].r = grid1[2][M - 1].r + grid1[2][M - 1].h1;
                        }
                        grid1[0][M - 1].r = grid1[1][M - 1].r + delta;
                    }
                    fprintf(file, "%6.2lf            ", grid1[i][M-1].r);
                }
                else if (j != N)
                {
                    fprintf(file, "   +    ");
                }
                else
                {
                    fprintf(file, "   +   \n");
                }
            }
        }
    }



    //��ӡ�ް�������
    fprintf(file, "\n\n�ް�������\n");
    //��ӡ��������
    fprintf(file, "\n���������Ϣ��\n");
    fprintf(file, "��������: %d ��\n", iteration_times_1);
    fprintf(file, "ÿ�ε����ļ������Ӧأ��������ȣ����в��ƽ���в�:\n");
    fprintf(file, "�����ִ�          �������Ӧ�          ��������          ƽ���в�:\n");
    Iteration_Process* p;
    p = head1;
    p = p->next;
    for (i = 0; i < iteration_times_1; i++)
    {
        fprintf(file, "%2d��%2d��             %.7lf         %.7lf         %.7lf\n", p->round,p->times, p->omega_r, p->max_res, p->avg_res);
        if (p->next != NULL)
        {
            p = p->next;
        }
    }

    //��ӡ���ϵ�ѹ
    if (mode == 1)
    {
        double z_1, v_1, z_0, v_0;
        double z_pace = z0 / INS;

        z_0 = 0;
        v_0 = 0;
        fprintf(file, "\n�ް���������ϵ��������λ��\n");
        fprintf(file, "����r             ��λV\n");
        fprintf(file, " %6.2lf     %11.7f \n", z_0, v_0);

        for (i = 0; i < INS - 1; i++) {
            z_1 = z_0 + z_pace;
            for (int j = 0; j < (N - 1); j++) {
                if ((zi[j] <= z_1) && (z_1 < zi[j + 1])) {
                    v_1 = (z_1 - zi[j]) * (z_1 - zi[j + 1]) / (zi[j - 1] - zi[j]) / (zi[j - 1] - zi[j + 1]) * grid1[M - 1][j - 1].voltage\
                        + (z_1 - zi[j - 1]) * (z_1 - zi[j + 1]) / (zi[j] - zi[j - 1]) / (zi[j] - zi[j + 1]) * grid1[M - 1][j].voltage\
                        + (z_1 - zi[j - 1]) * (z_1 - zi[j]) / (zi[j + 1] - zi[j - 1]) / (zi[j + 1] - zi[j]) * grid1[M - 1][j + 1].voltage;
                    fprintf(file, " %6.2lf     %11.7f \n", z_1, v_1);
                    z_0 = z_1;
                    v_0 = v_1;
                }
            }
        }
        fprintf(file, " %6.2lf     %11.7f \n", z0, grid1[0][N - 1].voltage);
    }
    else
    {
        double z_1, v_1;
        double z_0 = 0;
        for (i = 0; i < n + 1; i++) {
            z_0 = z_0 + dz[i];
        }
        double z_pace = z_0 / INS;
        double z_2 = 0;
        double v_2 = 0;

        fprintf(file, "\n�ް���������ϵ��������λ��\n");
        fprintf(file, "����r             ��λV\n");
        fprintf(file, " %6.2lf     %11.7f \n", z_pace, v_2);

        for (i = 0; i < INS - 1; i++) {
            z_1 = z_2 + z_pace;
            for (int j = 0; j < (N - 1); j++) {
                if ((zi[j] <= z_1) && (z_1 < zi[j + 1])) {
                    v_1 = (z_1 - zi[j]) * (z_1 - zi[j + 1]) / (zi[j - 1] - zi[j]) / (zi[j - 1] - zi[j + 1]) * grid1[M - 1][j - 1].voltage\
                        + (z_1 - zi[j - 1]) * (z_1 - zi[j + 1]) / (zi[j] - zi[j - 1]) / (zi[j] - zi[j + 1]) * grid1[M - 1][j].voltage\
                        + (z_1 - zi[j - 1]) * (z_1 - zi[j]) / (zi[j + 1] - zi[j - 1]) / (zi[j + 1] - zi[j]) * grid1[M - 1][j + 1].voltage;
                    fprintf(file, " %6.2lf     %11.7f \n", z_1, v_1);
                    z_2 = z_1;
                    v_2 = v_1;
                }
            }
        }
        fprintf(file, " %6.2lf     %11.7f \n", z_0, grid1[0][N - 1].voltage);
        //cout << z_0 << "��ѹ" << 100 << endl;
    }
    

    //��ӡ�����λ
    fprintf(file, "\n�ް��������λ��\n");
    for (int i = 0; i < M; i += NST)
    {
        for (int j = 0; j < N; j += NST)
        {
            fprintf(file, "%11.7f ", grid1[i][j].voltage); 
        }
        fprintf(file,"\n");
    }
       



    //��ӡ�а�������
    fprintf(file, "\n\n�а�������\n");
    //��ӡ��������
    fprintf(file, "\n���������Ϣ��\n");
    fprintf(file, "��������: %d ��\n", iteration_times_1);
    fprintf(file, "ÿ�ε����ļ������Ӧأ��������ȣ����в��ƽ���в�:\n");
    fprintf(file, "��������          �������Ӧ�          ��������          ƽ���в�:\n");
    Iteration_Process* q;
    q = head2;
    q = q->next;
    for (i = 0; i < iteration_times_2; i++)
    {
        fprintf(file, "%2d��%2d��             %.7lf         %.7lf         %.7lf\n", q->round, q->times, q->omega_r, q->max_res, q->avg_res);
        if (q->next != NULL)
        {
             q = q->next;
        }

    }

    //��ӡ���ϵ�ѹ
    if (mode == 1)
    {
        double z_1, v_1, z_0, v_0;
        double z_pace = z0 / INS;

        z_0 = 0;
        v_0 = 0;
        fprintf(file, "\n�ް���������ϵ��������λ��\n");
        fprintf(file, "����r             ��λV\n");
        fprintf(file, " %6.2lf     %11.7f \n", z_0, v_0);

        for (i = 0; i < INS - 1; i++) {
            z_1 = z_0 + z_pace;
            for (int j = 0; j < (N - 1); j++) {
                if ((zi[j] <= z_1) && (z_1 < zi[j + 1])) {
                    v_1 = (z_1 - zi[j]) * (z_1 - zi[j + 1]) / (zi[j - 1] - zi[j]) / (zi[j - 1] - zi[j + 1]) * grid2[M - 1][j - 1].voltage\
                        + (z_1 - zi[j - 1]) * (z_1 - zi[j + 1]) / (zi[j] - zi[j - 1]) / (zi[j] - zi[j + 1]) * grid2[M - 1][j].voltage\
                        + (z_1 - zi[j - 1]) * (z_1 - zi[j]) / (zi[j + 1] - zi[j - 1]) / (zi[j + 1] - zi[j]) * grid2[M - 1][j + 1].voltage;
                    fprintf(file, " %6.2lf     %11.7f \n", z_1, v_1);
                    z_0 = z_1;
                    v_0 = v_1;
                }
            }
        }
        fprintf(file, " %6.2lf     %11.7f \n", z0, grid2[0][N - 1].voltage);
    }
    else
    {
        double z_1, v_1;
        double z_0 = 0;
        for (i = 0; i < n + 1; i++) {
            z_0 = z_0 + dz[i];
        }
        double z_pace = z_0 / INS;
        double z_2 = 0;
        double v_2 = 0;

        fprintf(file, "\n�ް���������ϵ��������λ��\n");
        fprintf(file, "����r             ��λV\n");
        fprintf(file, " %6.2lf     %11.7f \n", z_pace, v_2);

        for (i = 0; i < INS - 1; i++) {
            z_1 = z_2 + z_pace;
            for (int j = 0; j < (N - 1); j++) {
                if ((zi[j] <= z_1) && (z_1 < zi[j + 1])) {
                    v_1 = (z_1 - zi[j]) * (z_1 - zi[j + 1]) / (zi[j - 1] - zi[j]) / (zi[j - 1] - zi[j + 1]) * grid2[M - 1][j - 1].voltage\
                        + (z_1 - zi[j - 1]) * (z_1 - zi[j + 1]) / (zi[j] - zi[j - 1]) / (zi[j] - zi[j + 1]) * grid2[M - 1][j].voltage\
                        + (z_1 - zi[j - 1]) * (z_1 - zi[j]) / (zi[j + 1] - zi[j - 1]) / (zi[j + 1] - zi[j]) * grid2[M - 1][j + 1].voltage;
                    fprintf(file, " %6.2lf     %11.7f \n", z_1, v_1);
                    z_2 = z_1;
                    v_2 = v_1;
                }
            }
        }
    }
    //��ӡ�����λ
    fprintf(file, "\n�а��������λ��\n");
    for (int i = 0; i < M; i += NST)
    {
        for (int j = 0; j < N; j += NST)
        {
            fprintf(file, "%11.7f ", grid2[i][j].voltage);
        }
        fprintf(file, "\n");
    }
}
