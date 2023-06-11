#define FILE_OPERATION_H  //�����ļ��Ķ�д
#include "Grid_Array.h"
#include<stdio.h>
#include<stdlib.h>




int readdata1(const char* filepath)
{
    int tmp = 0;
    int count = 0;
    int i;
    // ��ȡ����
    //fopen_s(&handle, filepath, "rt");

    fscanf_s(handle, "delta = %lf mm; n = %d;\n", &delta, &n); // �缫���

    dz = (double*)malloc((n+1) * sizeof(double)); // ���ڵ缫�����
    fscanf_s(handle, "Zi =");
    for (i = 0; i < n+1; i++)
    {
        fscanf_s(handle, " %lf", &dz[i]);
    }
    fscanf_s(handle, "mm;\n");

    _N = (int*)malloc(n * sizeof(int)); // ���ڵ缫�䲽��
    fscanf_s(handle, "Ni =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(handle, " %d", &_N[i]);
    }
    fscanf_s(handle, ";\n");

    V = (double*)malloc(n * sizeof(double)); // �缫��λ
    fscanf_s(handle, "Vi =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(handle, " %lf", &V[i]);
    }
    fscanf_s(handle, "V;\n");

    VI = (double*)malloc(n * sizeof(double)); // ������缫��λ
    fscanf_s(handle, "Vi =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(handle, " %lf", &VI[i]);
    }
    fscanf_s(handle, "V;\n");

    fscanf_s(handle, "r1 = %lf mm; M1 = %d; r2 = %lf mm; M2 = %d;\n", &r1, &M1, &r2, &M2);

    fscanf_s(handle, "epsilon =%lf V;NST = %d;INS = %d.\n", &epsilon, &NST, &INS);

    if (fscanf_s(handle, "%d\n", &tmp) == 1 && tmp == 3)
    {
        fscanf_s(handle, "dengweixian:%d V\n", &V1);
    }
    else
    {
        V1 = (int*)malloc((n + 10) * sizeof(int)); // ��λ��
        fscanf_s(handle, "dengweixian:");
        for (i = 0; i < n + 10; i++)
        {
            fscanf_s(handle, "%d", &V1[i]);
            count = count++;
            if (feof(handle))
            {
                break;
            }
        }
    }

    fclose(handle);

    // ��ӡ��ȡ�������ݣ��ɸ�����Ҫ���е���
    printf("delta = %.1lf mm\n", delta);

    printf("n = %d\n", n);

    printf("Zi = ");
    for (i = 0; i < n; i++)
    {
        printf("%.1lf ", dz[i]);
    }
    printf("mm\n");

    printf("Ni = ");
    for (i = 0; i < n; i++)
    {
        printf("%d ", _N[i]);
    }
    printf("\n");

    printf("Vi = ");
    for (i = 0; i < n; i++)
    {
        printf("%.1lfV ", V[i]);
    }
    printf("\n");

    printf("VI = ");
    for (i = 0; i < n; i++)
    {
        printf("%.1lfV ", VI[i]);
    }
    printf("\n");

    printf("r1 = %.1lf mm\n", r1);

    printf("M1 = %.d\n", M1);

    printf("r2 = %.1lf mm\n", r2);

    printf("M2 = %.d\n", M2);

    printf("epsilon = %.7lf V\n", epsilon);

    printf("NST = %d\n", NST);

    printf("INS = %d\n", INS);

    if (tmp == 3)
    {
        printf("��λ�߼����%d V\n", V1);
    }
    else
    {
        printf("��λ�߼��:");
        for (i = 0; i < count; i++)
        {
            printf(" %d", V1[i]);
        }
    }

    return 0;
}

int readdata2(const char* filepath)
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
    int tmp = 0;
    int count = 0;
    int i;
    //fopen_s(&handle, filepath, "rt");

    fscanf_s(handle, "delta = %lf mm; n = %d;\n", &delta, &n); // �缫���

    dz = (double*)malloc((n+1) * sizeof(double)); // ���ڵ缫�����
    fscanf_s(handle, "Zi =");
    for (i = 0; i < n+1; i++)
    {
        fscanf_s(handle, " %lf", &dz[i]);
    }
    fscanf_s(handle, "mm;\n");

    _N = (int*)malloc((n + 1) * sizeof(int)); // ���ڵ缫�䲽��
    fscanf_s(handle, "Ni =");
    for (i = 0; i < n + 1; i++)
    {
        fscanf_s(handle, " %d", &_N[i]);
    }
    fscanf_s(handle, ";\n");

    V = (double*)malloc(n * sizeof(double)); // �缫��λ
    fscanf_s(handle, "Vi =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(handle, " %lf", &V[i]);
    }
    fscanf_s(handle, "V;\n");

    VI = (double*)malloc(n * sizeof(double)); // ������缫��λ
    fscanf_s(handle, "Vi =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(handle, " %lf", &VI[i]);
    }
    fscanf_s(handle, "V;\n");


    dr = (double*)malloc((n) * sizeof(double)); // dr
    fscanf_s(handle, "r =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(handle, " %lf", &dr[i]);
    }
    fscanf_s(handle, "mm;\n");


    _M = (int*)malloc((n) * sizeof(int)); // dr
    fscanf_s(handle, "M =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(handle, " %d", &_M[i]);
    }
    fscanf_s(handle, ";\n");

    fscanf_s(handle, "epsilon =%lf V;NST = %d;INS = %d.\n", &epsilon, &NST, &INS);

    if (fscanf_s(handle, "%d\n", &tmp) == 1 && tmp == 3)
    {
        fscanf_s(handle, "dengweixian:%d V\n", &V1);
    }
    else
    {
        V1 = (int*)malloc((n + 10) * sizeof(int)); // ��λ��
        fscanf_s(handle, "dengweixian:");
        for (i = 0; i < n + 10; i++)
        {
            fscanf_s(handle, "%d", &V1[i]);
            count = count++;
            if (feof(handle))
            {
                break;
            }
        }
    }

    fclose(handle);

    // ��ӡ��ȡ�������ݣ��ɸ�����Ҫ���е���
    printf("delta = %.1lf mm\n", delta);
    printf("n = %d\n", n);
    printf("Zi = ");
    for (i = 0; i < n + 1; i++)
    {
        printf("%.1lf ", dz[i]);
    }
    printf("mm\n");

    printf("Ni = ");
    for (i = 0; i < n + 1; i++)
    {
        printf("%d ", _N[i]);
    }
    printf("\n");

    printf("Vi = ");
    for (i = 0; i < n; i++)
    {
        printf("%.1lf ", V[i]);
    }
    printf("V\n");

    printf("Vi = ");
    for (i = 0; i < n; i++)
    {
        printf("%.1lf ", VI[i]);
    }
    printf("V\n");


    printf("r = ");
    for (i = 0; i < n - 1; i++)
    {
        printf(" %.1lf", dr[i]);
    }
    printf("mm;\n");

    printf("M = ");
    for (i = 0; i < n - 1; i++)
    {
        printf(" %d", _M[i]);
    }
    printf(";\n");


    printf("epsilon = %.7lf V\n", epsilon);

    printf("NST = %d\n", NST);

    printf("INS = %d\n", INS);

    if (tmp == 3)
    {
        printf("��λ�߼����%d V\n", V1);
    }
    else
    {
        printf("��λ�߼��:");
        for (i = 0; i < count; i++)
        {
            printf(" %d", V1[i]);
        }
    }



    return 0;
}