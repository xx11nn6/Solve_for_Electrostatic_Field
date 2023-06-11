//ͷ�ļ�����
#include<stdio.h>
#include<stdlib.h>

//��������
int readdata1(const char* filepath);//һ����ܶ��ļ�
int readdata2(const char* filepath);//������ܶ��ļ�

//�ļ�����
FILE* handle;
errno_t err;

//ȫ�ֱ���
int type = 0;          //ȷ���������
int tmp = 0;           //ȷ����λ��ɨ�跽ʽ
double delta;         //�缫���(ͨ��)
int n;                //�缫����(ͨ��)
double* Zi;           //�缫���(ͨ��)
int* Ni;              //���ڵ缫�仮�ֲ���(ͨ��)
double* Vi;           // �缫��λ(ͨ��)
double* VI;           //������缫��λ(ͨ��)
double r1;            // �缫�ڿװ뾶(��һ��)
int M1;               // r1��Χ�ڵȲ������ֵ�������(��һ��)
double r2;            // �ӵ缫�ڿ׾����ص���ձ߽紦�ľ������(��һ��)
int M2;               //�ӵ缫�ڿ׾����ص���ձ߽紦�ľ������(��һ��)
int* M;               // r��Χ�ڵȲ������ֵ�������(�ڶ���)
double* r;            // �缫�ڿװ뾶(�ڶ���)
double epsilon;       // ��������(ͨ��)
int NST;              // �����ӡ�ռ��λʱ���������(ͨ��)
int INS;              // ���ϵ�λ���Ⱦ��ֵʱ������(ͨ��)
int* V1;              // Ҫ��ɨ��ȵ�λ�ߵĵ�λ������ߵ�λֵ(ͨ��)
int i;                //ѭ���ò���
int count = 0;        //һ�����������ж�ɨ���λ�ߵ�����

int readdata1(const char* filepath)
{


#if 0
    double delta;//�缫���
    int n;//�缫����
    double* Zi;//�缫���
    int* Ni;//���ڵ缫�仮�ֲ���
    double* Vi;// �缫��λ
    double* VI;//������缫��λ
    double r1;// �缫�ڿװ뾶
    int M1;// r1��Χ�ڵȲ������ֵ�������
    double r2;// �ӵ缫�ڿ׾����ص���ձ߽紦�ľ������
    int M2;//�ӵ缫�ڿ׾����ص���ձ߽紦�ľ������
    double epsilon;// ��������
    int NST;// �����ӡ�ռ��λʱ���������
    int INS;// ���ϵ�λ���Ⱦ��ֵʱ������
    int* V1;// Ҫ��ɨ��ȵ�λ�ߵĵ�λ���

    int i;
    int m_V = 0;// ɨ���λ����
#endif
    
    // ��ȡ����
    //fopen_s(&handle, filepath, "rt");
    
    fscanf_s(handle, "delta = %lf mm; n = %d;\n", &delta, &n); // �缫���

    Zi = (double*)malloc((n + 1) * sizeof(double)); // ���ڵ缫�����
    fscanf_s(handle, "Zi =");
    for (i = 0; i < n + 1; i++)
    {
        fscanf_s(handle, " %lf", &Zi[i]);
    }
    fscanf_s(handle, "mm;\n");

    Ni = (int*)malloc(n * sizeof(int)); // ���ڵ缫�䲽��
    fscanf_s(handle, "Ni =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(handle, " %d", &Ni[i]);
    }
    fscanf_s(handle, ";\n");

    Vi = (double*)malloc(n * sizeof(double)); // �缫��λ
    fscanf_s(handle, "Vi =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(handle, " %lf", &Vi[i]);
    }
    fscanf_s(handle, "V;\n");

    VI = (double*)malloc(n * sizeof(double)); // ������缫��λ
    fscanf_s(handle, "Vi =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(handle, " %lf", &VI[i]);
    }
    fscanf_s(handle, "V;\n");

    fscanf_s(handle, "r1 = %lf mm; M1 = %d; r2 = %lf mm; M2 = %d;\n", &r1, &M1, &r2,&M2);

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

    printf("n = %d\n",n);

    printf("Zi = ");
    for (i = 0; i < n; i++)
    {
        printf("%.1lf ", Zi[i]);
    }
    printf("mm\n");

    printf("Ni = ");
    for (i = 0; i < n; i++)
    {
        printf("%d ", Ni[i]);
    }
    printf("\n");

    printf("Vi = ");
    for (i = 0; i < n; i++)
    {
        printf("%.1lfV ", Vi[i]);
    }
    printf("\n");

    printf("Vi = ");
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

    // �ͷŶ�̬������ڴ�
    free(Zi);
    free(Ni);
    free(Vi);
    free(VI);
    free(V1);

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
    
    //fopen_s(&handle, filepath, "rt");

    fscanf_s(handle, "delta = %lf mm; n = %d;\n", &delta, &n); // �缫���

    Zi = (double*)malloc((n + 1) * sizeof(double)); // ���ڵ缫�����
    fscanf_s(handle, "Zi =");
    for (i = 0; i < n + 1; i++)
    {
        fscanf_s(handle, " %lf", &Zi[i]);
    }
    fscanf_s(handle, "mm;\n");

    Ni = (int*)malloc((n + 1) * sizeof(int)); // ���ڵ缫�䲽��
    fscanf_s(handle, "Ni =");
    for (i = 0; i < n + 1; i++)
    {
        fscanf_s(handle, " %d", &Ni[i]);
    }
    fscanf_s(handle, ";\n");

    Vi = (double*)malloc(n * sizeof(double)); // �缫��λ
    fscanf_s(handle, "Vi =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(handle, " %lf", &Vi[i]);
    }
    fscanf_s(handle, "V;\n");

    VI = (double*)malloc(n * sizeof(double)); // ������缫��λ
    fscanf_s(handle, "Vi =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(handle, " %lf", &VI[i]);
    }
    fscanf_s(handle, "V;\n");


    r = (double*)malloc((n) * sizeof(double)); // r
    fscanf_s(handle, "r =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(handle, " %lf", &r[i]);
    }
    fscanf_s(handle, "mm;\n");


    M = (int*)malloc((n) * sizeof(int)); // r
    fscanf_s(handle, "M =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(handle, " %d", &M[i]);
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
        printf("%.1lf ", Zi[i]);
    }
    printf("mm\n");

    printf("Ni = ");
    for (i = 0; i < n + 1; i++)
    {
        printf("%d ", Ni[i]);
    }
    printf("\n");

    printf("Vi = ");
    for (i = 0; i < n; i++)
    {
        printf("%.1lf ", Vi[i]);
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
        printf(" %.1lf", r[i]);
    }
    printf("mm;\n");

    printf("M = ");
    for (i = 0; i < n - 1; i++)
    {
        printf(" %d", M[i]);
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
    
    free(Zi);
    free(Ni);
    free(Vi);
    free(VI);
    free(r);
    free(M);
    free(V1);
    

    return 0;
}


    
int main()
{
    if ((err = fopen_s(&handle, "D:\\CADTest\\2023\\2\\����.txt", "rt")) != 0)
    {
        printf("�Ҳ����ļ��������룺%d\n", err);
        return -1;
    }
    
    
    fopen_s(&handle, "D:\\CADTest\\2023\\2\\����.txt", "rt");
    if (fscanf_s(handle, "%d\n", &type) == 1 && type == 2)
    {
        readdata2("D:\\CADTest\\2023\\2\\����.txt");
        fclose(handle);
    }
    else
    {
        readdata1("D:\\CADTest\\2023\\2\\����.txt");
        fclose(handle);
    }
    return 0;
}
