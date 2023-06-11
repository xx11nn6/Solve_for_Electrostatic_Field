#define FILE_OPERATION_H  //包括文件的读写
#include "Grid_Array.h"
#include <stdio.h>
#include <stdlib.h>




int readdata1(FILE* file)
{
    tmp = 0;
    //int count = 0;
    int i;
    // 读取数据
    //fopen_s(&file, filepath, "rt");

    fscanf_s(file, "delta = %lf mm; n = %d;\n", &delta, &n); // 电极宽度

    dz = (double*)malloc(n * sizeof(double)); // 相邻电极间距离
    fscanf_s(file, "Zi =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(file, " %lf", &dz[i]);
    }
    fscanf_s(file, "mm;\n");

    _N = (int*)malloc(n * sizeof(int)); // 相邻电极间步长
    fscanf_s(file, "Ni =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(file, " %d", &_N[i]);
    }
    fscanf_s(file, ";\n");

    V = (double*)malloc(n * sizeof(double)); // 电极电位
    fscanf_s(file, "Vi =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(file, " %lf", &V[i]);
    }
    fscanf_s(file, "V;\n");

    VI = (double*)malloc(n * sizeof(double)); // 含鞍点电极电位
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
        V1 = (int*)malloc((n + 10) * sizeof(int)); // 等位线
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

    // 打印读取到的数据，可根据需要自行调整
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
        printf("等位线间隔：%d V\n", V1[0]);
    }
    else
    {
        printf("等位线间隔:");
        for (i = 0; i < count1; i++)
        {
            printf(" %d", V1[i]);
        }
    }

    return 0;
}

int readdata2(FILE* file)
{

#if 0
    double a;        // 电极厚度
    int n;           // 电极个数
    double* z;       // 电极间距
    int* N;          // 相邻电极间划分步长
    double* V;       // 电极电位
    double* VI;      //含鞍点电极电位
    double* r;       // 电极内孔半径
    //double r2;       // 从电极内孔径边沿到封闭边界处的径向距离
    int* M;      // r1范围内等步长划分的网格数,从电极内孔径边沿到封闭边界处的径向距离
    int NST;         // 输出打印空间电位时网格点间隔数
    int INS;         // 轴上电位做等距插值时步长数
    double e;        // 迭代精度
    int m_V;         // 扫描电位个数
    //double* V_scan;  // 制定扫描点位时暂存扫描电位
    //int I,m = 0;    // m为格点列数，I为格点行数
    double V1;       // 要求扫描等电位线的电位间隔

    int i;
    m_V = 0;
#endif
    tmp = 0;
    //int count = 0;
    int i;
    //fopen_s(&file, filepath, "rt");

    fscanf_s(file, "delta = %lf mm; n = %d;\n", &delta, &n); // 电极宽度

    dz = (double*)malloc((n+1) * sizeof(double)); // 相邻电极间距离
    fscanf_s(file, "Zi =");
    for (i = 0; i < n+1; i++)
    {
        fscanf_s(file, " %lf", &dz[i]);
    }
    fscanf_s(file, "mm;\n");

    _N = (int*)malloc((n + 1) * sizeof(int)); // 相邻电极间步长
    fscanf_s(file, "Ni =");
    for (i = 0; i < n + 1; i++)
    {
        fscanf_s(file, " %d", &_N[i]);
    }
    fscanf_s(file, ";\n");

    V = (double*)malloc(n * sizeof(double)); // 电极电位
    fscanf_s(file, "Vi =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(file, " %lf", &V[i]);
    }
    fscanf_s(file, "V;\n");

    VI = (double*)malloc(n * sizeof(double)); // 含鞍点电极电位
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
        V1 = (int*)malloc((n + 10) * sizeof(int)); // 等位线
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

    // 打印读取到的数据，可根据需要自行调整
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
        printf("等位线间隔：%d V\n", V1[0]);
    }
    else
    {
        printf("等位线间隔:");
        for (i = 0; i < count1; i++)
        {
            printf(" %d", V1[i]);
        }
    }
    return 0;
}



void writedata(FILE* file, Grid_Array** grid1, Grid_Array** grid2, Iteration_Process* head1,Iteration_Process* head2)
{
    //输出基本信息
    int i;
    fprintf(file, "像管类型: %d类像管 \n", mode);
    fprintf(file, "电极尺寸δ = %.1lf mm\n", delta);
    fprintf(file, "电极个数n = %d\n", n);
    fprintf(file, "电极间距Zi = ");
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
    fprintf(file, "轴向划分要求Ni = ");
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

    fprintf(file, "电极电压Vi = ");
    for (i = 0; i < n; i++)
    {
        fprintf(file, "%.1lf ", V[i]);
    }
    fprintf(file, "V\n");

    fprintf(file,"含鞍点电压Vi = ");
    for (i = 0; i < n; i++)
    {
        fprintf(file, "%.1lf ", VI[i]);
    }
    fprintf(file, "V\n");

    if (mode == 2)
    {
        fprintf(file, "垂直间距r = ");
        for (i = 0; i < n - 1; i++)
        {
            fprintf(file, " %.1lf", dr[i]);
        }
        fprintf(file, "mm;\n");
    }
    else
    {
        fprintf(file, "垂直间距r1 = %lf mm, r2 = %lf mm\n", r1, r2);
    }
    

    if (mode == 2)
    {
        fprintf(file, "径向划分要求M = ");
        for (i = 0; i < n - 1; i++)
        {
            fprintf(file, " %d", _M[i]);
        }
        fprintf(file, ";\n");
    }
    else
    {
        fprintf(file, "径向划分要求M1 = %d, M2 = %d\n", M1, M2);
    }
    


    fprintf(file, "迭代精度ε = %.7lf V\n", epsilon);

    fprintf(file, "输出电位间隔数NST = %d\n", NST);

    fprintf(file, "轴上电位插值步长数INS = %d\n", INS);

    if (tmp == 3)
    {
        fprintf(file, "等宽扫描等位线间隔：%d V\n", V1[0]);
    }
    else
    {
        fprintf(file, "给定扫描等位线间隔:");
        for (i = 0; i < count1; i++)
        {
            fprintf(file, " %d", V1[i]);
        }
    }
    fprintf(file, "\n");
    fprintf(file, "总网格数M = %d, N = %d\n\n", M, N);

    //输出网格点坐标
    //计算轴向坐标z
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

    fprintf(file, "网格点坐标：\n         轴向坐标 ");
    for (i = 0; i < N; i++)
    {
        fprintf(file, "%6.2lf  ", z[i]);
    }
    fprintf(file, "\n径向坐标\n");
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

    //打印无鞍点网格
    fprintf(file, "\n\n无鞍点网格：\n");
    //打印迭代链表
    fprintf(file, "\n迭代相关信息：\n");
    fprintf(file, "迭代次数: %d 次\n", iteration_times_1);
    fprintf(file, "每次迭代的加速因子ω，迭代精度（最大残差）和平均残差:\n");
    fprintf(file, "迭代轮次          加速因子ω          迭代精度          平均残差:\n");
    Iteration_Process* p;
    p = head1;
    p = p->next;
    for (i = 0; i < iteration_times_1; i++)
    {
        fprintf(file, "%2d轮%2d次             %.7lf         %.7lf         %.7lf\n", p->round,p->times, p->omega_r, p->max_res, p->avg_res);
        if (p->next != NULL)
        {
            p = p->next;
        }
        
    }


    //打印网格电位
    fprintf(file, "\n网格电位：\n");
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            fprintf(file, "%11.7f ", grid1[i][j].voltage); 
        }
        fprintf(file,"\n");
    }
            
    //打印有鞍点网格
    fprintf(file, "\n\n有鞍点网格：\n");
    //打印迭代链表
    fprintf(file, "\n迭代相关信息：\n");
    fprintf(file, "迭代次数: %d 次\n", iteration_times_1);
    fprintf(file, "每次迭代的加速因子ω，迭代精度（最大残差）和平均残差:\n");
    fprintf(file, "迭代次数          加速因子ω          迭代精度          平均残差:\n");
    Iteration_Process* q;
    q = head2;
    q = q->next;
    for (i = 0; i < iteration_times_2; i++)
    {
        fprintf(file, "%2d轮%2d次             %.7lf         %.7lf         %.7lf\n", q->round, q->times, q->omega_r, q->max_res, q->avg_res);
        if (q->next != NULL)
        {
             q = q->next;
        }

    }


    //打印网格电位
    fprintf(file, "\n网格电位：\n");
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            fprintf(file, "%11.7f ", grid1[i][j].voltage);
        }
        fprintf(file, "\n");
    }
}
