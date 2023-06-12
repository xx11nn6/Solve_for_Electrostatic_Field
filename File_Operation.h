#define FILE_OPERATION_H  //包括文件的读写，读文件由葛军韬负责，写文件由魏千怀负责，扫描点电压和坐标由李煜翔负责
#include "Grid_Array.h"
#include <stdio.h>
#include <stdlib.h>


int rescan_1(int V, Grid_Array** grid);
void scan_all_1(int n, int M1, int M2, int M, int N, Grid_Array** grid, double z0, double r0, double* dz, double delta, int tmp, int* V1, int count1, FILE* file);
int rescan_2(int V, Grid_Array** grid, int n, double* dz, double* dr, double delta);
void scan_all_2(int n, int M, int N, Grid_Array** grid, double* dz, double delta, int* _M, double* dr, int* _N, int tmp, int* V1, int count1, FILE* file);

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
    return 0;
}

int readdata2(FILE* file)
{
    tmp = 0;
    //int count = 0;
    int i;
    //fopen_s(&file, filepath, "rt");

    fscanf_s(file, "delta = %lf mm; n = %d;\n", &delta, &n); // 电极宽度

    dz = (double*)malloc((n + 1) * sizeof(double)); // 相邻电极间距离
    fscanf_s(file, "Zi =");
    for (i = 0; i < n + 1; i++)
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
    return 0;
}



void writedata(FILE* file, Grid_Array** grid1, Grid_Array** grid2, Iteration_Process* head1, Iteration_Process* head2)
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

    fprintf(file, "含鞍点电压Vi = ");
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
        z[i] = z[i - 1] + grid1[M - 1][i - 1].h2;
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
                    fprintf(file, "%6.2lf            ", grid1[i][M - 1].r);
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
        fprintf(file, "%2d轮%2d次             %.7lf         %.7lf         %.7lf\n", p->round, p->times, p->omega_r, p->max_res, p->avg_res);
        if (p->next != NULL)
        {
            p = p->next;
        }
    }

    //打印轴上电压
    if (mode == 1)
    {
        double z_1, v_1, z_0, v_0;
        double z_pace = z0 / INS;

        z_0 = 0;
        v_0 = 0;
        fprintf(file, "\n无鞍点像管轴上点坐标与电位：\n");
        fprintf(file, "坐标r             电位V\n");
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

        fprintf(file, "\n无鞍点像管轴上点坐标与电位：\n");
        fprintf(file, "坐标r             电位V\n");
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
        //cout << z_0 << "电压" << 100 << endl;
    }


    //打印网格电位
    fprintf(file, "\n无鞍点网格电位：\n");
    for (int i = 0; i < M; i += NST)
    {
        for (int j = 0; j < N; j += NST)
        {
            fprintf(file, "%11.7f ", grid1[i][j].voltage);
        }
        fprintf(file, "\n");
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

    //打印轴上电压
    if (mode == 1)
    {
        double z_1, v_1, z_0, v_0;
        double z_pace = z0 / INS;

        z_0 = 0;
        v_0 = 0;
        fprintf(file, "\n无鞍点像管轴上点坐标与电位：\n");
        fprintf(file, "坐标r             电位V\n");
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

        fprintf(file, "\n无鞍点像管轴上点坐标与电位：\n");
        fprintf(file, "坐标r             电位V\n");
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
    //打印网格电位
    fprintf(file, "\n有鞍点网格电位：\n");
    for (int i = 0; i < M; i += NST)
    {
        for (int j = 0; j < N; j += NST)
        {
            fprintf(file, "%11.7f ", grid2[i][j].voltage);
        }
        fprintf(file, "\n\n");
    }


    //打印扫描点电压及坐标
    if (mode == 1)
    {
        fprintf(file, "扫描点电压及坐标:\n");
        fprintf(file, "无鞍点时:\n");
        scan_all_1(n, M1, M2, M, N, grid1, z0, r0, dz, delta, tmp, V1, count1, file);
        fprintf(file, "\n有鞍点时:\n");
        scan_all_1(n, M1, M2, M, N, grid2, z0, r0, dz, delta, tmp, V1, count1, file);
    }
    else
    {
        fprintf(file, "扫描点电压及坐标:\n");
        fprintf(file, "无鞍点时:\n");
        scan_all_2(n, M, N, grid1, dz, delta, _M, dr, _N, tmp, V1, count1, file);
        fprintf(file, "\n有鞍点时:\n");
        scan_all_2(n, M, M, grid2, dz, delta, _M, dr, _N, tmp, V1, count1, file);
    }

}



// 扫描电压
void scan_all_1(int n, int M1, int M2, int M, int N, Grid_Array** grid, double z0, double r0, double* dz, double delta, int tmp, int* V1, int count1, FILE* file)//画电位线,
{
    int i, j, k, num;
    double z_1, v_1, z_0, v_0, V_temp1, V_temp2, t;
    V_z = (double*)malloc(n * (M1 + M2) * sizeof(double));//存扫描电Z坐标
    V_r = (double*)malloc(n * (M1 + M2) * sizeof(double));//存扫描电R坐标
    zi = (double*)malloc(N * sizeof(double));//存每列到r轴距离
    ri = (double*)malloc(M * sizeof(double));//存每行径向距离
    zi[0] = 0;
    for (i = 1; i < N; i++)//计算z1
    {
        zi[i] = zi[i - 1] + grid[0][i - 1].h2;
    }
    for (i = 0; i < M; i++)//计算
    {
        ri[i] = grid[i][0].r;
    }

    for (i = 0; i < n * (M1 + M2); i++)//将所有扫描点的坐标初始化为0
    {
        V_z[i] = 0;
        V_r[i] = 0;
    }

    if (tmp == 3) {
        for (double sss = 0; sss < 100; sss += V1[0]) {
            k = rescan_1(sss, grid);//输入要扫描的电位，进行指定电位的扫描，返回k为扫描点总数
            for (int j = 0; j < k; j++)
            {
                fprintf(file, "扫描点电压：%d V     扫描点坐标：z= %f    r=%f \n", sss, V_z[j] / 10, V_r[j] / 10);
            }
        }
    }
    else {
        for (i = 0; i < count1; i++) {
            k = rescan_1(V1[i], grid);//输入要扫描的电位，进行指定电位的扫描，返回k为扫描点总）
            for (int j = 0; j < k; j++)
            {
                fprintf(file, "扫描点电压：%d V     扫描点坐标：z= %f    r=%f \n", V1[i], V_z[j] / 10, V_r[j] / 10);
            }
        }
    }
}

int rescan_1(int V, Grid_Array** grid)//扫描//
{
    int i, j, k = 0;
    double r, z;
    for (i = 0; i < M; i++) //行扫描
    {
        for (j = 0; j < N - 1; j++)
        {

            if (grid[i][j].voltage == V) {//若点电位与扫描电位相等，将该点的坐标保存，10为图像的放大倍数
                V_z[k] = 10 * zi[j];
                V_r[k] = 10 * ri[i];
                k++;
            }
            if ((grid[i][j].voltage <= V && grid[i][j + 1].voltage > V) || (grid[i][j].voltage >= V && grid[i][j + 1].voltage < V) || (grid[i][j].voltage > V && grid[i][j + 1].voltage <= V)\
                || (grid[i][j].voltage < V && grid[i][j + 1].voltage >= V))//若扫描点与后一点的电位不相等，则进行插值
            {
                z = (V - grid[i][j].voltage) / (grid[i][j + 1].voltage - grid[i][j].voltage) * (zi[j + 1] - zi[j]) + zi[j];//线性插值
                V_z[k] = 10 * z;
                V_r[k] = 10 * ri[i];
                k++;
            }
        }
    }
    for (j = 0; j < N; j++)//列扫描，具体步骤同上
    {
        for (i = 0; i < M - 1; i++)
        {
            if (grid[i][j].voltage == V) {//若点电位与扫描电位相等，将该点的坐标保存，10为图像的放大倍数
                V_z[k] = 10 * zi[j];
                V_r[k] = 10 * ri[i];
                k++;
            }
            if ((grid[i][j].voltage <= V && grid[i + 1][j].voltage > V) || (grid[i][j].voltage >= V && grid[i + 1][j].voltage < V) || (grid[i][j].voltage > V && grid[i + 1][j].voltage <= V)\
                || (grid[i][j].voltage < V && grid[i + 1][j].voltage >= V))//若扫描点与后一点的电位不相等，则进行插值
            {


                r = (V - grid[i][j].voltage) / (grid[i + 1][j].voltage - grid[i][j].voltage) * (ri[i + 1] - ri[i]) + ri[i];//线性插值
                V_z[k] = 10 * zi[j];
                V_r[k] = 10 * r;
                k++;

            }
        }
    }
    return k;
}
void scan_all_2(int n, int M, int N, Grid_Array** grid, double* dz, double delta, int* _M, double* dr, int* _N, int tmp, int* V1, int count1, FILE* file)//画电位线,
{
    int i, j, k, num;
    double z_temp1, z_temp2, V_temp1, V_temp2, r_0, z_0, t, s, v_2, z_2, z_1, v_1;
    V_z = (double*)malloc(n * (M) * sizeof(double));//存扫描电Z坐标
    V_r = (double*)malloc(n * (M) * sizeof(double));//存扫描电R坐标
    zi = (double*)malloc(N * sizeof(double));//存每列到r轴距离
    ri = (double*)malloc(M * sizeof(double));//存每行径向距离
    zi[0] = 0;
    grid[M - 1][0].h2 = dz[0] / _N[0];
    for (i = 1; i < N; i++)//计算z1
    {
        zi[i] = zi[i - 1] + grid[M - 1][i - 1].h2;
    }
    for (i = 0; i < M; i++)//计算r1
    {
        ri[i] = grid[i][N - 2].r;
    }
    for (i = 0; i < n * (M); i++)//将所有扫描点的坐标初始化为0
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

    if (tmp == 3) {
        for (i = V1[0]; i < 100; i = i + V1[0]) {

            k = rescan_2(i, grid, n, dz, dr, delta);//输入要扫描的电位，进行指定电位的扫描，返回k为扫描点总数
            for (int j = 0; j < k; j++)
            {
                fprintf(file, "扫描点电压：%d V     扫描点坐标：z= %f    r=%f \n", i, V_z[j] / 6, V_r[j] / 6);
            }

        }
    }
    else {

        for (i = 0; i < count1; i++) {

            k = rescan_2(V1[i], grid, n, dz, dr, delta);//输入要扫描的电位，进行指定电位的扫描，返回k为扫描点总数
            for (int j = 0; j < k; j++)
            {
                fprintf(file, "扫描点电压：%d V     扫描点坐标：z= %f    r=%f \n", i, V_z[j] / 6, V_r[j] / 6);
            }
        }
    }
    system("pause");
    closegraph();



}
int rescan_2(int V, Grid_Array** grid, int n, double* dz, double* dr, double delta)//扫描//
{
    int i, j, k = 0;
    int p, q;
    double r, z, t, s;

    for (i = 0; i < M; i++) //行扫描
    {
        for (j = 0; j < N - 1; j++)
        {
            if (grid[i][j].voltage == V) {//若点电位与扫描电位相等，将该点的坐标保存，10为图像的放大倍数
                V_z[k] = zi[j];
                V_r[k] = ri[i];
                k++;

            }
            if ((grid[i][j].voltage <= V && grid[i][j + 1].voltage > V) || (grid[i][j].voltage >= V && grid[i][j + 1].voltage < V) || (grid[i][j].voltage > V && grid[i][j + 1].voltage <= V)\
                || (grid[i][j].voltage < V && grid[i][j + 1].voltage >= V))//若扫描点与后一点的电位不相等，则进行插值
            {
                z = (V - grid[i][j].voltage) / (grid[i][j + 1].voltage - grid[i][j].voltage) * (zi[j + 1] - zi[j]) + zi[j];//线性插值
                V_z[k] = z;
                V_r[k] = ri[i];
                k++;

            }
        }
    }
    for (j = 0; j < N; j++)//列扫描，具体步骤同上
    {
        for (i = 0; i < M - 1; i++)
        {
            if (grid[i][j].voltage == V) {//若点电位与扫描电位相等，将该点的坐标保存，10为图像的放大倍数
                V_z[k] = zi[j];
                V_r[k] = ri[i];
                k++;
            }
            if ((grid[i][j].voltage <= V && grid[i + 1][j].voltage > V) || (grid[i][j].voltage >= V && grid[i + 1][j].voltage < V) || (grid[i][j].voltage > V && grid[i + 1][j].voltage <= V)\
                || (grid[i][j].voltage < V && grid[i + 1][j].voltage >= V))//若扫描点与后一点的电位不相等，则进行插值
            {
                r = (V - grid[i][j].voltage) / (grid[i + 1][j].voltage - grid[i][j].voltage) * (ri[i + 1] - ri[i]) + ri[i];//线性插值
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
                }
                t = t + dz[0];
            }
            else if (p == 1) {
                if ((V_z[i] <= (t + dz[p])) && (V_r[i] <= dr[0]) && (V_z[i] > t))
                {
                    V_z[q] = 8 * V_z[i];
                    V_r[q] = 8 * V_r[i];
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
                    q++;
                }
            }

            else {
                if ((V_z[i] <= ((t + dz[p]))) && (V_r[i] <= ((s + dr[p - 1]))) && (V_z[i] > (t)))
                {
                    V_z[q] = 8 * V_z[i];
                    V_r[q] = 8 * V_r[i];
                    q++;

                }
                t = t + dz[p];
                s = s + dr[p - 1] + delta;
            }
        }
    }

    k = q;
    return k;//返回所扫描到的坐标个数
}

