//头文件调用
#include<stdio.h>
#include<stdlib.h>

//函数声明
int readdata1(const char* filepath);//一类像管读文件
int readdata2(const char* filepath);//二类像管读文件

//文件声明
FILE* handle;
errno_t err;

//全局变量
int type = 0;          //确定像管类型
int tmp = 0;           //确定等位线扫描方式
double delta;         //电极厚度(通用)
int n;                //电极个数(通用)
double* Zi;           //电极间距(通用)
int* Ni;              //相邻电极间划分步长(通用)
double* Vi;           // 电极电位(通用)
double* VI;           //含鞍点电极电位(通用)
double r1;            // 电极内孔半径(第一类)
int M1;               // r1范围内等步长划分的网格数(第一类)
double r2;            // 从电极内孔径边沿到封闭边界处的径向距离(第一类)
int M2;               //从电极内孔径边沿到封闭边界处的径向距离(第一类)
int* M;               // r范围内等步长划分的网格数(第二类)
double* r;            // 电极内孔半径(第二类)
double epsilon;       // 迭代精度(通用)
int NST;              // 输出打印空间电位时网格点间隔数(通用)
int INS;              // 轴上电位做等距插值时步长数(通用)
int* V1;              // 要求扫描等电位线的电位间隔或者电位值(通用)
int i;                //循环用参数
int count = 0;        //一个计数器，判断扫描等位线的条数

int readdata1(const char* filepath)
{


#if 0
    double delta;//电极厚度
    int n;//电极个数
    double* Zi;//电极间距
    int* Ni;//相邻电极间划分步长
    double* Vi;// 电极电位
    double* VI;//含鞍点电极电位
    double r1;// 电极内孔半径
    int M1;// r1范围内等步长划分的网格数
    double r2;// 从电极内孔径边沿到封闭边界处的径向距离
    int M2;//从电极内孔径边沿到封闭边界处的径向距离
    double epsilon;// 迭代精度
    int NST;// 输出打印空间电位时网格点间隔数
    int INS;// 轴上电位做等距插值时步长数
    int* V1;// 要求扫描等电位线的电位间隔

    int i;
    int m_V = 0;// 扫描电位个数
#endif
    
    // 读取数据
    //fopen_s(&handle, filepath, "rt");
    
    fscanf_s(handle, "delta = %lf mm; n = %d;\n", &delta, &n); // 电极宽度

    Zi = (double*)malloc((n + 1) * sizeof(double)); // 相邻电极间距离
    fscanf_s(handle, "Zi =");
    for (i = 0; i < n + 1; i++)
    {
        fscanf_s(handle, " %lf", &Zi[i]);
    }
    fscanf_s(handle, "mm;\n");

    Ni = (int*)malloc(n * sizeof(int)); // 相邻电极间步长
    fscanf_s(handle, "Ni =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(handle, " %d", &Ni[i]);
    }
    fscanf_s(handle, ";\n");

    Vi = (double*)malloc(n * sizeof(double)); // 电极电位
    fscanf_s(handle, "Vi =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(handle, " %lf", &Vi[i]);
    }
    fscanf_s(handle, "V;\n");

    VI = (double*)malloc(n * sizeof(double)); // 含鞍点电极电位
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
        V1 = (int*)malloc((n + 10) * sizeof(int)); // 等位线
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

    // 打印读取到的数据，可根据需要自行调整
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
        printf("等位线间隔：%d V\n", V1);
    }
    else
    {
        printf("等位线间隔:");
        for (i = 0; i < count; i++)
        {
            printf(" %d", V1[i]);
        }
    }

    // 释放动态分配的内存
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
    
    //fopen_s(&handle, filepath, "rt");

    fscanf_s(handle, "delta = %lf mm; n = %d;\n", &delta, &n); // 电极宽度

    Zi = (double*)malloc((n + 1) * sizeof(double)); // 相邻电极间距离
    fscanf_s(handle, "Zi =");
    for (i = 0; i < n + 1; i++)
    {
        fscanf_s(handle, " %lf", &Zi[i]);
    }
    fscanf_s(handle, "mm;\n");

    Ni = (int*)malloc((n + 1) * sizeof(int)); // 相邻电极间步长
    fscanf_s(handle, "Ni =");
    for (i = 0; i < n + 1; i++)
    {
        fscanf_s(handle, " %d", &Ni[i]);
    }
    fscanf_s(handle, ";\n");

    Vi = (double*)malloc(n * sizeof(double)); // 电极电位
    fscanf_s(handle, "Vi =");
    for (i = 0; i < n; i++)
    {
        fscanf_s(handle, " %lf", &Vi[i]);
    }
    fscanf_s(handle, "V;\n");

    VI = (double*)malloc(n * sizeof(double)); // 含鞍点电极电位
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
        V1 = (int*)malloc((n + 10) * sizeof(int)); // 等位线
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

    // 打印读取到的数据，可根据需要自行调整
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
        printf("等位线间隔：%d V\n", V1);
    }
    else
    {
        printf("等位线间隔:");
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
    if ((err = fopen_s(&handle, "D:\\CADTest\\2023\\2\\测试.txt", "rt")) != 0)
    {
        printf("找不到文件！错误码：%d\n", err);
        return -1;
    }
    
    
    fopen_s(&handle, "D:\\CADTest\\2023\\2\\测试.txt", "rt");
    if (fscanf_s(handle, "%d\n", &type) == 1 && type == 2)
    {
        readdata2("D:\\CADTest\\2023\\2\\测试.txt");
        fclose(handle);
    }
    else
    {
        readdata1("D:\\CADTest\\2023\\2\\测试.txt");
        fclose(handle);
    }
    return 0;
}
