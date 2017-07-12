#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h> 

#define HROW           640  //校验矩阵行数行数 未知  
#define M              960  //列数             未知
#define K              320   //生成矩阵行数     未知
#define ROWLEN     sizeof(struct rownode)
#define COLLEN     sizeof(struct colnode)
#define BLOCKNUMBER        1000
#define PI                   3.1415926f
#define ITERATIMES           5

struct rownode
{
	int rownum,colnum;
	double q;
	struct rownode * rownext;
};

struct colnode
{
    int rownum,colnum;
	double r;
	struct colnode * colnext;
	
};

double snr[10]={1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2};             //设定信噪比
int l;
int hmatrix[HROW][M]={0};    //定义为全局变量
int gmatrix[K][M]={0};
struct rownode * rowhead[HROW];
struct colnode * colhead[M];

void main()
{
   void generate(int gseq[]);
   void addgaussnoise(double addcode[]);
   void decoder(double designal[]);
   int bercounter(int array1[],int array2[]);
   int percounter(int array1[],int array2[]);
   double Gauss_rand();
   double contanh(double y);
  
   int i,j,n;
   int benum=0;
   int penum=0;
   double ber[10]={0};
   double per[10]={0};
   int trancode[M];           //传送的双极性序列
   int decode[M];             //译码得出的0-1序列，就是将译码器返回的double型的signal序列转成int型
   double signal[M];           //传送和接受的实际信号 
   struct rownode * p1,* p2;
   struct colnode * q1,* q2;
   FILE * fph,* fpg,* fp_ber,* fp_per;

   //****************************从文件读入生成矩阵和校验矩阵*******************************************//
   fpg=fopen("Matrix_g2.txt","r");
   //fpg=fopen("gmatrix.txt","r");
   fph=fopen("Matrix_h2.txt","r");
   //fph=fopen("hmatrix.txt","r");
   fp_ber=fopen("ber.txt","w");
   fp_per=fopen("per.txt","w");
   if(fpg==NULL)
   { 
	   printf("生成矩阵不在指定的文件中!\n");
	   exit(0);
   }
   if(fph==NULL)
   { 
	   printf("校验矩阵不在指定的文件中!\n");
	   exit(0);
   } 

   for (i=0;i<K;i++)
   {
	   for (j=0;j<M;j++)
	   {
		   fscanf(fpg,"%d",&gmatrix[i][j]);
	   }
   }
   for (i=0;i<HROW;i++)
   {
	   for (j=0;j<M;j++)
	   {
		   fscanf(fph,"%d",&hmatrix[i][j]);
	   }
   }

   for (j=0;j<HROW;j++)                 //根据行建立链表
   {
	   n=0;
	   p1=(struct rownode *)malloc(ROWLEN);
	   p2=p1;
	   for (i=0;i<M;i++)
	   {
		   if (hmatrix[j][i]==1)
		   {
			   p1->rownum=j;
			   p1->colnum=i;
			   //p1->q=p[i];             //初始化信息
			   n=n+1;
			   if (n==1)rowhead[j]=p1;
			   else p2->rownext=p1;
			   p2=p1;
			   p1=(struct rownode *)malloc(ROWLEN);
		   }
	   }
	   p2->rownext=NULL;
	   free(p1);
    }                                  //建立行链表，第j行的头指针为rowhead[j]

   for (i=0;i<M;i++)                  //按列建立链表
   {
	   n=0;
	   q1=(struct colnode *)malloc(COLLEN);
	   q2=q1;
	   for (j=0;j<HROW;j++)
	   {
		   if (hmatrix[j][i]==1)
		   {
			   q1->rownum=j;
			   q1->colnum=i;
			   n=n+1;
			   if (n==1)colhead[i]=q1;
			   else q2->colnext=q1;
			   q2=q1;
			   q1=(struct colnode *)malloc(COLLEN);
		   }
	   }
	   q2->colnext=NULL;
	   free(q1);
	}                                   //建立列链表，第i列的头指针为colhead[i]

   srand((int)time(0));     //产生随机数种子
  
   for (l=0;l<10;l++)
   {

   benum=0;
   penum=0;
   
   for(j=0;j<BLOCKNUMBER;j++)
   {
      //************************产生一组随机的发送码序列bpsk调制信号***************************************//
      if ((j%20)==0)
      {
		  printf("haha!");
      }
	  generate(trancode);
      for (i=0;i<M;i++)
	  {
	      trancode[i]=trancode[i]*2-1;  //bpsk
	      signal[i]=(double)trancode[i];       //将信号转成double型
	      //printf("%f ",signal[i]);
	  }
      //printf("\n以上是产生的bpsk信号.\n");

      //********************************add gaussnoise******************************************//
      addgaussnoise(signal);
      /*for (i=0;i<M;i++)
	  {
	      printf("%f ",signal[i]);
	  }
      printf("\n以上是加上高斯白噪声后的信道信息\n");*/
      //*************************************decode*********************************************//
      decoder(signal);

      //***********************************count BER and PER************************************//
      //for (i=0;i<M;i++)
      //{
      //	   printf("%f ",signal[i]);
      //}
      //printf("\n");

      for (i=0;i<M;i++)
	  {
	   
	      decode[i]=(int)signal[i];//译码结果转成整型
	  }
      
      /*for (i=0;i<M;i++)
	  {
	      printf("%d",decode[i]);
	  }
      printf("\n以上是译码比特\n");*/
      benum=benum+bercounter(trancode,decode);
      penum=penum+percounter(trancode,decode);
      //printf("\nbenum=%d, penum=%d\n",benum,penum);
      //printf("\n");
   }
   ber[l]=(double)benum/(double)(BLOCKNUMBER*M);
   per[l]=(double)penum/(double)BLOCKNUMBER;
   //printf("\n误比特率为：%f\n误分组率为：%f\n",ber,per);
   fprintf(fp_ber,"%f,%f  ",ber[l],snr[l]);
   fprintf(fp_per,"%f,%f  ",per[l],snr[l]);
   //printf("%d\n",l);
   }
   fclose(fpg);
   fclose(fph);
   fclose(fp_ber);
   fclose(fp_per);
}

void generate(int gseq[])
{
	int ranseq[K]={0};
	int i,j;
	
	
	for (i=0;i<K;i++)
	{
		ranseq[i]=rand()%2;
		//printf("%d",ranseq[i]);
	}
	//printf("\n以上是产生的随机信源序列:\n");
	
	for (i=0;i<M;i++)
	{
		gseq[i]=0;
		
		for (j=0;j<K;j++)
		{
			gseq[i]=gseq[i]+ranseq[j]*gmatrix[j][i];
		}
		gseq[i]=gseq[i]%2;
		//transeq[i]=gseq[i]*2-1;//pam to bpsk
		//printf("%d",gseq[i]);
	}
	/*printf("\nNow here are the bpsk signal:\n");
	for (i=0;i<M;i++)
	{
		printf("%d ",transeq[i]);
	}*/
	//printf("\n以上是生成的码字\n");
}

void addgaussnoise(double addcode[])
{
     int i;         
     double var=1/snr[l];
	 double noise;                                 
	                                               
	 for (i=0;i<M;i++)
	 {
		 noise=(double)sqrt(var)*Gauss_rand();     //产生均值为0，方差为var的白噪声
		 addcode[i]=addcode[i]+noise;             //加上高斯白噪声
	 }
	                                         
}

void decoder(double designal[])
{ 
    int n,k,j,i;
	double temp=1;
	double var=1/snr[l];
	double p[M];
	double deq[M];

	struct rownode * p1;
	struct colnode * q1;
	//int check[HROW][M]={0};
	//double rcheck[HROW][M]={0};

	//for (i=0;i<5;i++)
	//{
	//	var[i]=1/snr[i];              //由于这里是双极性，不考虑调制波形
	//}

    for (i=0;i<M;i++)
    {
		p[i]=0-2*designal[i]/var;     //信道信息
		//printf("%f ",p[i]);
    }

	for (i=0;i<HROW;i++)              //初始化信息
	{
		p1=rowhead[i];
		while(p1!=NULL)
		{
			p1->q=p[p1->colnum];
			p1=p1->rownext;
		}
	}
	/*for (j=0;j<HROW;j++)                //检验校验矩阵是否构造正确
	{
		p1=rowhead[j];
		while(p1!=NULL)
		{
			check[p1->rownum][p1->colnum]=1;
			p1=p1->rownext;
		}

	}
	printf("\n");
	for (j=0;j<HROW;j++)
	{
		for (i=0;i<M;i++)
		{
			printf("%d ",check[j][i]);
		}
		printf("\n");
	}*/
    


	/*for (i=0;i<M;i++)                 //输出校验矩阵，检验校验矩阵是否构造正确
	{
		q1=colhead[i];
		while (q1!=NULL)
		{
			check[q1->rownum][q1->colnum]=1;
			q1=q1->colnext;
		}
	}
    printf("\n");
	for (j=0;j<HROW;j++)
	{
		for (i=0;i<M;i++)
		{
			printf("%d ",check[j][i]);
		}
		printf("\n");
	}*/

/******************开始迭代处理****************************/
    for (k=0;k<ITERATIMES;k++)
	{
	
	
	for (i=0;i<M;i++)                               //校验节点消息处理
	{
		q1=colhead[i];
		while (q1!=NULL)
		{
            p1=rowhead[q1->rownum];
			temp=1;
			while (p1!=NULL)
			{   
				if (p1->colnum!=i)
				{
					temp=temp*tanh(0.5*(p1->q));
				}				
                p1=p1->rownext;
			}

			q1->r=2*contanh(temp);
			
			q1=q1->colnext;
		}
        
	}                                                //校验节点消息处理完成

/*	for (i=0;i<M;i++)                                //显示各个节点的r值
	{
		q1=colhead[i];
		while (q1!=NULL)
		{
			rcheck[q1->rownum][q1->colnum]=q1->r;
			q1=q1->colnext;
		}
	}
    printf("\n");
	for (j=0;j<HROW;j++)
	{
		for (i=0;i<M;i++)
		{
			printf("%f ",rcheck[j][i]);
		}
		printf("\n");
	}
*/
	for (j=0;j<HROW;j++)                                //变量节点消息处理
	{
		p1=rowhead[j];
		while (p1!=NULL)
		{
            q1=colhead[p1->colnum];
			temp=0;
			while (q1!=NULL)
			{
				if (q1->rownum!=j)
				{
					temp=temp+(q1->r);
				}
				q1=q1->colnext;
			}
            p1->q=p[p1->colnum]+temp;
			p1=p1->rownext;
		}
	}                                              //变量节点消息处理完成
    //printf("\n");
    for (i=0;i<M;i++)                              //译码判决
    {
		q1=colhead[i];
		temp=0;
		while (q1!=NULL)
		{
			temp=temp+(q1->r);
			q1=q1->colnext;
		}
		deq[i]=p[i]+temp;
		//printf("%f ",deq[i]);
		if (deq[i]>0)
		{
			designal[i]=0;
		} 
		else
		{
			designal[i]=1;
		}
		
    }                                              //译码判决完成

	n=0;

	for (j=0;j<HROW;j++)
	{
		temp=0;
		for (i=0;i<M;i++)
		{
			temp=temp+hmatrix[j][i]*designal[i];
		}
		if (temp)
		{
			n=n+1;
		}
	}
	if(n==0)break;                        //如果校正子都为0，则译码结束，跳出迭代
    }


}

/*void decoder(double designal[])
{
     int i,j,k;
	 int n,m;
	 double temp1=1;
	 double temp2=1;
	 double snr=2;
	 double var=1/snr;
	 int hmatrix[15][15]={0,0,0,0,0,0,0,1,1,0,1,0,0,0,1,1,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,1,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,1,0};
     double q0[15][15];
	 double q1[15][15];
	 double r0[15][15];
	 double r1[15][15];
     double p0[15],p1[15];
	 double deq0[15],deq1[15];
     
	 
	 printf("\n");
     for (i=0;i<15;i++)
     {
		 for (j=0;j<15;j++)
		 {
			 printf("%d ",hmatrix[i][j]);
		 }
		 printf("\n");
     }
	 printf("\n");

     for (i=0;i<15;i++)                                //初始化
     {
		 //p0[i]=1/(1+exp(2*designal[i]/var));
		 //p1[i]=1-p0[i];

		 p0[i]=1/(exp((designal[i]+1)*(designal[i]+1)/(2*var))*sqrt(2*PI*var));
		 p1[i]=1/(exp((designal[i]-1)*(designal[i]-1)/(2*var))*sqrt(2*PI*var));
		 p0[i]=p0[i]/(p0[i]+p1[i]);
		 p1[i]=1-p0[i];
		 printf("%f %f	",p0[i],p1[i]);
		 for (j=0;j<15;j++)
		 {
			 q0[i][j]=p0[i];
			 q1[i][j]=p1[i];
		 }
     }
	 printf("\n");

     for (j=0;j<15;j++)
     {
		 for (i=0;i<15;i++)
		 {
			 if (hmatrix[j][i]==0)
			 {
				 q0[i][j]=0;
				 q1[i][j]=0;
				 r0[j][i]=1;
				 r1[j][i]=1;
			 }
		 }
     }
	 
	 for (k=0;k<ITERATIMES;k++)
	 {
		 for (j=0;j<15;j++)                            //校验节点处理开始
		 {
			 for (i=0;i<15;i++)
			 {
                 if ((r0[j][i]!=1)&&(r1[j][i]!=1))
                 {
					 temp1=1;
					 for (m=0;m<15;m++)
					 {
						 if (m!=i)
						 {
                            temp1=temp1*(1-2*q1[m][j]);
						 }
					 }
					 r0[j][i]=0.5*temp1+0.5;
					 r1[j][i]=1-r0[j][i];
                 }
			 }
		 }                                             //校验节点处理完成
	     

	     for (i=0;i<15;i++)                            //变量节点处理开始
	     {
			 for (j=0;j<15;j++)
			 {
				 if ((q0[i][j]!=0)&&(q1[i][j]!=0))
				 {
					 temp1=1;
					 temp2=1;
					 for (m=0;m<15;m++)
					 {
						 if (m!=j)
						 {
							 temp1=temp1*r0[m][i];
							 temp2=temp2*r1[m][i];
						 }
					 }
					 q0[i][j]=p0[i]*temp1;
					 q1[i][j]=p1[i]*temp2;
					 q0[i][j]=q0[i][j]/(q0[i][j]+q1[i][j]);
					 q1[i][j]=1-q0[i][j];
				 }
			 }
		 }                                              //变量节点处理完成
         
		 
		 for (i=0;i<15;i++)                             //译码判决
         {
			 temp1=1;
			 temp2=1;
			 for (j=0;j<15;j++)
			 {
				 temp1=temp1*r0[j][i];
                 temp2=temp2*r1[j][i];
			 }
			 deq0[i]=p0[i]*temp1;
			 deq1[i]=p1[i]*temp2;
			 deq0[i]=deq0[i]/(deq0[i]+deq1[i]);
			 deq1[i]=1-deq0[i];
         }

		 printf("\n"); 
		 for (i=0;i<15;i++)
		 {
			 printf("%f %f	",deq0[i],deq1[i]);
		 }
         printf("\n以上是硬判决信息，左为0右为1\n"); 

		 for (i=0;i<15;i++)
		 {
			 if (deq1[i]>deq0[i])
			 {
				 designal[i]=1;
			 }
			 else
				 designal[i]=0;
		 }

         n=0;
		 
		 for (i=0;i<15;i++)
		 {
			 temp1=0;
		     for (j=0;j<15;j++)
		     {
                 temp1=temp1+hmatrix[i][j]*designal[j];
		     }
			 if (temp1)
			 {
				 n=n+1;
			 }
		     	 
		 }
		 if (n==0)break;                      //如果校正子都为0，则译码结束，跳出迭代
	 }
 
}*/

int bercounter(int array1[],int array2[])
{
    int n=0;
	int i;
    for (i=0;i<M;i++)
    {
        array1[i]=(array1[i]+1)/2;
		if (array1[i]!=array2[i])
		{
            n=n+1;
		}
    }
	return n;
}

int percounter(int array1[],int array2[])
{
    int n=0;
	int i;
	for (i=0;i<M;i++)
	{
		if (array1[i]!=array2[i])
		{
			n=1;
		}
	}
	return n;
}

double Gauss_rand()
{
	double random_num;
	double u1,u2;
	u1=((double)rand()+1.f)/32768.f;
	u2=((double)rand()+1.f)/32768.f;
	random_num=(double)(sqrt(-2.*log(u1))*cos(2*PI*u2));
	return random_num;
}

double contanh(double y)
{
	double value;
	value=0.5*log((1+y)/(1-y));
	return value;
}
