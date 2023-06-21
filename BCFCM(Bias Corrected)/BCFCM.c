#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#define fmp 2
#define Sr 9.0
#define alpha 2.0
#define e 0.5

int main()
{
	
	int i,j,k,n,d,r,c,p;srand(time(0));char str[10];
	int m,*l=NULL,pos,f=0; float *avg=NULL,*center=NULL;int fm;
	float fmvsum=0.0,sum =0.0,mini;	int flag = 1;int win[10]; float Xpsum=0.0;int q;int in,jn;
    int lx; int count = 0;float Sig,Sij;int h;float scar=0.0;
	float sigma=0.0;

    int s=0; int g;
	int **a=NULL,**b=NULL;
	float ***dist=NULL;
	float ***temp=NULL;
	
	FILE *fp1 = fopen("m21_gray_gaus_hist.txt","r");
	FILE *fp2 = fopen("m21_gray_gaus_histbcfcm.pgm","wb");
	
	fscanf(fp1,"%s %d %d %d",str,&c,&r,&p);
	fprintf(fp2,"%s %d %d %d\n",str,c,r,p);
	
	printf("\n\nEnter the value of k: ");
	scanf("%d",&k);                           //Input no. of clusters
	
	
	a=(int **)calloc(r+2,sizeof(int *));
	for(i=0;i<r+2;i++)
	{
		a[i]=(int *)calloc(c+2,sizeof(int));
	}
	
	b=(int **)calloc(r+2,sizeof(int *));
	for(i=0;i<r+2;i++)
	{
		b[i]=(int *)calloc(c+2,sizeof(int));
	}
	
	dist = (float***)calloc(r,sizeof(float**));
	for(i=1;i<=r;i++)
	{
		dist[i]=(float**)calloc(c,sizeof(float*));
		for(j=1;j<=c;j++)
		{
			dist[i][j]=(float*)calloc(k,sizeof(float));
		}
	}
  
  	temp = (float***)calloc(r,sizeof(float**));
	for(i=1;i<=r;i++)
	{temp[i]=(float**)calloc(c,sizeof(float*));
		for(j=1;j<=c;j++)
		{
			temp[i][j]=(float*)calloc(k,sizeof(float));
		}
	}
	
    //*************************Input matrix of image************************************//	
    
   for(i=1;i<=r;i++)	
	{
		for(j=1;j<=c;j++)	                    //to read the input image
		{
			fscanf(fp1,"%d",&a[i][j]);
        }
	}
	//*************************matrix padding************************************//
	
	for(i=0;i<r+2;i++)	
	{
		a[i][0]=a[i][1];
        a[i][c+1]=a[i][c];
	}
	for(j=0;j<c+2;j++)
	{
		a[0][j] = a[1][j];
		a[r+1][j] = a[r][j];
	}
	//***********************1st time membership value calculation**************************************//
    float temps=0.0;
	float fmv[p+1][k]; //declaration
	for(i=1;i<=r;i++)
    {
	   for(j=1;j<=c;j++)                          
         {                             //membership function calculation
	        for(d=0;d<k;d++)
	          { 
				fmv[a[i][j]][d] =((float)rand()/(float)((k)*RAND_MAX));
//				temps+=fmv[a[i][j]][d];
//				printf("%.2f ",fmv[a[i][j]][d]);
			  }
//               fmv[a[i][j]][k-1]=temps;          						//This prints exactly between 0.00 & 0.80
		 }
	}
	
	
	//********************random center initialize first time*****************************************//
	
	center=(float *)calloc(k,sizeof(float));
    avg=(float *)calloc(k,sizeof(float));
    
    for(i=0;i<k;i++)
		{	
		    center[i]=rand()%(p+1);             //0-255
            printf("%.2f ",center[i]);
		}
	printf("\n");
	
    float temp6=0.0,temp5=0.0,temp7=0.0;
				
     //************************start of while****************************************//
   g = 0 ;
	while(flag)
	{   
	  for(d=0;d<k;d++)
            {   
			    avg[d]=center[d];
            	center[d]=0.0;
			}		   
			   
	  for(d=0;d<k;d++)
	  {
		 fmvsum=0.0;
			for(i=1;i<=r;i++)
	     	{   
		   		for(j=1;j<=c;j++)
		    	{
		    		for(in=-1;in<2;in++) 
						{
							for(jn=-1;jn<2;jn++) 
								{
								   Xpsum+=(a[i+in][j+jn]);
								}
//								printf("%.2f ",Xpsum);
					    }
                	center[d]+=(float)(pow(fmv[a[i][j]][d],fmp)*(a[i][j]+(float)((alpha/Sr)*Xpsum)));
		   			fmvsum+=(float)(pow(fmv[a[i][j]][d],fmp));
              							//fmvsum is the summation (uij)^m in calculation of center[d]
              		Xpsum=0.0;
                }
                Xpsum=0.0;
            }
             center[d]=center[d]/((1.0+alpha)*(fmvsum));
      }

		printf("\n\nCenters are: ");
			
		for(d=0;d<k;d++)	
			printf("%.2f ",center[d]);

		//************************Distance calculation from the center*************************************// 

	    for(i=1;i<=r;i++)
           { Xpsum=0.0;
              for(j=1;j<=c;j++)
              { 
                 Xpsum=0.0;
                 for(d=0;d<k;d++)
                  {
					dist[i][j][d]=(float)(fabs(a[i][j]-center[d]));
                    
                    for(in=-1;in<2;in++)
                    	{ 
						   for(jn=-1;jn<2;jn++) 
                    	     {
							   scar+=pow((a[i+in][j+jn]-center[d]),2);
							 }
//							 printf("%.2f ",scar);
	                    }
					temp[i][j][d]=(float)scar;
					scar=0.0;
					Xpsum=0.0;
//	                printf("%.2f ",temp[i][j][d]);
                  }
              }
           }
           
        //***************************update membership values**********************************//
        float temp2 = 0.0;                               
        int x;
        for(i=1;i<=r;i++)
        {
	      for(j=1;j<=c;j++)
	        {
	          for(d=0;d<k;d++)
	              {  
	                for(x=0;x<k;x++)
	                    {
						  temp2+=((pow(dist[i][j][d],2)+((float)(alpha/Sr)*temp[i][j][d]))/(pow(dist[i][j][x],2)+((float)(alpha/Sr)*temp[i][j][x])));
	                    }
					fmv[a[i][j]][d]=pow(temp2,-1);
//					printf("%.2f ",fmv[a[i][j]][d]);     
					temp2 = 0.0;
				  }
	        }
	    }       
		flag=0;
	   
		     for(d=0;d<k;d++)
			   {	
			      if(fabs(avg[d]-center[d])>e)
				  {   
					flag=1;
					break;
				  }
			
			   }
	}
	
	//***********************************Partition Coefficient Calculation(Vpc)**************************//
	float Vpc=0.0,Vpe=0.0;
	for(d=0;d<k;d++) for(i=1;i<=r;i++) for(j=1;j<=c;j++) Vpc+=(float)pow(fmv[a[i][j]][d],2);
	
	Vpc=(float)((Vpc)/(r*c));
	printf("\n\nThe Partition Coefficient for BCFCM(Vpc): %.2f ",Vpc);
	
	//***********************************Partition Entropy Calculation(Vpe)**************************//
	for(d=0;d<k;d++) for(i=1;i<=r;i++) for(j=1;j<=c;j++) Vpe+=(float)(fmv[a[i][j]][d]*log(fmv[a[i][j]][d]));
	       	  	
	Vpe=-((Vpe)/(r*c));
	printf("\n\nThe Partition Entropy for BCFCM(Vpe): %.2f ",Vpe);
	
	 //***************************assigning the values to the output matrix**********************************// 
   
   for(i=1;i<=r;i++)
			{	for(j=1;j<=c;j++)
				{	mini=fabs(a[i][j]-center[0]);
					pos=0;
					for(m=0;m<k;m++)
					{
						if(fabs(a[i][j]-center[m])<mini)
						{
							mini=fabs(a[i][j]-center[m]);
							pos=m;
						}
					}
					b[i][j]=(int)center[pos];
				}
			}
			
	for(i=1;i<=r;i++)	
	{
		for(j=1;j<=c;j++)	//to write the output image
		{
			fprintf(fp2," %d",b[i][j]);
        }
        fprintf(fp2,"\n");
	}
	fclose(fp1);fclose(fp2);
	return 0;
}
