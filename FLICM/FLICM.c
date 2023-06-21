#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#define fmp 2
#define Sr 9.0
#define alpha 2.0
#define e 0.05
//********************FLICM was introduced by Krinidis and does not need to manually insert  parameters*******//
int main()
{
	
	int i,j,k,n,d,r,c,p;srand(time(0));char str[10];
	int m,*l=NULL,pos,f=0; float *avg=NULL,*center=NULL;int fm;
	float fmvsum=0.0,sum =0.0,mini;	int flag = 1;int win[10]; float Xpsum=0.0;int q;int in,jn;
    int lx; int count = 0;float Sig,Sij;int h;float scar=0.0;
	float sigma=0.0;

    int s=0;
	int **a=NULL,**b=NULL;
	float ***dist=NULL;
	float **g=NULL;
	float d_value=0.0;
	
	FILE *fp1 = fopen("original_lena.txt","r");
	FILE *fp2 = fopen("original_lena_elifcm.txt","w");
	
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
  
//  	g = (float***)calloc(r,sizeof(float**));
//	for(i=1;i<=r;i++)
//	{g[i]=(float**)calloc(c,sizeof(float*));
//		for(j=1;j<=c;j++)
//		{
//			g[i][j]=(float*)calloc(k,sizeof(float));
//		}
//	}
	
	g=(float **)calloc(r+2,sizeof(float *));
	for(i=0;i<r+2;i++)
	{
		g[i]=(float *)calloc(c+2,sizeof(float));
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
                	center[d]+=(float)(pow(fmv[a[i][j]][d],fmp)*(a[i][j]));
		   			fmvsum+=(float)(pow(fmv[a[i][j]][d],fmp));
              							//fmvsum is the summation (uij)^m in calculation of center[d]
                }
              
            }
             center[d]=center[d]/(fmvsum);
      }

		printf("\n\nCenters are: ");
			
		for(d=0;d<k;d++)	
			printf("%.2f ",center[d]);
		printf("\n");

		//************************Distance calculation from the center*************************************// 
	    //******************************Objective Function Calculation************************************//
	    int x1,y1,jn1,in1;
	    for(i=1;i<=r;i++)
           { Xpsum=0.0; 
              for(j=1;j<=c;j++)
              { 
                 Xpsum=0.0;scar =0.0;
                 for(d=0;d<k;d++)
                  {	 scar =0.0;
                 for(in1=-1;in1<2;in1++)
		           		{    
		           	  		for(jn1=-1;jn1<2;jn1++)  // R X R --> Search Window and S X S --> Neighbourhood Window
		           	  			{ 	 
                                   for(in=-1;in<2;in++)
		           					{    
						   				for(jn=-1;jn<2;jn++) 
		           	  					{ 
											x1=i+in+in1;
											y1=j+jn+jn1;
											if(x1<0) x1=0;
											if(x1>r) x1=r;
											if(y1<0) y1=0;
											if(y1>c) y1=c;
//											printf("%d ",a[i+in+in1][j+jn+jn1]-a[i+in][j+jn]);
		           	  					   d_value=(float)abs(a[x1][y1]-a[i+in][j+jn]);
//		           	  					   printf("%.2f ",d_value);
							               scar+=((1.0/((d_value)+1.0))*pow((1.0-fmv[a[i][j]][d]),fmp)*pow((a[i+in][j+jn]-center[d]),2));							              
									    }								
									}
								}
						}
                 
					dist[i][j][d]=fabs(a[i][j]-center[d]);
					
                  }
                  g[i][j]=scar;
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
						  temp2+=(pow(dist[i][j][d],2)+((float)g[i][j])/(pow(dist[i][j][x],2)+((float)g[i][j])));
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
	printf("\n\nThe Partition Coefficient for FCM_S1(Vpc): %.2f ",Vpc);
	
	//***********************************Partition Entropy Calculation(Vpe)**************************//
	for(d=0;d<k;d++) for(i=1;i<=r;i++) for(j=1;j<=c;j++) Vpe+=(float)(fmv[a[i][j]][d]*log(fmv[a[i][j]][d]));
	       	  	
	Vpe=-((Vpe)/(r*c));
	printf("\n\nThe Partition Entropy for FCM_S1(Vpe): %.2f ",Vpe);
	
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
