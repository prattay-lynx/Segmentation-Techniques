#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#define fmp 2
#define Sr 9
#define B 1
//#define max(a,b)    (((a>b))?(a):(b))
#define lamdaS 3
#define lamdaG 6
float max(int a,int b);
int main()
{
	int i,j,k,n,d,r,c,p;srand(time(0));char str[10];
	int m,*l=NULL,pos,f=0; float *avg=NULL,*center=NULL;int fm;
	float fmvsum=0.0,sum =0.0,mini;	int flag = 1;int win[10]; float Xpsum=0.0;int q;int in,jn;
	int lx; int count = 0;float Sig,Sij;int h;
	float sigma=0.0; float e1=0.0,e2=0.0;
	int *hist = NULL; float ***fmv=NULL;
//	float **E=NULL; 
    int s=0; int g;
	int **a=NULL,**b=NULL;
	float ***dist=NULL;
	float ***temp=NULL;
//    float EE[Sr][Sr];

	
	FILE *fp1 = fopen("pepper.ascii.txt","r");
	FILE *fp2 = fopen("pepper.ascii_fgfcm.pgm","wb");
	
	fscanf(fp1,"%s %d %d %d",str,&c,&r,&p);
	fprintf(fp2,"%s %d %d %d\n",str,c,r,p);
	
	printf("\n\nEnter the value of k: ");
	scanf("%d",&k);                           //Input no. of clusters
	
	
	a=(int **)calloc(r+2,sizeof(int *));
	for(i=0;i<r+2;i++)
	{
		a[i]=(int *)calloc(c+2,sizeof(int));
	}
	
	b=(int **)calloc(r,sizeof(int *));
	for(i=0;i<r;i++)
	{
		b[i]=(int *)calloc(c,sizeof(int));
	}
	
//	int hist[p+1];
	
	dist = (float***)calloc(r,sizeof(float**));
	for(i=0;i<r;i++)
	{
		dist[i]=(float**)calloc(c,sizeof(float*));
		for(j=0;j<c;j++)
		{
			dist[i][j]=(float*)calloc(k,sizeof(float));
		}
	}
	
	fmv = (float***)calloc(r,sizeof(float**));
	for(i=0;i<r;i++)
	{
		fmv[i]=(float**)calloc(c,sizeof(float*));
		for(j=0;j<c;j++)
		{
			fmv[i][j]=(float*)calloc(k,sizeof(float));
		}
	}
  
	hist=(int*)calloc(p+1,sizeof(int));

    //*************************Input matrix of image************************************//	
    
   for(i=1;i<=r;i++)	
	{
		for(j=1;j<=c;j++)	                    //to read the input image
		{
			fscanf(fp1,"%d",&a[i][j]);
		    
        }
	}
	float temps=0.0,temp5=0.0,temp7=0.0;
	//*************************matrix padding************************************//
	
	for(i=0;i<r+2;i++)	
	{
		a[i][0]=a[i][1];
        a[i][c+1]=a[i][c];
	}
	for(j=0;j<c+2;j++)
	{
		a[0][j] = a[1][j] ;
		a[r+1][j] = a[r][j];
	}
	//***********************1st time membership matrix calculation**************************************//

//	float fmv[p+1][k]; //declaration
	for(i=0;i<r;i++)
    {
	   for(j=0;j<c;j++)                          
         {                             //membership function calculation
	        for(d=0;d<k-1;d++)
	          { 
				fmv[i][j][d]=((float)rand()/(float)((k)*RAND_MAX));
				temps+=fmv[i][j][d];
//				printf("%.2f ",fmv[i][j][d]);
			  }
			  fmv[i][j][k-1]=1.0-temps;
               //printf("%.2f ",sum);           //This prints exactly between 0.00 & 0.80
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

    
	//************************assignment of EE[a[i][j]][d]**************************************//
	for(i=1;i<=r;i++)
	     	{   temp5=0.0;
		   		for(j=1;j<=c;j++)
		    	{
	                 e1=0.0;e2=0.0;temp5=0.0;
				     for(in=-1;in<2;in++)
		           			{
		           	  			for(jn=-1;jn<2;jn++) 
		           	  		         temp5+=(float)(fabs(a[i+in][j+jn]-a[i][j]));            //3 X 3 filter window
				        	}
					        	sigma=(float)sqrt((1+temp5)/(float)Sr);
					        	temp5=0.0;
//		        	 printf("\nsigma:\n");
//		        	 printf("%.2f ",sigma);
		        for(in=-1;in<2;in++)
           			{
           	  			for(jn=-1;jn<2;jn++)
           	  			{ 
					            if(in != 0 || jn != 0)
					            {
						  			Sij=(float)-(max((fabs(in)),(fabs(jn))))/(float)(lamdaS);
//						  			printf("%.2f ",Sij);	
					            	Sig=(float)-((fabs(a[i][j]-a[i+in][j+jn]))/((float)(lamdaG)*(float)sigma));
//					            	printf("%.2f ",Sig);
						  			e2+=exp(Sij+Sig);
						  			e1+=(a[i+in][j+jn]*exp(Sij+Sig));				  
							    }
							}   
						}		
					temp5=0.0;
					b[i-1][j-1]=(int)(e1/e2);
					if(b[i-1][j-1]>p) b[i-1][j-1]=p;
					if(b[i-1][j-1]<0) b[i-1][j-1]=0;
					}
				}
				
	//*****************************Frequency of image calculation******************************************//
	int x;
	for(i=0;i<r;i++) for(j=0;j<c;j++) hist[b[i][j]]+=1;    //frequency calculation of each grey level vector in the image
	 
	printf("Histogram levels of each pixel:\n");
	for(i=0;i<p+1;i++)
	{  
	  printf("%d ",hist[i]);
	}
			
     //************************start of while****************************************//
   g = 0 ;
	while(1)
	{   
	  for(d=0;d<k;d++)
            {   
			    avg[d]=center[d];
            	center[d]=0.0;
			}		   
			   
	  for(d=0;d<k;d++)
	  {
		 fmvsum=0.0;
			for(i=0;i<r;i++)
	     	{   
		   		for(j=0;j<c;j++)
		    	{
                	center[d]+=(float)(pow(fmv[i][j][d],fmp)*((float)hist[b[i][j]])*((float)b[i][j]));
		   			fmvsum+=(float)(pow(fmv[i][j][d],fmp)*hist[b[i][j]]);
				}   
				
			}
			 center[d]=(center[d]/(fmvsum));	//fmvsum is the summation (uij)^m in calculation of center[d]
      }	

		printf("\n\nCenters are: ");
			
		for(d=0;d<k;d++)	
		{
			printf("%.2f ",center[d]);
		}
		//************************Distance calculation from the center*************************************// 
	
	    for(i=0;i<r;i++) for(j=0;j<c;j++) for(d=0;d<k;d++) dist[i][j][d]=(float)(fabs(b[i][j]-center[d]));
           
        //***************************update membership values**********************************//
        float temp2 = 0.0;                               
      
        for(i=0;i<r;i++)
        {
	      for(j=0;j<c;j++)
	        {
	          for(d=0;d<k;d++)
	              {  
	                for(x=0;x<k;x++)
	                    {
						  temp2+=(pow(dist[i][j][d],2))/(pow(dist[i][j][x],2));
	                    }
					fmv[i][j][d]=pow(temp2,-1);
//					printf("%f ",fmv[i][j][d]);     
					temp2 = 0.0;
				  }
	        }
	    }       
		flag=0;
	   
	     for(d=0;d<k;d++)
		   {	
		      if(fabs(avg[d]-center[d])>0.5)
			  {   
				flag=1;
				break;
			  }
		
		   }
	if(flag==0)
	{
	 //***************************assigning the values to the output matrix**********************************// 
   
      for(i=0;i<r;i++)
			{	for(j=0;j<c;j++)
				{	mini=fabs(b[i][j]-center[0]);
					pos=0;
					for(m=0;m<k;m++)
					{
						if(fabs(b[i][j]-center[m])<mini)
						{
							mini=fabs(b[i][j]-center[m]);
							pos=m;
						}
					}
					fprintf(fp2," %d",(int)center[pos]);
				}
					fprintf(fp2,"\n");
			}
			break;
		}
	}

    //***********************************Partition Coefficient Calculation(Vpc)**************************//
	float Vpc=0.0,Vpe=0.0;
	for(d=0;d<k;d++) for(i=0;i<r;i++) for(j=0;j<c;j++) Vpc+=(float)(pow(fmv[i][j][d],2));
	
	Vpc=(float)((Vpc)/(r*c));
	printf("\n\nThe Partition Coefficient for fgFCM(Vpc): %.2f ",Vpc);
	
	//***********************************Partition Entropy Calculation(Vpe)**************************//
	for(d=0;d<k;d++) for(i=0;i<r;i++) for(j=0;j<c;j++) Vpe+=(float)(fmv[i][j][d]*log(fmv[i][j][d]));
	       	  	
	Vpe=-((Vpe)/(r*c));
	printf("\n\nThe Partition Entropy for fgFCM(Vpe): %.2f ",Vpe);
	
	fclose(fp1);fclose(fp2);
	return 0;
}
float max(int a,int b)
{
	if(a>b) return a;
	return b;
}
