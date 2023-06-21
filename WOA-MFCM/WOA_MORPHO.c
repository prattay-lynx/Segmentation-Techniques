#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<limits.h>
#define fmp 2
#define iter 0.05
#define bdash 1
#define pi 3.145
#define C1 6.41
#define C2 59.2

int f_max = INT_MIN;
int f_min = INT_MAX;
int Fbest=INT_MAX;
int best_whalepop=INT_MAX;
float temp8=0.0;

int r,c;
int k;
int in=2;
float At,CV;double P=0.0;int ldash=0;
 void vector_update(int f,int A,int R,int ldash)
 {  srand(time(0));
 	int A2=0;
    //***********************************Updating vectors At,CV,P****************************************************//
	while(in>=0)
	{   
		A=in; 						//linearly selected and decreased from 2 to 0
		in--;
		break;
	}
	if(in<0) in=2;
//	A = 2*(1-f/fmax);
//	A2 = -1+f*((-1)/fmax);
//	ldash = (A2-1)*rand()%(1)+1;
	ldash = rand()%(2);
	printf("\nA:%d ",A);
	printf("\nldash:%d ",ldash);
	R=rand()%(2);
	At = 2*A*R-A;
	CV = 2*A;
	printf("\nCV: %.2f ",CV);
	P = ((double)rand())/((double)RAND_MAX)/2.0 + 0.5 ;
	printf("\nP: %.2lf ",P);
 }
//**************************************Sorting the Array to satisfy the criteria k1<k2<k3<k4....********************//
void sort(float center[k],int n)         
{ 
	int temp,j,i;
	for(i=1;i<n;i++)
	{
		temp=center[i]; 					   //Insertion Sort
		j=i-1;
		while(j>=0 && center[j]>temp)
		{
			center[j+1]=center[j];
			j--;
		}
		center[j+1]=temp;
	}
}

int main()
{
	//*************************************Declarations**************************************************//
	int i,j,n,d,p,**a=NULL,**b=NULL,**binaryimg=NULL;srand(time(0));char str[10];
	int m,**group=NULL,*l=NULL,pos,f=0; float *avg=NULL,*center=NULL;
	float mini,fmvsum=0.0;int flag = 1;	float *fitness_obj=NULL;		
	int q=0; int *whale_pop=NULL,*D=NULL,*D1=NULL,*Xnew=NULL; 
	float ***dist=NULL;int dist_loc;
	float temp=0.0;int *hist = NULL;
	int x,A,R;int bestiter;
	
	
	FILE *fp1 = fopen("original_lena.txt","r");
	FILE *fp4 = fopen("original_lena_morpho.pgm","wb");
	FILE *fp2 = fopen("before_lena_morpho.pbm","wb");
	FILE *fp3 = fopen("after_lena_morpho.pbm","wb");

	

	fscanf(fp1,"%s %d %d %d",str,&c,&r,&p);
	fprintf(fp4,"%s %d %d %d",str,c,r,p);
	fprintf(fp2,"%s %d %d %d\n","P1",c,r,p);
	fprintf(fp3,"%s %d %d %d\n","P1",c,r,p);
	
	printf("\n\nEnter the value of k: ");
	scanf("%d",&k);
	
	a=(int **)calloc(r,sizeof(int *));
	for(i=0;i<r;i++)
	{
		a[i]=(int *)calloc(c,sizeof(int));
	}
	
	binaryimg=(int **)calloc(r,sizeof(int *));
	for(i=0;i<r;i++)
	{
		binaryimg[i]=(int *)calloc(c,sizeof(int));
	}
	
	b=(int **)calloc(r,sizeof(int *));
	for(i=0;i<r;i++)
	{
		b[i]=(int *)calloc(c,sizeof(int));
	}
	
	dist = (float***)calloc(r,sizeof(float**));
	for(i=0;i<r;i++)
	{
		dist[i]=(float**)calloc(c,sizeof(float*));
	}
	for(i=0;i<r;i++)
	{
		for(j=0;j<c;j++)
		{
			dist[i][j]=(float*)calloc(k,sizeof(float));
		}
	}
	
	for(i=0;i<r;i++)	
	{
		for(j=0;j<c;j++)	//to read the input image
		{
			fscanf(fp1,"%d",&a[i][j]);
        }
	}
	

   float fmv[p+1][k+1];
   //***********************1st time Membership Value calculation***************************************//	
   fmvsum=0.0;
   for(i=0;i<r;i++)
    	{
		for(j=0;j<c;j++)                            
            {
			for(d=0;d<k;d++)
	          {
		        fmv[a[i][j]][d] =((float)rand()/(float)((k)*RAND_MAX));			//membership function calculation
//              	printf("%.2f ",fmv[a[i][j]][d]);
			  }
			}
		}
	
    center=(float *)calloc(k,sizeof(float));
    avg=(float *)calloc(k,sizeof(float));
	group=(int **)calloc(k,sizeof(int *));
	 for(i=0;i<k;i++)	
		group[i]=(int *)calloc(r*c,sizeof(int));
	
	int prev_i=0,prev_j=0;
	for(d=0;d<k;d++)
    {
//     srand(time(0));
	   i=(int)rand()%(p+1);
	   j=(int)rand()%(p+1);
	   if(i>r) i=r-1;
	   if(j>c) j=c-1;                						//randomly selected centers form the dataset
		center[d]=a[i][j];
	}
		
		sort(center,k);
			
	for(i=0;i<k;i++)
	{
		printf("%0.2f ",center[i]);
	}
		
	//************************************Grouping all the initial population together according to center cluster!!*********************//	
		l=(int *)calloc(k,sizeof(int));
		for(i=0;i<r;i++)
		{	for(j=0;j<c;j++)
			{	mini=fabs(a[i][j]-center[0]);	
			    pos=0;
				for(m=0;m<k;m++)
				{
					if(mini>fabs(a[i][j]-center[m]))
					{
						mini=fabs(a[i][j]-center[m]);
						pos=m;
					}
				}
				group[pos][l[pos]]=a[i][j];
				++l[pos];
			}
		}
		
		for(m=0;m<k;m++)
		{
			printf("\nSet: %d and Count: l[%d]: %d\n",m+1,m,l[m]);
			for(i=0;i<l[m];i++)	
			{	
				avg[m]+=group[m][i];                  	//l[m] = number of elements in a cluster j
//				printf("%d ",group[m][i]);4       
			}
			avg[m]/=l[m];
		}
		
		
		m=rand()%4;
		printf("\nm: %d\n",m);
		whale_pop=(int*)calloc(k,sizeof(int));
		D=(int*)calloc(k,sizeof(int));
		D1=(int*)calloc(k,sizeof(int));
		Xnew=(int*)calloc(k,sizeof(int));
		int Xrand;

		for(i=0;i<k;i++)
		{
			whale_pop[i]=(int)avg[i];                 //random whale population generated
			if(whale_pop[i]<best_whalepop)
				best_whalepop=whale_pop[i];
		}
		printf("\nWhale population: ");
		for(i=0;i<k;i++)
		{
			printf("%d ",whale_pop[i]);
		}
		
	    int Xbest;
//		Xbest=best_whalepop;
//		printf("\nXbest: %d\n",Xbest);
			
//		for(i=0;i<k;i++) printf("%0.2f ",avg[i]);
		
	//*************************************Objective function Minimization(fitness value)***************************//
		fitness_obj=(float*)calloc(k,sizeof(float));
		for(d=0;d<k;d++)
		{
			for(i=0;i<r;i++)
			{
				for(j=0;j<c;j++)
				{
					fitness_obj[d]+=(float)fabs(pow(fmv[a[i][j]][d],fmp)*a[i][j]-pow(avg[d],2));	
//					printf("%.2f ",fitness_obj[i][d]);
				}
			}
//			fitness_obj[d]/=(float)(r*c);
//			printf("%.2f ",fitness_obj[d]);
			if(fitness_obj[d]<Fbest)
				{ 
					Fbest =(int)fitness_obj[d];
				  	Xbest = (int)avg[d];
				}
				           
		}
		printf("\nXbest: %d\n",Xbest);


	//****************************************Start of while*******************************************************//
	while(flag)
		{
			for(i=0;i<k;i++)
			{
				center[i]=avg[i];avg[i]=0.0;
				l[i]=0;
			}
				
//			sort(center,k);

			vector_update(f,A,R,ldash);                     //Updating the vector quantities
			
			if(P<0.5)
			{
				if(fabs(At)<1)
				{
				   for(j=0;j<k;j++)
				   {
					   D[j] = (int)fabs(CV*Xbest-whale_pop[j]);
					   printf("\ndistloc: %d ",D[j]);            			//shrinking mechanism
					   Xnew[j] = fabs(whale_pop[j]-At*D[j]);
					   printf("\nXnew: %d ",Xnew[j]); 
				   }
				}
				else
				{
					for(j=0;j<k;j++)
					{
						i=(int)rand()%(4);
						Xrand=center[i];
						printf("\nXrand: %d ",D[j]);							    	//spiral mechanism
						D[j] = (int)fabs(CV*Xrand-whale_pop[j]);
						printf("\ndistloc: %d ",D[j]); 
						Xnew[j]=(int)abs(Xrand-At*D[j]);
						printf("\nXnew: %d ",Xnew[j]);
					}
				}
			}
			else
			{
				for(j=0;j<k;j++)
					{
						D1[j] = (int)abs(Xbest-whale_pop[j]);
//						ldash=rand()%(2);
						printf("\ndistloc: %d ",D1[j]); 
						Xnew[j] = abs(D1[j]*exp(bdash*ldash)*cos(2*pi*ldash)+Xbest);
						printf("\nXnew: %d ",Xnew[j]); 
					}
			}
			for(j=0;j<k;j++)
			{
				if(Xnew[j]>center[j])
					whale_pop[j]=center[j];
				else
				{
//					if(Xnew[j]>p) Xnew[j]=p;
					whale_pop[j]=Xnew[j];
				}
				if(whale_pop[j]<best_whalepop)
					best_whalepop=whale_pop[j];
			}
//			Xbest=best_whalepop;
//			printf("\nXbest: %d",abs(Xbest));
			printf("\n");
			
		printf("\nWhale population: ");
		for(i=0;i<k;i++)
		{
			printf("%d ",whale_pop[i]);
		}
			printf("\nIteration Completed: %d\n",f++);
			printf("\nCenters are: ");
			for(i=0;i<k;i++)	
			{
				printf("%0.2f ",center[i]);
			}
	
	//*********************Grouping all the initial population together according to center cluster!!****************//		
			for(i=0;i<r;i++)
			{	for(j=0;j<c;j++)
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
				group[pos][l[pos]]=a[i][j];
				++l[pos];
				}
			}
			    			
			for(m=0;m<k;m++)
			{
				printf("\nSet: %d and Count: l[%d]: %d\n",m+1,m,l[m]);
				for(i=0;i<l[m];i++)	
				{	
					avg[m]+=group[m][i];
//					printf("%d ",group[m][i]);
				}
				avg[m]/=l[m];
			}
	//************************distance calculation from the center*************************************// 
			
	       for(i=0;i<r;i++)
              {
			    for(j=0;j<c;j++)
			    {
                  for(d=0;d<k;d++)
                 {
	                dist[i][j][d]=(float)(fabs(a[i][j]-avg[d]));
				 }
				}
			  }
			
	//***********************************update membership values*************************************//  
        for(i=0;i<r;i++)
        {
	      for(j=0;j<c;j++)
	        {
	          for(d=0;d<k;d++)
	              {  temp = 0.0;
	                for(x=0;x<k;x++)
	                    {
						  temp+=(pow(dist[i][j][d],2)/pow(dist[i][j][x],2));
	                    }
					fmv[a[i][j]][d]=(float)pow(temp,-1);  
//					printf("%f ",fmv[a[i][j]][d]);   
				  }
	        }
	    }
			
		//*************************************Updating (fitness value)***************************//
		for(d=0;d<k;d++)
		{
			fitness_obj[d]=0.0;
		}
		for(d=0;d<k;d++)
		{
			for(i=0;i<r;i++)
			{
				for(j=0;j<c;j++)
				{
					fitness_obj[d]+=(float)fabs(pow(fmv[a[i][j]][d],fmp)*a[i][j]-pow(avg[d],2));	
				}
			}
			fitness_obj[d]/=(float)(r*c);
//			printf("%.2f ",fitness_obj[d]);
			if(fitness_obj[d]<Fbest)
				{
					Fbest =(int)fitness_obj[d]; 
					bestiter=f; 
					Xbest=(int)whale_pop[d];							//fitness of the next iterated population 
				}           
		}
		
		printf("\nXbest: %d",abs(Xbest));
		
	   //**************************************Breaking condition***************************************//
			flag=0;
			for(i=0;i<k;i++)
			{	if(fabs(center[i]-avg[i])>iter)
				{
					flag=1;
					break;
				}
			}
		}
		printf("\nThe best fitness value is: %d and obtained in iteration: %d\n ",Fbest,bestiter);
		
		for(i=0;i<k;i++)
		{
			center[i]=whale_pop[i];
		}
	//**************************Frequency calculation of each grey level vector in the image*****************//
	hist=(int*)calloc(p+1,sizeof(int));
	int count=0;
//	printf("\nHistogram levels of each pixel:\n");
	for(i=0;i<r;i++) for(j=0;j<c;j++) hist[a[i][j]]++;
	
	for(i=0;i<r;i++)
	{
		for(j=0;j<c;j++)
		{
			if(a[i][j]>f_max) f_max =  a[i][j];
			if(a[i][j]<f_min) f_min =  a[i][j];
		}
	}
	printf("\nfmax: %d and fmin: %d\n",f_max,f_min);
	for(i=0;i<p+1;i++)
	{  
	  if(hist[i]!=0) count++;
//	  printf("%d ",hist[i]);
    }
//    printf("\nCount: %d",count);
	
	//***************************assigning the values to the output matrix**********************************// 
   int *count1=NULL;
   count1=(int*)calloc(k,sizeof(int));
   for(i=0;i<r;i++)
			{	for(j=0;j<c;j++)
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
					count1[pos]++;
				}
			}
		//************************************image after segmentation***********************************//
    for(i=0;i<r;i++)	
	{
		for(j=0;j<c;j++)	//to write the segmented image
		{
			fprintf(fp4," %d",b[i][j]);
        }
        fprintf(fp4,"\n");
	}		
		
	//**************************************Uniformity calculation**************************************//
		float Uni=0.0;
	
		int ls = b[r-1][c-1];
			for(i=0;i<r;i++)
			{
				for(j=0;j<c;j++)
				{
				  if(b[i][j]!=ls)
				  {
				  	ls = b[i][j];
				  }
				  if(b[i][j]==b[0][0])
				    temp8++;
				}
			}
//		U=fabs((1-2*k*temp8)/((r*c)*pow((f_min),2)));
        Uni=(float)(temp8)/(float)(r*c);                                     //must lie between 0 & 1
		printf("\nUniformity of the segmented image: %.2f ",Uni);
		
	//**************************************Structural Similarity Index**************************************//
	
	float a_mean=0.0;
	float b_mean=0.0;
	
	for(i=0;i<r;i++)
	{
		for(j=0;j<c;j++)
		{
			a_mean+=a[i][j];
			b_mean+=b[i][j];
		}
	}
	a_mean/=(float)(r*c);
	b_mean/=(float)(r*c);
	
	//************* standard deviation ****************//
    float a_std = 0, b_std = 0, img12_cov = 0;
    for(i = 0; i < r; i++)
	{
		for(j = 0; j < c; j++)
	  {
	        a_std += pow((a[i][j] - a_mean),2);
	        b_std += pow((b[i][j] - b_mean),2);
	        img12_cov += (a[i][j] - a_mean) * (b[i][j] - b_mean);
      }
    }
	a_std = sqrt(a_std/(float)(r*c));
	b_std = sqrt(b_std/(float)(r*c));
	 /************* ssim ****************/
	double ssim = (2*a_mean*b_mean + C1)*(2*img12_cov + C2);
	/************* deno ****************/
    double deno = (a_mean*a_mean + b_mean*b_mean - C1) *(a_std*a_std + b_std*b_std - C2);
    ssim=(double)ssim / (double)deno;
    printf("\nssim: %0.2f ",ssim);
    
    
	//***********************************Partition Coefficient Calculation(Vpc)**************************//
	printf("\n----------------------------------------------------------------------------------------------");
	float Vpc=0.0,Vpe=0.0;
	for(d=0;d<k;d++) for(i=0;i<r;i++) for(j=0;j<c;j++) Vpc+=(float)pow(fmv[a[i][j]][d],2);
	
	Vpc=(float)((Vpc)/(r*c));
	printf("\n\nThe Partition Coefficient for WOA-MFCM(Vpc): %.2f ",Vpc);
	
	//***********************************Partition Entropy Calculation(Vpe)**************************//
	printf("\n----------------------------------------------------------------------------------------------");
	for(d=0;d<k;d++) for(i=0;i<r;i++) for(j=0;j<c;j++) Vpe+=(float)(fmv[a[i][j]][d]*log(fmv[a[i][j]][d]));
	       	  	
	Vpe=-((Vpe)/(r*c));
	printf("\n\nThe Partition Entropy for WOA-MFCM(Vpe): %.2f ",Vpe);	
    //***********************************Conversion of the segmented image to binary image********************************************//
    float sum1=0.0;
    for(i=0;i<r;i++)	
	{
		for(j=0;j<c;j++)
		{
			sum1=sum1+b[i][j];
		}
	}
	float threshold=0.0;
	threshold=sum1/(r*c);            //Calculation of threshold value
	
	for(i=0;i<r;i++)	
	{
		for(j=0;j<c;j++)	
		{
			if(b[i][j] >= threshold) binaryimg[i][j]=1;
			else binaryimg[i][j]=0;
		}
	}
	//***********************************************image before morpho****************************************************************//
    for(i=0;i<r;i++)	
	{
		for(j=0;j<c;j++)	//to write the output image
		{
			fprintf(fp2," %d",binaryimg[i][j]);
        }
        fprintf(fp2,"\n");
	}
	
	//**********************************************Erosion******************************************************//
	int se[3][3] = {{1,1,1},
	{1,1,1},
	{1,1,1}};
	
	int **binimg=NULL,in,jn;
	
	binimg=(int **)calloc(r,sizeof(int *));
	for(i=0;i<r;i++)
	{
		binimg[i]=(int *)calloc(c,sizeof(int));
	}
    
    for(i=0;i<r;i++)
	{
		for(j=0;j<c;j++)
		{
			binimg[i][j]=binaryimg[i][j];
		}
	}
	
	
	for(i=0;i<r;i++)
	{
		for(j=0;j<c;j++)
		{
			for(in=-1;in<=1;in++)
			{
				for(jn=-1;jn<=1;jn++)
				{
					if(se[i+in][j+jn])
					{
						if((i+in>=0) && (i+in < r) && (j+jn>= 0) && (j+jn< c))
						{
							if(binimg[i+in][j+jn]<binaryimg[i][j]) binaryimg[i][j]=binimg[i+in][j+jn];
						}
					}
				}
			}
		}
	}
	
	//**********************************************Dilation******************************************************//
	for(i=0;i<r;i++)
	{
		for(j=0;j<c;j++)
		{
			for(in=-1;in<=1;in++)
			{
				for(jn=-1;jn<=1;jn++)
				{
					if(se[i+in][j+jn])
					{
						if((i+in>=0) && (i+in < r) && (j+jn>= 0) && (j+jn< c))
						{
							if(binimg[i+in][j+jn]>binaryimg[i][j]) binaryimg[i][j]=binimg[i+in][j+jn];
						}
					}
				}
			}
		}
	}
	//********************************************image after morpho*******************************************************************//
    for(i=0;i<r;i++)	
	{
		for(j=0;j<c;j++)	//to write the output image
		{
			fprintf(fp3," %d",binaryimg[i][j]);
        }
        fprintf(fp3,"\n");
	}
	fclose(fp1);fclose(fp2);fclose(fp3);fclose(fp4);
				
	return 0;
}
