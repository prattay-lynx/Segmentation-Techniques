#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<limits.h>
#define fmp 2
#define iter 0.5

int MINI_OBJ=INT_MAX;
int k;
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
	int i,j,n,d,r,c,p,**a=NULL,**b=NULL;srand(time(0));char str[10];
	int m,**group=NULL,*l=NULL,pos,f=0; float *avg=NULL,*avg1=NULL,*center=NULL;
	float mini,fmvsum=0.0;int flag = 1;float fitness_obj=0.0;				
	int q=0,best_iter=0;
	float ***dist=NULL;
	float temp=0.0;
	int x;
	
	
	FILE *fp1 = fopen("m7_gray.txt","r");
	FILE *fp2 = fopen("m7_grayoriginal_lena_mfcmc.pgm","wb");
	

	fscanf(fp1,"%s %d %d %d",str,&r,&c,&p);
	fprintf(fp2,"%s %d %d %d\n",str,r,c,p);
	
	printf("\n\nEnter the value of k: ");
	scanf("%d",&k);
	
	a=(int **)calloc(r,sizeof(int *));
	for(i=0;i<r;i++)
	{
		a[i]=(int *)calloc(c,sizeof(int));
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
	
	//****************************Declarations************************************************************//
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
    avg1=(float *)calloc(k,sizeof(float));
	group=(int **)calloc(k,sizeof(int *));
	 for(i=0;i<k;i++)	
		group[i]=(int *)calloc(r*c,sizeof(int));
	
	int prev_i=0,prev_j=0;
	for(d=0;d<k;d++)
    {
//     srand(time(0));
	   i=(int)rand()%(p-0+1);
	   j=(int)rand()%(p-0+1);
	   if(i>r) i=r-1;
	   if(j>c) j=c-1;
		center[d]=a[i][j];
	}
		
		sort(center,k);
			
	for(i=0;i<k;i++)
	{
		printf("%0.2f ",center[i]);
	}
		
	//******************Grouping all the initial population together according to center cluster!!****************//	
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
			printf("\nSet: %d and Count : l[%d]: %d\n",m+1,m,l[m]);
			for(i=0;i<l[m];i++)	
			{	
				avg[m]+=group[m][i];                  	//l[m] = number of elements in a cluster j
//				printf("%d ",group[m][i]);      
			}
			avg[m]/=l[m];
		}
		
		//*************************************Objective function Minimization(fitness value)***************************//
		fitness_obj=0.0;
		for(d=0;d<k;d++)
		{
			for(i=0;i<r;i++)
			{
				for(j=0;j<c;j++)
				{
					fitness_obj+=(float)fabs(fabs(pow(fmv[a[i][j]][d],fmp)*a[i][j])-(avg[d]*avg[d]));
//					fitness_obj+=(float)fabs(pow((a[i][j]-avg[d]),2));
//					printf("%.2f ",fitness_obj);
				}
			}
		}
		fitness_obj/=(float)(r*c);
		printf("%0.2f \n ",fitness_obj);
		
		if(fitness_obj<MINI_OBJ)
			{
				MINI_OBJ = (int)fitness_obj;
//				best_iter=f;
				for(d=0;d<k;d++)
					avg1[d]=avg[d];             //assigning the cluster center with minimum objective function!!
		    }
			
		for(i=0;i<k;i++) printf("%0.2f ",avg[i]);

	while(flag)
		{
			for(i=0;i<k;i++)
			{
				center[i]=avg[i];
				avg[i]=0.0;
				l[i]=0;
			}
				
//			sort(center,k);
			
			printf("\nIteration Completed: %d\n",f++);
			printf("\nCenters are: ");
			for(i=0;i<k;i++)	
			{
				printf("%0.2f ",center[i]);
			}
	
	//****************Grouping all the initial population together according to center cluster!!****************//		
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
	    
	    //*************************************Objective function Minimization(fitness value)***************************//
		fitness_obj=0.0;
		for(d=0;d<k;d++)
		{
			for(i=0;i<r;i++)
			{
				for(j=0;j<c;j++)
				{
					fitness_obj+=(float)fabs(fabs(pow(fmv[a[i][j]][d],fmp)*a[i][j])-(avg[d]*avg[d]));
//					fitness_obj+=(float)fabs(pow((a[i][j]-avg[d]),2));
//					printf("%.2f ",fitness_obj);
				}
			}
		}
		fitness_obj/=(float)(r*c);
		printf("%0.2f after iteration %d\n ",fitness_obj,f);
		
		if(fitness_obj<MINI_OBJ)
			{
				MINI_OBJ = (int)fitness_obj;
				best_iter=f;
				for(d=0;d<k;d++)
					avg1[d]=avg[d];             //assigning the cluster center with minimum objective function!!
		    }
	   
	   //**************************************Breaking condition*******************************************************//		
			flag=0;
			for(i=0;i<k;i++)
			{	if(fabs(center[i]-avg[i])>iter)
				{
					flag=1;
					break;
				}
			}
			
		}
	for(d=0;d<k;d++)	
	   center[d]=avg1[d];
	printf("\t*****************************************************\t");
	printf("\nThe best fitness value is: %d obtained in the iteration: %d\n ",MINI_OBJ,best_iter);
	printf("\t*****************************************************\t");
	//***********************************Partition Coefficient Calculation(Vpc)**************************//
	float Vpc=0.0,Vpe=0.0;
	for(d=0;d<k;d++) for(i=0;i<r;i++) for(j=0;j<c;j++) Vpc+=(float)pow(fmv[a[i][j]][d],2);
	
	Vpc=(float)((Vpc)/(r*c));
	printf("\n\nThe Partition Coefficient for MFCMC(Vpc): %.2f ",Vpc);
	
	//***********************************Partition Entropy Calculation(Vpe)**************************//
	for(d=0;d<k;d++) for(i=0;i<r;i++) for(j=0;j<c;j++) Vpe+=(float)(fmv[a[i][j]][d]*log(fmv[a[i][j]][d]));
	       	  	
	Vpe=-((Vpe)/(r*c));
	printf("\n\nThe Partition Entropy for MFCMC(Vpe): %.2f ",Vpe);
	
	//***************************assigning the values to the output matrix**********************************// 

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
				}
			}
		
    for(i=0;i<r;i++)	
	{
		for(j=0;j<c;j++)	//to write the output image
		{
			fprintf(fp2," %d",b[i][j]);
        }
        fprintf(fp2,"\n");
	}
	fclose(fp1);fclose(fp2);
				
	return 0;
}
