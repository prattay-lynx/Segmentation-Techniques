#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#define k 4
#define fpm 2
#define e 0.05
#define d 0.8
#define w 0.8

int main()
{
	int i,j,r,c,m,l,p,**a=NULL,index,count=0,flag,b_flag,c_flag,first_index,second_index;	
	char str[25];	srand(time(0));
	float *avg=NULL,*sum=NULL,*cen=NULL,***mem=NULL,min,temp,temps,vpc,vpe;
	
	FILE *fp1 = fopen("m7_gray.txt","rb");
	FILE **fp2 = (FILE **)calloc(k+1,sizeof(FILE *));
	for(i=0;i<k+1;i++)	
	{
		sprintf(str,"m7_gray_new_%d.pgm",i);
		fp2[i]=fopen(str,"w");
	}
	
	fscanf(fp1,"%s %d %d %d",str,&c,&r,&m);
	for(i=0;i<k+1;i++)	fprintf(fp2[i],"%s %d %d %d",str,c,r,m);
	
	a=(int **)calloc(r,sizeof(int *));
	for(i=0;i<r;i++)	a[i]=(int *)calloc(c,sizeof(int));
	
	for(i=0;i<r;i++)	
		for(j=0;j<c;j++)	
			fscanf(fp1,"%d",&a[i][j]);
				
	avg=(float *)calloc(k,sizeof(float));
	sum=(float *)calloc(k,sizeof(float));
	cen=(float *)calloc(k,sizeof(float));
	
	mem=(float ***)calloc(r,sizeof(float **));
	{
		for(i=0;i<r;i++)
		{
			mem[i]=(float **)calloc(c,sizeof(float *));
			for(j=0;j<c;j++)
				mem[i][j]=(float *)calloc(k,sizeof(float));
		}	
	}
	
	for(i=0;i<k;i++)	avg[i]=(float)rand()/RAND_MAX*(m+1);
	
	do{
		for(i=0;i<k;i++)
			cen[i]=avg[i];
			
		printf("\n\n\nIteration : %d",++count);
		printf("\n\nCen : ");	for(i=0;i<k;i++)	printf("%0.4f ",cen[i]);
	
		for(i=0;i<r;i++)
		{	
			for(j=0;j<c;j++)
			{	
				for(l=0;l<k;l++)
				{	
					temps=0.0;
					for(p=0;p<k;p++)
						temps+=pow(a[i][j]-cen[l],2)/pow(a[i][j]-cen[p],2);
						
					mem[i][j][l]=1.0/(pow(temps,1/(fpm-1)));
				}
				
				
				(mem[i][j][0]>mem[i][j][1])?(first_index=0,second_index=1):(first_index=1,second_index=0);
				for(l=2;l<k;l++)
				{
					if(mem[i][j][l]>mem[i][j][first_index])
					{
						second_index=first_index;
						first_index=l;
					}
					else if(mem[i][j][l]>mem[i][j][second_index])
						second_index=l;
				}
				
				if((mem[i][j][first_index]-mem[i][j][second_index])>d)
				{
					for(l=0;l<k;l++)
						if(l==first_index)
							mem[i][j][l]=1.0;
						else
							mem[i][j][l]=0.0;
				}
			}
		}
		
		for(l=0;l<k;l++)
		{	
			b_flag=0;		
			for(i=0;i<r;i++)
			{
				for(j=0;j<c;j++)
				{	
					if(mem[i][j][l]!=0.0 && mem[i][j][l]!=1.0)
					{
						b_flag=1;		//boundary set
						break;
					}
				}
				if(b_flag)	break;
			}
			
			if(b_flag)
			{
				c_flag=0;
				for(i=0;i<r;i++)
				{
					for(j=0;j<c;j++)
					{	
						if(mem[i][j][l]==1.0)
						{
							c_flag=1;		//cluster set
							break;
						}
					}
					if(c_flag)	break;
				}
			}
			else	c_flag=1;
			
			if(b_flag && c_flag)	//mixed		
			{
				sum[l]=0.0;avg[l]=0.0;
				for(i=0;i<r;i++)
				{
					for(j=0;j<c;j++)
					{
						if(mem[i][j][l]==1.0)
						{
							sum[l]++;
							avg[l]+=a[i][j];
						}
					}
				}
				temps=avg[l]/sum[l];
				sum[l]=0.0;avg[l]=0.0;
				for(i=0;i<r;i++)
				{
					for(j=0;j<c;j++)
					{
						if(mem[i][j][l]!=0.0 && mem[i][j][l]!=1.0)
						{
							temp=pow(mem[i][j][l],fpm);
							sum[l]+=temp;
							avg[l]+=temp*a[i][j];
						}
					}
				}
				avg[l]=w*temps+(1-w)*(avg[l]/sum[l]);
			}
			else 
			{
				sum[l]=0.0;avg[l]=0.0;
				for(i=0;i<r;i++)
				{
					for(j=0;j<c;j++)
					{
						temp=pow(mem[i][j][l],fpm);
						sum[l]+=temp;
						avg[l]+=temp*a[i][j];
					}
				}
				avg[l]/=sum[l];
			}
		}
			
		printf("\n\nAvg : ");	for(i=0;i<k;i++)	printf("%0.4f ",avg[i]);
		
		flag=0;
		for(i=0;i<k;i++)
		{	if(fabs(cen[i]-avg[i])>e)
			{
				flag=1;break;
			}
		}
		if(flag==0)
		{
			for(i=0;i<r;i++)
			{	
				for(j=0;j<c;j++)
				{	
					min=fabs(a[i][j]-cen[0]);
					index=0;
					for(l=1;l<k;l++)
					{
						if(fabs(a[i][j]-cen[l])<min)
						{
							min=fabs(a[i][j]-cen[l]);
							index=l;
						}
					}
					for(l=1;l<k+1;l++)
					{
						if(l-1==index)	fprintf(fp2[l]," %d",(int)cen[index]);
						else			fprintf(fp2[l]," %d",0);
					}
					fprintf(fp2[0]," %d",(int)cen[index]);
				}
			}
			vpc=0.0;vpe=0.0;
			for(i=0;i<r;i++)
			{	
				for(j=0;j<c;j++)
				{	
					for(l=0;l<k;l++)
					{	
						if(mem[i][j][l])
						{
							vpc+=pow(mem[i][j][l],2);
							vpe+=(mem[i][j][l]*log(mem[i][j][l]));
						}
					}
				}
			}
			vpc = vpc/(r*c);
			vpe = -vpe/(r*c);
			printf("\n\nPartition coefficient : %0.6f ",vpc);
			printf("\nPartition entropy : %0.6f ",vpe);
		}
		
	}while(flag);
	
	for(i=0;i<r;i++)	free(a[i]);
	free(a);
	free(avg);free(sum);free(cen);
	
	for(i=0;i<r;i++)
	{
		for(j=0;j<c;j++)
			free(mem[i][j]);

		free(mem[i]);
	}
	free(mem);
	
	fclose(fp1);
	for(i=0;i<k+1;i++)	fclose(fp2[i]);
	return 0;
}
