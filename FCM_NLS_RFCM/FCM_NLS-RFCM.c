#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define fmp 2
#define threshold 0.01
#define Wi 0.9
#define delta 0.3
#define Sr 9.0
#define B 1.0
#define h 99.0 //[should lie between 5,8.....100]
#define var 0.03  //[should lie between 0.002 ..... 0.03] --> Gaussian Variance

int main()
{
	int i, j, k, n, d, r, c, p, **a = NULL,**ai=NULL, **b = NULL;
	int in =0,jn=0;
	int in1=0,jn1=0;
	float len=0.0;
    int si=0,hi=0,x1,y1;
	srand(time(0));
	char str[10];
	int m, pos, f = 0;
	float *avg = NULL, *center = NULL;
	int fm;
	//	int lowerApp_flag = 0, boundary_flag = 0;
	float fmvsum = 0.0, mini;
	int flag = 1;
	int x;
	float *sum = NULL;
	float ***dist=NULL,***dist1=NULL,**Z=NULL;
	float temp_fmv[100];
	int q = 0;
	int l;
	float sigma = 0.0,temp5=0.0;
	float **w=NULL;
	

	FILE *fp1 = fopen("m12_gray.txt", "r");
	FILE *fp2 = fopen("m12_wocvanls_4.pgm", "w");

	fscanf(fp1, "%s %d %d %d", str, &c, &r, &p);
	fprintf(fp2, "%s %d %d %d\n", str, c, r, p);

	printf("\n\nEnter the value of k: ");
	scanf("%d", &k); // Take no. of clusters
	
	ai=(int **)calloc(r+2,sizeof(int *));
	for(i=0;i<r+2;i++)
	{
		ai[i]=(int *)calloc(c+2,sizeof(int));
	}
	
	a = (int **)calloc(r + 2, sizeof(int *));
	for (i = 0; i < r + 2; i++)
	{
		a[i] = (int *)calloc(c + 2, sizeof(int));
	}

	b = (int **)calloc(r + 2, sizeof(int *));
	for (i = 0; i < r + 2; i++)
	{
		b[i] = (int *)calloc(c + 2, sizeof(int));
	}

	Z=(float **)calloc(r,sizeof(float *));
	for(i=1;i<=r;i++)
	{
		Z[i]=(float *)calloc(c,sizeof(float));
	}

	dist = (float ***)calloc(r, sizeof(float **));
	for (i = 1; i <= r; i++)
	{
		dist[i] = (float **)calloc(c, sizeof(float *));
		for (j = 1; j <= c; j++)
		{
			dist[i][j] = (float *)calloc(k, sizeof(float));
		}
	}
	
	dist1 = (float***)calloc(r,sizeof(float**));
	for(i=1;i<=r;i++)
	{
		dist1[i]=(float**)calloc(c,sizeof(float*));
		for(j=1;j<=c;j++)
		{
			dist1[i][j]=(float*)calloc(k,sizeof(float));
		}
	}
	
	w=(float **)calloc(r,sizeof(float *));
	for(i=1;i<=r;i++)
	{
		w[i]=(float *)calloc(c,sizeof(float));
	}
	

	/**********************************************************************************************************/
	for (i = 1; i <= r; i++)
	{
		for (j = 1; j <= c; j++) // to read the input image
		{
			fscanf(fp1, "%d", &a[i][j]);
		}
	}

	//*************************matrix padding************************************//

	for (i = 0; i < r + 2; i++)
	{
		a[i][0] = a[i][1];
		a[i][c + 1] = a[i][c];
	}
	for (j = 0; j < c + 2; j++)
	{
		a[0][j] = a[1][j];
		a[r + 1][j] = a[r][j];
	}
	//**************************************************************************//

	float temp = 0.0;
	float temps = 0.0;
	float temp1 = 0.0;
	float temp33=0.0;
	float fmv[p + 1][k + 1];
	float Xpsum = 0.0;
	
	//********************random center initialize first time*****************************************//

	center = (float *)calloc(k, sizeof(float));
	avg = (float *)calloc(k, sizeof(float));
	sum = (float *)calloc(k, sizeof(float));

	//************************************************************************************************//

	for (i = 0; i < k; i++)
	{
		avg[i] = ((double)rand()/RAND_MAX) * (p + 1); // 0-255
		printf("%.2f ", avg[i]);
	}
	printf("\n");

	//************************Calculation Of The Non local Spatial Information for Each Pixel**************************************//
	for(i=1;i<=r;i++)
	     	{     
		   		for(j=1;j<=c;j++)
		    	{  temp33=0.0;
		    	
		    	        	for(in=-1;in<2;in++)
		           					{    
						   				for(jn=-1;jn<2;jn++) 
		           	  					{
           	  		      					temp5+=(float)(fabs(a[i+in][j+jn]-a[i][j]));
		                				}
		                			}
		        	sigma=(float)sqrt((1.0+temp5)/(float)Sr);
		        	temp5=0.0;
//		        	printf("%.2f ",sigma);
		        	
					for(in1=-1;in1<2;in1++)
		           		{    
		           	  		for(jn1=-1;jn1<2;jn1++)  // R X R --> Search Window and S X S --> Neighbourhood Window
		           	  			{ 	 
								len=0.0; 
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
		           	  					   
							               len+=(float)pow((a[x1][y1]-a[i+in][j+jn]),2);
//							               hi++;
//							               printf("%.2f ",len);
									    }
//									    si++;
									}
									Z[i][j]=(float)exp(-(float)(var*len)/(pow(h,2)));	
									temp33+=Z[i][j];
								}                  			// Normalized Factor: Z , h is the decay Parameter
						}                                  //a[i][j]-a[i+in][j+jn] 
				   	   	
//				printf("%.2f ",Z);
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
										w[i][j]=(float)(1.0/(float)(temp33))*(float)exp(-(var*pow((a[x1][y1]-a[i+in][j+jn]),2))/(pow(h,2)));  
//										count+=w[i][j];     						//Wij must lie between 0 and 1
						    		}
						    		
				   					}
//				   					printf("%.2f ",w[i][j]); 
				   				}
				   			}  
				}
				
			}
			
				for(i=1;i<=r;i++) 
					   	 {
							for(j=1;j<=c;j++)    
					    	  {
							  	for(in=-1;in<2;in++)
			           			{
			           	  			for(jn=-1;jn<2;jn++)
			           	  			{ 
							   			ai[i][j]+=(float)((float)w[i][j]*(float)a[i+in][j+jn]);	  //Non local spatial information for each pixel
						            }
						            	if(ai[i][j]>p) ai[i][j]=p;
										if(ai[i][j]<0) ai[i][j]=0;
				                }
				            }
				        }
	int lowerApp_flag = 0, boundary_flag = 0;
	//*********************start of while****************************************//
	int first_index, second_index;
	while (flag)
	{
		for (d = 0; d < k; d++)
		{
			center[d] = avg[d];
			//center[d]=0.0;
		}
		for (i = 1; i <= r; i++)
		{
			for (j = 1; j <= c; j++)
			{
				for (d = 0; d < k; d++)
				{
					dist[i][j][d] = (float)(fabs(a[i][j] - center[d]));
					dist1[i][j][d]=(float)(fabs(ai[i][j]-center[d]));
//					printf("%.2f ",dist[i][j][d]);
				}
			}
		}

		for (i = 1; i <= r; i++)
		{
			for (j = 1; j <= c; j++)
			{
				for (d = 0; d < k; d++)
				{
					temps = 0.0;
					for (x = 0; x < k; x++)
					{
						temps+=(float)(pow(dist[i][j][d],2)+(B*pow(dist1[i][j][d],2)))/(float)((pow(dist[i][j][x],2))+(B*pow(dist1[i][j][x],2)));
					}
//						printf("%.2f ",temp1);
					fmv[a[i][j]][d] = (double)pow(temps, -1);
//					printf("%.2f ",fmv[a[i][j]][d]);
				}

				temp1 = 0;
				(fmv[a[i][j]][0] > fmv[a[i][j]][1]) ? (first_index = 0, second_index = 1) : (first_index = 1, second_index = 0);
				for (l = 2; l < k; l++)
				{
					if (fmv[a[i][j]][l] > fmv[a[i][j]][first_index])
					{
						second_index = first_index;
						first_index = l;
					}
					else if (fmv[a[i][j]][l] > fmv[a[i][j]][second_index])
						second_index = l;
				}

				if ((fmv[a[i][j]][first_index] - fmv[a[i][j]][second_index]) > delta)
				{
					for (l = 0; l < k; l++)
						if (l == first_index)
							fmv[a[i][j]][l] = 1.0;
						else
							fmv[a[i][j]][l] = 0.0;
				}
			}
		}

		for (l = 0; l < k; l++)
		{
			boundary_flag = 0;
			for (i = 1; i <= r; i++)
			{
				for (j = 1; j <= c; j++)
				{
					if (fmv[a[i][j]][l] != 0.0 && fmv[a[i][j]][l] != 1.0)
					{
						boundary_flag = 1; // boundary set
						break;
					}
				}
				if (boundary_flag)
					break;
			}

			if (boundary_flag)
			{
				lowerApp_flag = 0;
				for (i = 1; i <=r; i++)
				{
					for (j = 1; j <= c; j++)
					{
						if (fmv[a[i][j]][l] == 1.0)
						{
							lowerApp_flag = 1; // cluster set
							break;
						}
					}
					if (lowerApp_flag)
						break;
				}
			}
			else
				lowerApp_flag = 1;

			if (boundary_flag && lowerApp_flag) // mixed
			{
				sum[l] = 0.0;
				avg[l] = 0.0;
				for (i = 1; i <= r; i++)
				{
					for (j = 1; j <= c; j++)
					{
						if (fmv[a[i][j]][l] == 1.0)
						{
							sum[l]++;
							avg[l] += a[i][j];
						}
					}
				}
				temps = avg[l] / sum[l];
				sum[l] = 0.0;
				avg[l] = 0.0;
				for (i = 1; i <= r; i++)
				{
					for (j = 1; j <= c; j++)
					{
						if (fmv[a[i][j]][l] != 0.0 && fmv[a[i][j]][l] != 1.0)
						{
							temp = pow(fmv[a[i][j]][l], fmp);
							sum[l] += temp*(1+B);
							avg[l] += temp * ((float)a[i][j]+(float)(B*ai[i][j]));
						}
					}
				}
				avg[l] = Wi * temps + (1 - Wi) * (avg[l] / sum[l]);
			}
			else
			{
				sum[l] = 0.0;
				avg[l] = 0.0;
				for (i = 1; i <= r; i++)
				{
					for (j = 1; j <= c; j++)
					{
						temp = pow(fmv[a[i][j]][l], fmp);
						sum[l] += temp*(1+B);
						avg[l] += temp * ((float)a[i][j]+(float)(B*ai[i][j]));
					}
				}
				avg[l] /= sum[l];
			}
		}

		printf("\n\nCenters : ");
		for (i = 0; i < k; i++)
			printf("%0.2f ", avg[i]);

		flag = 0;

		for (d = 0; d < k; d++)
		{
			if (fabs(avg[d] - center[d]) > 0.2)
			{
				flag = 1;
				break;
			}
		}
	}

	//***********************************Partition Coefficient Calculation(Vpc)**************************//
		float vpc=0.0,vpe=0.0;
		for(i=0;i<r;i++)
			{	
				for(j=0;j<c;j++)
				{	
					for(l=0;l<k;l++)
					{	
						if(fmv[a[i][j]][l])
						{
							vpc+=pow(fmv[a[i][j]][l],2);
							vpe+=(fmv[a[i][j]][l]*log(fmv[a[i][j]][l]));
						}
					}
				}
			}
			vpc = vpc/(r*c);
			vpe = -vpe/(r*c);
			printf("\n\nPartition coefficient : %0.2f ",vpc);
			printf("\nPartition entropy : %0.2f ",vpe);
	
	//***************************assigning the values to the output matrix**********************************//

	for (i = 1; i <= r; i++)
	{
		for (j = 1; j <= c; j++)
		{
			mini = fabs(a[i][j] - center[0]);
			pos = 0;
			for (m = 0; m < k; m++)
			{
				if (fabs(a[i][j] - center[m]) < mini)
				{
					mini = fabs(a[i][j] - center[m]);
					pos = m;
				}
			}
			b[i][j] = (int)center[pos];
		}
	}

	for (i = 1; i <= r; i++)
	{
		for (j = 1; j <= c; j++) // to write the output image
		{
			fprintf(fp2, " %d", b[i][j]);
		}
		fprintf(fp2, "\n");
	}
	fclose(fp1);
	fclose(fp2);
	/*************************************************Memoery deallocation***************************************/
	free(a);
	free(b);
	free(ai);
	free(w);
	free(Z);
	free(avg);
	free(center);
	return 0;
}
