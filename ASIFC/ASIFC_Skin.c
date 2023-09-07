#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#define fmp 2
#define Sr 9
#define max(a, b) (((a > b)) ? (a) : (b))
// #define k 4

int main()
{

	int i, j, k, n, d, r, c, p;
	srand(time(0));
	char str[10];
	int m, pos, f = 0;
	float *avg = NULL, *center = NULL;
	int fm;
	float fmvsum = 0.0, sum = 0.0, mini;
	int flag = 1;
	int win[10];
	float Xpsum = 0.0;
	int q;
	int in, jn;
	int lx;
	int count = 0;
	float mu = 0.0, tao = 0.0;
	float **eta = NULL, **lamdaS = NULL, **lamdaG = NULL;
	int h;
	float maxi = -999999.999;
	float **sigma = NULL;
	float temp5 = 0.0, deno = 0.0;
	;
	float **C = NULL;
	int x;
	float arr[5];
	//	int *hist = NULL;
	float **E = NULL;
	int s = 0;
	int g;
	int **a = NULL, **b = NULL, **a1 = NULL, **tempimg = NULL, **binaryimg = NULL;
	float ***dist = NULL;
	float scar = 0.0;
	float ***temp = NULL;

	FILE *fp1 = fopen("m7_gray.txt", "r");
	FILE *fp3 = fopen("m7_7_gray.txt", "r");
	FILE *fp2 = fopen("m7_gray_ASFIC.pgm", "wb");
	FILE *fp4 = fopen("m7_gray_binary.pbm", "wb");
	FILE *fp5 = fopen("m7_gray_cva.pgm","wb");

	fscanf(fp1, "%s %d %d %d", str, &c, &r, &p);
	fscanf(fp3, "%s %d %d %d", str, &c, &r, &p);
	fprintf(fp2, "%s %d %d %d\n", str, c, r, p);
	fprintf(fp4, "P1 %d %d %d", c, r, p);
	fprintf(fp5, "%s %d %d %d", str,c, r, p);

	printf("\n\nEnter the value of k: ");
	scanf("%d", &k); // Input no. of clusters

	a = (int **)calloc(r + 2, sizeof(int *));
	for (i = 0; i < r + 2; i++)
	{
		a[i] = (int *)calloc(c + 2, sizeof(int));
	}

	a1 = (int **)calloc(r + 2, sizeof(int *));
	for (i = 0; i < r + 2; i++)
	{
		a1[i] = (int *)calloc(c + 2, sizeof(int));
	}

	tempimg = (int **)calloc(r + 2, sizeof(int *));
	for (i = 0; i < r + 2; i++)
	{
		tempimg[i] = (int *)calloc(c + 2, sizeof(int));
	}

	binaryimg = (int **)calloc(r + 2, sizeof(int *));
	for (i = 0; i < r + 2; i++)
	{
		binaryimg[i] = (int *)calloc(c + 2, sizeof(int));
	}

	b = (int **)calloc(r + 2, sizeof(int *));
	for (i = 0; i < r + 2; i++)
	{
		b[i] = (int *)calloc(c + 2, sizeof(int));
	}

	E = (float **)calloc(r, sizeof(float *));
	for (i = 1; i <= r; i++)
	{
		E[i] = (float *)calloc(c, sizeof(float));
	}

	C = (float **)calloc(r, sizeof(float *));
	for (i = 1; i <= r; i++)
		C[i] = (float *)calloc(c, sizeof(float));

	lamdaS = (float **)calloc(r, sizeof(float *));
	for (i = 1; i <= r; i++)
		lamdaS[i] = (float *)calloc(c, sizeof(float));

	lamdaG = (float **)calloc(r, sizeof(float *));
	for (i = 1; i <= r; i++)
		lamdaG[i] = (float *)calloc(c, sizeof(float));
	sigma = (float **)calloc(r, sizeof(float *));
	for (i = 1; i <= r; i++)
		sigma[i] = (float *)calloc(c, sizeof(float));

	eta = (float **)calloc(r, sizeof(float *));
	for (i = 1; i <= r; i++)
		eta[i] = (float *)calloc(c, sizeof(float));

	dist = (float ***)calloc(r, sizeof(float **));
	for (i = 1; i <= r; i++)
	{
		dist[i] = (float **)calloc(c, sizeof(float *));
		for (j = 1; j <= c; j++)
		{
			dist[i][j] = (float *)calloc(k, sizeof(float));
		}
	}

	temp = (float ***)calloc(r, sizeof(float **));
	for (i = 1; i <= r; i++)
	{
		temp[i] = (float **)calloc(c, sizeof(float *));
		for (j = 1; j <= c; j++)
		{
			temp[i][j] = (float *)calloc(k, sizeof(float));
		}
	}
	//*************************Input matrix of image************************************//

	for (i = 1; i <= r; i++)
		for (j = 1; j <= c; j++)
			fscanf(fp1, "%d", &a[i][j]);
	// to read the input image

	//*************************Input matrix of image1************************************//

	for (i = 1; i <= r; i++)
	{
		for (j = 1; j <= c; j++) // to read the input image
		{
			fscanf(fp3, "%d", &a1[i][j]);
		}
	}

//		for (i = 1; i <= r; i++)
//		{
//			for (j = 1; j <= c; j++) // to read the input image
//			{
//				printf("%d ",a1[i][j]);
//			}
//		}

	//*************************CVA of temp1************************************//

	for (i = 1; i <= r; i++)
	{
		for (j = 1; j <= c; j++)
		{
			tempimg[i][j] = sqrt(pow((a[i][j] - a1[i][j]), 2));
		}
	}

		for (i = 1; i <= r; i++)
		{
			for (j = 1; j <= c; j++) // to write the output image
			{
				fprintf(fp5," %d ",tempimg[i][j]);
			}
			fprintf(fp5,"\n");
		}

	float sum2=0.0;
    for(i=0;i<r;i++)	
	{
		for(j=0;j<c;j++)
		{
			sum2=sum2+tempimg[i][j];
		}
	}
	float threshold_temp=0.0;
	threshold_temp=sum2/(r*c);
	printf("\n%.2f\n ",threshold_temp); 
	
//	for (i = 1; i <= r; i++)
//	{
//		for (j = 1; j <= c; j++)
//		{
//			if(tempimg[i][j]==0)
//			   tempimg[i][j]=(int)(threshold_temp);
//		}
//	}
	//*************************matrix padding************************************************//

	for (i = 0; i < r + 2; i++)
	{
		tempimg[i][0] = tempimg[i][1];
		tempimg[i][c + 1] = tempimg[i][c];
	}
	for (j = 0; j < c + 2; j++)
	{
		tempimg[0][j] = tempimg[1][j];
		tempimg[r + 1][j] = tempimg[r][j];
	}

	//***********************1st time membership value calculation**************************************//
	float temp_value = 0.0;
	float fmv[p + 1][k]; // declaration
	for (i = 1; i <= r; i++)
		for (j = 1; j <= c; j++)
		{ // membership function calculation
			for (d = 0; d < k; d++)
			{
				fmv[tempimg[i][j]][d] = ((float)rand() / (float)((k)*RAND_MAX));
				//				temp_value+=fmv[tempimg[i][j]][d];
				//				printf("%.2f ",fmv[tempimg[i][j]][d]);
			}
			//			  fmv[tempimg[i][j]][k-1]=1.0-temp_value;
		}
	//********************random center initialize first time*****************************************//

	center = (float *)calloc(k, sizeof(float));
	avg = (float *)calloc(k, sizeof(float));

	center[0] = 50.0;
	center[1] = 100.0;
	center[2] = 150.0;
	center[3] = 170.0;
	for (d = 0; d < k; d++)
	{
		i = (int)rand() % (p + 1);
		j = (int)rand() % (p + 1);
		if (i > r)
			i = r - 1;
		if (j > c)
			j = c - 1; // randomly selected centers form the dataset
		center[d] = tempimg[i][j];
		printf("%.2f ", center[d]);
	}

	//****************************************calculation of lamda smoothing factor****************************************************//

	for (i = 1; i <= r; i++)
	{
		for (j = 1; j <= c; j++)
		{
			E[i][j] = 1 + ((float)1 / (float)(r * c));
			//			   printf("%.2f ",E[i][j]);
		}
	}

	float t1 = 0.0, t2 = 0.0, t3 = 0.0;
	for (i = 1; i <= r; i++)
	{
		for (j = 1; j <= c; j++)
		{
			for (in = -1; in < 2; in++)
			{
				for (jn = -1; jn < 2; jn++)
				{
					t1 += fabs(tempimg[i][j] - tempimg[i + in][j + jn]);
				}
			}
			mu = ((float)1 / (float)Sr) * t1;
			t1 = 0.0;

			for (in = -1; in < 2; in++)
				for (jn = -1; jn < 2; jn++)
					t2 += pow(((tempimg[i][j] - tempimg[i + in][j + jn]) - mu), 2);
			eta[i][j] = (float)(0.5 + pow((((float)1 / (float)Sr) * t2), 0.5));
			t2 = 0.0;
			//           		printf("%.2f ",eta[i][j]);
		}
	}

	for (i = 1; i <= r; i++)
	{
		for (j = 1; j <= c; j++)
		{
			for (in = -1; in < 2; in++)
			{
				for (jn = -1; jn < 2; jn++)
				{
					t3 += fabs(tempimg[i + in][j + jn]);
				}
			}
			tao = ((float)1 / (float)Sr) * t3;
			//           	               	    printf("%f ",tao);
			t3 = 0.0;

			for (in = -1; in < 2; in++)
			{
				for (jn = -1; jn < 2; jn++)
				{
					sigma[i][j] += pow((((float)1 / (float)Sr) * tao), 0.5);
					//						  printf("%.2f ",sigma[i][j]);
				}
			}

			lamdaS[i][j] = (float)((eta[i][j]) / (1.0 + sigma[i][j]));
			//           	   	printf("%.2f ",lamdaS[i][j]);
			if (lamdaS[i][j] > maxi)
				maxi = lamdaS[i][j];
		}
	}
	printf("\n%.2f", (float)maxi);
	for (i = 1; i <= r; i++)
	{
		for (j = 1; j <= c; j++)
		{
			lamdaG[i][j] = (float)(0.1 + (lamdaS[i][j] / maxi));
			//		    	printf("%.2f ",lamdaG[i][j]);
		}
	}

	float temp3 = 0.0, temp4 = 0.0;

	//************************start of while****************************************//
	while (flag)
	{
		for (d = 0; d < k; d++)
		{
			avg[d] = center[d];
			center[d] = 0.0;
		}
		Xpsum = 0.0;

		for (d = 0; d < k; d++)
		{
			fmvsum = 0.0;
			for (i = 1; i <= r; i++)
			{
				for (j = 1; j <= c; j++)
				{

					for (in = -1; in < 2; in++)
					{
						for (jn = -1; jn < 2; jn++)
						{
							Xpsum += 0.5 + (tempimg[i + in][j + jn]);
						}
					}
					//           		printf("%.2f ",Xpsum);

					center[d] += (float)(pow(fmv[tempimg[i][j]][d], fmp) * ((float)E[i][j]) * (lamdaG[i][j] * tempimg[i][j] + (((float)fabs((1.0 - lamdaG[i][j])) / (float)Sr) * Xpsum)));
					fmvsum += (float)(pow(fmv[tempimg[i][j]][d], fmp) * E[i][j]);
					Xpsum = 0.0;
				}
				Xpsum = 0.0;
			}
			center[d] = (center[d] / (fmvsum));
			//	        printf("%.2f ",center[d]);
		}

		printf("\n\nCenters are: ");

		for (d = 0; d < k; d++)
		{
			printf("%.2f ", center[d]);
		}

		temp3 = 0.0, temp4 = 0.0;
		deno = 0.0;
		for (d = 0; d < k; d++)
		{
			for (i = 1; i <= r; i++)
			{
				for (j = 1; j <= c; j++)
				{
					deno += (float)(fmv[tempimg[i][j]][d] * E[i][j]);
				}
			}
			//	     		printf("%.2f ",deno);
			arr[d] = deno;
			deno = 0.0;
		}
		//	 	deno=0.0;
		//	 	printf("%.2f ",temp4);

		for (i = 1; i <= r; i++)
		{
			for (j = 1; j <= c; j++)
			{
				for (d = 0; d < k; d++)
				{
					temp3 += (fmv[tempimg[i][j]][d] * (-1) * log((float)fmv[tempimg[i][j]][d] / (float)arr[d]));
				}

				//	        	printf("%.2f ",temp3);
				C[i][j] = (float)(exp(temp3));
				//				printf("%.2f ",C[i][j]);
				temp3 = 0.0;
				temp5 += C[i][j];
			}
		}

		//	 	for(i=1;i<=r;i++)
		//        {
		//	      for(j=1;j<=c;j++)
		//	        {
		//	           temp5+=C[i][j];
		//	        }
		//	    }
		//	    printf("%.2f ",temp5);

		for (i = 1; i <= r; i++)
		{
			for (j = 1; j <= c; j++)
			{
				E[i][j] = (float)(1 + ((E[i][j] * C[i][j]) / (float)temp5));
				//	        	printf("%.2f ",E[i][j]);
			}
		}
		temp5 = 0.0;
		//************************Distance calculation from the center*************************************//

		for (i = 1; i <= r; i++)
			for (j = 1; j <= c; j++)
				for (d = 0; d < k; d++)
					dist[i][j][d] = (float)(fabs(tempimg[i][j] - center[d]));

		//***************************update membership values**********************************//

		float temp2 = 0.0;

		for (i = 1; i <= r; i++)
		{
			for (j = 1; j <= c; j++)
			{
				for (d = 0; d < k; d++)
				{
					for (x = 0; x < k; x++)
					{
						temp2 += (pow(dist[i][j][d], 2)) / (pow(dist[i][j][x], 2));
					}
					fmv[tempimg[i][j]][d] = pow(temp2, -1);
					//					printf("%f ",fmv[a[i][j]][d]);
					temp2 = 0.0;
				}
			}
		}

		flag = 0;

		for (d = 0; d < k; d++)
		{
			if (fabs(avg[d] - center[d]) > 0.5)
			{
				flag = 1;
				break;
			}
		}
	}

	//***********************************Partition Coefficient Calculation(Vpc)**************************//
	float Vpc = 0.0, Vpe = 0.0;
	for (d = 0; d < k; d++)
		for (i = 1; i <= r; i++)
			for (j = 1; j <= c; j++)
				Vpc += (float)pow(fmv[tempimg[i][j]][d], 2);

	Vpc = (float)((Vpc) / (r * c));
	printf("\n\nThe Partition Coefficient for ASIFC(Vpc): %.2f ", Vpc);

	//***********************************Partition Entropy Calculation(Vpe)**************************//
	for (d = 0; d < k; d++)
		for (i = 1; i <= r; i++)
			for (j = 1; j <= c; j++)
				Vpe += (float)(fmv[tempimg[i][j]][d] * log(fmv[tempimg[i][j]][d]));

	Vpe = -((Vpe) / (r * c));
	printf("\n\nThe Partition Entropy for ASIFC(Vpe): %.2f ", Vpe);

	//***************************assigning the values to the output matrix**********************************//

	for (i = 1; i <= r; i++)
	{
		for (j = 1; j <= c; j++)
		{
			mini = fabs(tempimg[i][j] - center[0]);
			pos = 0;
			for (m = 0; m < k; m++)
			{
				if (fabs(tempimg[i][j] - center[m]) < mini)
				{
					mini = fabs(tempimg[i][j] - center[m]);
					pos = m;
				}
			}
			if (center[pos] < 0)
				b[i][j] = (int)((-1) * center[pos]);
			else
				b[i][j] = (int)(center[pos]);
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

	//*******************************************************//
	float sum1 = 0.0;
	for (i = 0; i < r; i++)
	{
		for (j = 0; j < c; j++)
		{
			sum1 = sum1 + tempimg[i][j];
		}
	}
	float threshold = 0.0;
	threshold = sum1 / (r * c);
	printf("\n%.2f ", threshold); // Calculation of threshold value

	for (i = 1; i <= r; i++)
	{
		for (j = 1; j <= c; j++)
		{
			if (b[i][j] >= threshold)
				binaryimg[i][j] = 0;
			else
				binaryimg[i][j] = 1;
		}
	}
	//***********************************************image****************************************************************//
	for (i = 1; i <= r; i++)
	{
		for (j = 1; j <= c; j++) // to write the output image
		{
			fprintf(fp4, " %d", binaryimg[i][j]);
		}
		fprintf(fp4, "\n");
	}

	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);
	fclose(fp5);
	free(a);
	free(a1);
	return 0;
}
