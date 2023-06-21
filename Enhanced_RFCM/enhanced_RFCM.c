#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define fmp 2
#define threshold 0.01
#define w 0.9
#define delta 0.3
#define Sr 9.0
#define B 0.2

int main()
{
	int i, j, k, n, d, r, c, p, **a = NULL, **b = NULL;
	float **E = NULL;
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
	float ***dist = NULL;
	float temp_fmv[100];
	int q = 0;
	int l;
	int *hist = NULL;

	FILE *fp1 = fopen("original_lena.txt", "r");
	FILE *fp2 = fopen("original_lena_enrfcm1.pgm", "w");

	fscanf(fp1, "%s %d %d %d", str, &c, &r, &p);
	fprintf(fp2, "%s %d %d %d\n", str, c, r, p);

	printf("\n\nEnter the value of k: ");
	scanf("%d", &k); // Take no. of clusters

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

	E = (float **)calloc(r, sizeof(float *));
	for (i = 1; i <= r; i++)
	{
		E[i] = (float *)calloc(c, sizeof(float));
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
	hist = (int *)calloc(p + 1, sizeof(int));

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
	float fmv[p + 1][k + 1];
	int in, jn;
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

	for (i = 1; i <= r; i++)
	{
		for (j = 1; j <= c; j++)
		{
			for (in = -1; in < 2; in++)
			{
				for (jn = -1; jn < 2; jn++)
				{
					Xpsum += (a[i + in][j + jn]);
					//printf("%f ",Xpsum);
				}
			}
			E[i][j] = ((float)1 / (float)(1 + B)) * (a[i][j] + ((float)B / (float)Sr) * Xpsum);
			// arr[l]=(float)E[i][j];
			//l++;			printf("%.2f ",E[i][j]);
			Xpsum = 0.0;
		}
		Xpsum = 0.0;
	}
	int lowerApp_flag = 0, boundary_flag = 0;

	for (i = 1; i <= r; i++)
		for (j = 1; j <= c; j++)
			hist[a[i][j]]++;
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
					dist[i][j][d] = (float)(fabs(E[i][j] - center[d]));
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
						temps += (double)(pow(dist[i][j][d], 2.0) / (pow(dist[i][j][x], 2.0)));
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
							sum[l] += temp * ((float)hist[a[i][j]]);
							avg[l] += temp * ((float)hist[a[i][j]]) * ((float)E[i][j]);
						}
					}
				}
				avg[l] = w * temps + (1 - w) * (avg[l] / sum[l]);
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
						sum[l] += temp * ((float)hist[a[i][j]]);
						avg[l] += temp * ((float)hist[a[i][j]]) * ((float)E[i][j]);
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
			if (fabs(avg[d] - center[d]) > 0.5)
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
	free(avg);
	free(E);
	free(center);
	return 0;
}
