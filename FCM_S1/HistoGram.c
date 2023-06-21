#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

int a[2000][2000];
int b[2000][2000];

int main()
{
	int i, j, r, c, m;
	int k = 0;
	char str[10];
	FILE *fp1 = fopen("m2_2_Gray_S1.txt", "r");
	FILE *fp2 = fopen("m2_2_Gray_S1_hist.pgm", "wb");

	fscanf(fp1, "%s %d %d %d", str, &c, &r, &m);
	fprintf(fp2, "%s %d %d %d\n", str, c, r, m);

	for (i = 0; i < r; i++)
	{
		for (j = 0; j < c; j++) // to read the input image
		{
			fscanf(fp1, "%d", &a[i][j]);
		}
	}

	int pix = r * c;

	int hist[m + 1];
	for (i = 0; i < m + 1; i++)
	{
		hist[i] = 0;
	}
	int cdf[m + 1];
	int eh[m + 1];
	for (i = 0; i < r; i++)
	{
		for (j = 0; j < c; j++)
		{
			hist[a[i][j]]++;
		}
	}

	for (i = 0; i < m + 1; i++)
	{
		printf("%d ", hist[i]);
	}
	// To check the number of pixels of each intensity

	// Cumulative freqency distribution

	int cdfmax = m, cdfmin;
	int count = 0;
	for (i = 0; i < m + 1; i++)
	{

		count += hist[i];
		cdf[i] = count;
	}

	printf("\n\n\n");

	// equalized histogram

	for (i = 0; i < m + 1; i++)
	{
		eh[i] = (m * (cdf[i])) / (pix);
		printf("%d ", eh[i]);
	}
	// final image

	for (i = 0; i < r; i++)
	{
		for (j = 0; j < c; j++)
		{
			b[i][j] = eh[a[i][j]];
		}
	}

	for (i = 0; i < r; i++)
	{
		for (j = 0; j < c; j++) // to write the output image
		{
			fprintf(fp2, " %d", b[i][j]);
		}
		fprintf(fp2, "\n");
	}
	fclose(fp1);
	fclose(fp2);
	return 0;
}
