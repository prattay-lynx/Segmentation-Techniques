#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#define fmp 2
#define Sr 9.0
#define B 1.0
#define h 100.0 //[should lie between 5,8.....100]
#define var 0.03  //[should lie between 0.002 ..... 0.03] --> Gaussian Variance
#define iter 0.5  //for iteration

//**************************************No. of clusters will be taken as an input during runtime*********************//
int main()
{
	int i,j,k,n,d,r,c,p;srand(time(0));char str[10];
	int m,*l=NULL,pos,f=0; float *avg=NULL,*center=NULL;int fm;
	float fmvsum=0.0,sum =0.0,mini;	int flag = 1;int win[10]; float Xpsum=0.0;int q;int in,jn,in1,jn1;
	int lx; float count = 0.0;
    int x;float sigma = 0.0,temp5=0.0;
    
	int s=0; int g;
	int **a=NULL,**b=NULL;
	float **w=NULL;int **ai=NULL;
	float ***dist=NULL,***dist1=NULL,**Z=NULL;
//	float Z=0.0;
    float len=0.0;
    int si=0,hi=0,x1,y1;
	
	FILE *fp1 = fopen("rubella-german-measles_m1.txt","r");
	FILE *fp2 = fopen("rubella-german-measles_m1_fcm_nls.pgm","wb");
	
	fscanf(fp1,"%s %d %d %d",str,&c,&r,&p);
	fprintf(fp2,"%s %d %d %d\n",str,c,r,p);
	
	printf("\n\nEnter the value of k: ");
	scanf("%d",&k);                           //Input no. of clusters
	
	
	a=(int **)calloc(r+2,sizeof(int *));
	for(i=0;i<r+2;i++)
	{
		a[i]=(int *)calloc(c+2,sizeof(int));
	}
	//****************************************Array initializations******************************//
	ai=(int **)calloc(r+2,sizeof(int *));
	for(i=0;i<r+2;i++)
	{
		ai[i]=(int *)calloc(c+2,sizeof(int));
	}
	
	b=(int **)calloc(r+2,sizeof(int *));
	for(i=0;i<r+2;i++)
	{
		b[i]=(int *)calloc(c+2,sizeof(int));
	}
	
	w=(float **)calloc(r,sizeof(float *));
	for(i=1;i<=r;i++)
	{
		w[i]=(float *)calloc(c,sizeof(float));
	}
	
	Z=(float **)calloc(r,sizeof(float *));
	for(i=1;i<=r;i++)
	{
		Z[i]=(float *)calloc(c,sizeof(float));
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
	
	dist1 = (float***)calloc(r,sizeof(float**));
	for(i=1;i<=r;i++)
	{
		dist1[i]=(float**)calloc(c,sizeof(float*));
		for(j=1;j<=c;j++)
		{
			dist1[i][j]=(float*)calloc(k,sizeof(float));
		}
	}
	
	for(i=1;i<=r;i++) 
	   for(j=1;j<=c;j++) 
	       fscanf(fp1,"%d",&a[i][j]);     //To read the input image
	       
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
	
	//***********************1st time membership value calculation**************************************//
    float temp_value=0.0;
	float fmv[p+1][k]; 								//declaration
	for(i=1;i<=r;i++)
	   for(j=1;j<=c;j++)
	   {                          					//membership function calculation
	        for(d=0;d<k;d++)
	          { 
				fmv[a[i][j]][d] =((float)rand()/(float)((k)*RAND_MAX));
//				temp_value+=fmv[a[i][j]][d];
//				printf("%.2f ",fmv[a[i][j]][d]);
			  }
//			  fmv[a[i][j]][k-1]=1.0-temp_value;
		}
												          //This prints exactly between 0.00 & 0.80
		
	//***********************random center initialize first time*****************************************//
	
	center=(float *)calloc(k,sizeof(float));
    avg=(float *)calloc(k,sizeof(float));
    
    for(i=0;i<k;i++)
		{	
		    center[i]=rand()%(p+1);             //0-255
            printf("%2f ",center[i]);
		}
	printf("\n");

	//************************Calculation Of The Non local Spatial Information for Each Pixel**************************************//
    float temp33=0.0;
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
//		        			        	 printf("%.2f ",sigma);
		        	
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
					   	 for(j=1;j<=c;j++)    
					    	  for(in=-1;in<2;in++)
			           			{
			           	  			for(jn=-1;jn<2;jn++)
			           	  			{ 
							   			ai[i][j]+=(float)((float)w[i][j]*(float)a[i+in][j+jn]);	  //Non local spatial information for each pixel
						            }
						            	if(ai[i][j]>p) ai[i][j]=p;
										if(ai[i][j]<0) ai[i][j]=0;
				                }
//	for(i=1;i<=r;i++)
//	     	{   
//		   	 for(j=1;j<=c;j++)
//		    	{ 
//		    		printf("%d ",ai[i][j]);
//		        }
//		    }
				
     //****************************************start of while****************************************//
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
                	center[d]+=(float)(pow(fmv[a[i][j]][d],fmp)*((float)a[i][j]+(float)(B*ai[i][j])));
		   			fmvsum+=(float)(pow(fmv[a[i][j]][d],fmp));	   
				}
				//fmvsum is the summation (uij)^m  //in calculation of center[d]
            }
             center[d]=(float)(center[d]/((float)(1+B)*fmvsum));
     }

		printf("\n\nCenters are: ");
			
		for(d=0;d<k;d++) printf("%.2f ",center[d]);
//		printf("\n");

		//************************Distance calculation from the center*************************************// 
	    for(i=1;i<=r;i++)
              for(j=1;j<=c;j++)
                 for(d=0;d<k;d++)
                  {
					dist[i][j][d]=(float)(fabs(a[i][j]-center[d]));
					dist1[i][j][d]=(float)(fabs(ai[i][j]-center[d]));
                  }
           
        //***************************update membership values**********************************//
        float temp2 = 0.0;                               
      
        for(i=1;i<=r;i++)
        {
	      for(j=1;j<=c;j++)
	        {
	          for(d=0;d<k;d++)
	              {  
	                for(x=0;x<k;x++)
	                    {
						  temp2+=(float)(pow(dist[i][j][d],2)+(B*pow(dist1[i][j][d],2)))/(float)((pow(dist[i][j][x],2))+(B*pow(dist1[i][j][x],2)));
	                    }
					fmv[a[i][j]][d]=pow(temp2,-1);
//					printf("%f ",fmv[a[i][j]][d]);     
					temp2 = 0.0;
				  }
	        }
	    }       
		flag=0;
	   
		     for(d=0;d<k;d++)
			   {	
			      if(fabs(avg[d]-center[d])>iter)
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
	printf("\n\nThe Partition Coefficient for FCM_NLS(Vpc): %.2f ",Vpc);
	
	//***********************************Partition Entropy Calculation(Vpe)**************************//
	for(d=0;d<k;d++) for(i=1;i<=r;i++) for(j=1;j<=c;j++) Vpe+=(float)(fmv[a[i][j]][d]*log(fmv[a[i][j]][d]));
	       	  	
	Vpe=-((Vpe)/(r*c));
	printf("\n\nThe Partition Entropy for FCM_NLS(Vpe): %.2f ",Vpe);
	
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
