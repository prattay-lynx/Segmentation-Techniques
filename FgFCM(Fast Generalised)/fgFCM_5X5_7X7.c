#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#define fmp 2
#define Sr 9
#define B 1
#define max(a,b)    (((a>b))?(a):(b))
#define lamdaS 3
#define lamdaG 6


int showmenu()
{   
    int ch;
	printf("\n1.3 X 3\n2.5 X 5\n3.7 X 7\n\nEnter your choice:");
	scanf("%d",&ch);
	return ch;
}

int main()
{
	
	int i,j,k,n,d,r,c,p;srand(time(0));char str[10];
	int m,*l=NULL,pos,f=0; float *avg=NULL,*center=NULL;int fm;
	float fmvsum=0.0,sum =0.0,mini;	int flag = 1;int win[10]; float Xpsum=0.0;int q;int in,jn;
	float sum_avg=0.0;int lx; int count = 0;float Sig,Sij;int h;
	float sigma=0.0;int ch;
	int *hist = NULL;
	float **E=NULL; int s=0; int g;
	int **a=NULL,**b=NULL;
	float ***dist=NULL;float scar = 0.0;
	float ***temp=NULL;
	float **EE=NULL;
//    float EE[Sr][Sr];

	
	FILE *fp1 = fopen("pepper.ascii.txt","r");
	FILE *fp2 = fopen("pepper.ascii_fgfcm.txt","w");
	
	fscanf(fp1,"%s %d %d %d",str,&c,&r,&p);
	fprintf(fp2,"%s %d %d %d\n",str,c,r,p);
	
	printf("\n\nEnter the value of k: ");
	scanf("%d",&k);                           //Input no. of clusters
	
	ch=showmenu();
	
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
	
	E=(float **)calloc(r,sizeof(float *));
    for(i=1;i<=r;i++)
	{
		E[i]=(float *)calloc(c,sizeof(float));
	}
    
	EE=(float **)calloc(r+2,sizeof(float *));
	for(i=0;i<r+2;i++)
	{
		EE[i]=(float *)calloc(c+2,sizeof(float));
	}
	
//	int hist[p+1];
	
	dist = (float***)calloc(r,sizeof(float**));
	for(i=1;i<=r;i++)
	{
		dist[i]=(float**)calloc(c,sizeof(float*));
		for(j=1;j<=c;j++)
		{
			dist[i][j]=(float*)calloc(k,sizeof(float));
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

	float fmv[p+1][k]; //declaration
	for(i=1;i<=r;i++)
    {
	   for(j=1;j<=c;j++)                          
         {                             //membership function calculation
	        for(d=0;d<k;d++)
	          { 
		        
				fmv[a[i][j]][d] =((float)rand()/(float)((k)*RAND_MAX));
				//printf("%.2f ",fmv[a[i][j]][d]);
			  }
               //printf("%.2f ",sum);           //This prints exactly between 0.00 & 0.80
		 }
	}
	
	
	//********************random center initialize first time*****************************************//
	
	center=(float *)calloc(k,sizeof(float));
    avg=(float *)calloc(k,sizeof(float));
    
    for(i=0;i<k;i++)
		{	
		    center[i]=rand()%(p+1);             //0-255
            printf("%2f ",center[i]);
		}
	printf("\n");
	//*****************************Frequency of image calculation******************************************//
	int x;
	
	for(i=1;i<=r;i++)
	{
		for(j=1;j<=c;j++)
		{
			hist[a[i][j]]++;
		}		
	}	
	printf("\nThe histogram levels are:\n");
	for(i=0;i<p+1;i++)
	{  
	  printf("%d ",hist[i]);
	}
	
    float temp6=0.0,temp5=0.0,temp7=0.0;
	//************************assignment of EE[a[i][j]][d]**************************************//

	if(ch==1)
	{
	for(i=1;i<=r;i++)
	     	{   temp5=0.0;
		   		for(j=1;j<=c;j++)
		    	{	      
		       for(in=-1;in<2;in++)
           			{
           	  			for(jn=-1;jn<2;jn++)
           	  			{
           	  		      temp5+=(float)(fabs(a[i+in][j+jn]-a[i][j]));
		                }
		            }
		        	sigma=(float)sqrt((1+temp5)/(float)Sr);
		        	temp5=0.0;
//		        	 printf("\nsigma:\n");
//		        	 printf("%.2f ",sigma);
		  
//	            printf("\nEE:"); 
		        for(in=-1;in<2;in++)
           			{
           	  			for(jn=-1;jn<2;jn++)
           	  			{ 
					            if(in == 0 && jn == 0)
					            {
					                EE[i+in][j+jn]=0.0;
//					                printf("%.2f ",EE[i+in][j+jn]);	
						  		}
						  		else
						  		{
						  			Sij=(float)(fabs(max(in,jn)))/(float)(lamdaS);
//						  			printf("%.2f ",Sij);	
					            	Sig=(float)(fabs(a[i][j]-a[i+in][j+jn]))/((float)(lamdaG)*(float)sigma);
//					            	printf("%.2f ",Sig);
						  			EE[i+in][j+jn]=exp(-Sij-Sig);
//						  			printf("%.2f ",EE[i+in][j+jn]);				  
							    }
//							    printf("%.2f ",EE[i+in][j+jn]);	
							}   
						}		
					temp5=0.0;
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
							temp6+=((float)EE[i+in][j+jn]*(float)a[i+in][j+jn]); 
							temp7+=EE[i+in][j+jn];
			            }
			        }
			    	E[i][j]=(float)temp6/(float)temp7;
//			    	printf("%.2f ",E[i][j]);
			    	temp6=0.0;
			    	temp7=0.0;
			   }
	        }
    } 
    //************************************** 5 X 5 *******************************//
    else if(ch==2)
    {
    	for(i=1;i<=r;i++)
	     	{   temp5=0.0;
		   		for(j=1;j<=c;j++)
		    	{     
//					win[0]=a[i-1][j-2];
//					win[1]=a[i-2][j-1];    //5 * 5  filter window 
//					win[2]=a[i-2][j+1];
//					win[3]=a[i-1][j+2];
//					win[4]=a[i+1][j+2];
//					win[5]=a[i+2][j+1];
//					win[6]=a[i+2][j-1];
//					win[7]=a[i+1][j-2];
				
           	  		 for(in=-2;in<3;in++)
           			{
           	  			for(jn=-2;jn<3;jn++)
           	  			{
           	  				temp5+=(float)(abs(a[i+in][j+jn]-a[i][j]));	
		        	    }
		        	}
		        	sigma=(float)sqrt((1.0+temp5)/(float)Sr);
		        	temp5=0.0;
		        	 printf("\nsigma:\n");


		        for(in=-2;in<3;in++)
           			{
           	  			for(jn=-2;jn<3;jn++)
           	  			{ 
					            if(in == 0 && jn == 0)
					            {
					                EE[i+in][j+jn]=0.0;
//					                printf("%.2f ",EE[i+in][j+jn]);	
						  		}
						  		else
						  		{
						  			Sij=(float)(fabs(max(in,jn)))/(float)(lamdaS);
//						  			printf("%.2f ",Sij);	
					            	Sig=(float)(fabs(a[i][j]-a[i+in][j+jn]))/((float)(lamdaG)*(float)sigma);
//					            	printf("%.2f ",Sig);
						  			EE[i+in][j+jn]=exp(-Sij-Sig);
//						  			printf("%.2f ",EE[i+in][j+jn]);				  
							    }						    
							} 
							temp5=0.0;  
						}		
					
					}
				}
				

	for(i=1;i<=r;i++)
	     	{   
		   	 for(j=1;j<=c;j++)
		    	{     
				   	win[0]=a[i-1][j-2];
					win[1]=a[i-2][j-1];    //5 * 5  filter window 
					win[2]=a[i-2][j+1];
					win[3]=a[i-1][j+2];
					win[4]=a[i+1][j+2];
					win[5]=a[i+2][j+1];
					win[6]=a[i+2][j-1];
					win[7]=a[i+1][j-2];
           	  			for(q=0;q<Sr-1;q++)
           	  			{ 	
							temp6+=((float)EE[i+in][j+jn]*(float)win[q]); 
							temp7+=EE[i+in][j+jn];
			            }
			        
			    	E[i][j]=(float)temp6/(float)temp7;
//			    	printf("%.2f ",E[i][j]);
			    	temp6=0.0;
			    	temp7=0.0;
			   }
	        }
    } 
    
   else
    {
    	for(i=1;i<=r;i++)
	     	{   temp5=0.0;
		   		for(j=1;j<=c;j++)
		    	{	      
		     for(in=-3;in<4;in++)
           		{ 
           	  		for(jn=-3;jn<4;jn++)
           	  			{ 
           	  		      temp5+=(float)(abs(a[i+in][j+jn]-a[i][j]));
		                }
		        }
		        	sigma=(float)sqrt((1+temp5)/(float)Sr);
		        	temp5=0.0;
//		        	 printf("\nsigma:\n");
//		        	 printf("%.2f ",sigma);
					  
//	            printf("\nEE:"); 
		        for(in=-3;in<4;in++)
           			{
           	  			for(jn=-3;jn<4;jn++)
           	  			{ 
					            if(in == 0 && jn == 0)
					            {
					                EE[i+in][j+jn]=0.0;
//					                printf("%.2f ",EE[i+in][j+jn]);	
						  		}
						  		else
						  		{
						  			Sij=(float)(fabs(max(in,jn)))/(float)(lamdaS);
//						  			printf("%.2f ",Sij);	
					            	Sig=(float)(fabs(a[i][j]-a[i+in][j+jn]))/((float)(lamdaG)*(float)sigma);
//					            	printf("%.2f ",Sig);
						  			EE[i+in][j+jn]=exp(-Sij-Sig);
//						  			printf("%.2f ",EE[i+in][j+jn]);				  
							    }
//							    printf("%.2f ",EE[i+in][j+jn]);	
							}   
						}		
					temp5=0.0;
					}
				}

	for(i=1;i<=r;i++)
	     	{   
		   	 for(j=1;j<=c;j++)
		    	{     
				   for(in=-3;in<4;in++)
           			{
           	  			for(jn=-3;jn<4;jn++)
           	  			{ 	
							temp6+=((float)EE[i+in][j+jn]*(float)a[i+in][j+jn]); 
							temp7+=EE[i+in][j+jn];
			            }
			        }
			    	E[i][j]=(float)temp6/(float)temp7;
//			    	printf("%.2f ",E[i][j]);
			    	temp6=0.0;
			    	temp7=0.0;
			   }
	        }
    } 

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
                	center[d]+=(float)(pow(fmv[a[i][j]][d],fmp)*((float)hist[a[i][j]])*((float)E[i][j]));
		   			fmvsum+=(float)(pow(fmv[a[i][j]][d],fmp)*hist[a[i][j]]);
				}
				             						//fmvsum is the summation (uij)^m  //in calculation of center[d]
            }
         
             center[d]=(center[d]/(fmvsum));
      }
			
//	        printf("%.2f ",center[d]);

		printf("\n\nCenters are: ");
			
		for(d=0;d<k;d++)	
		{
			printf("%.2f ",center[d]);
		}
//		printf("\n");
		//************************Distance calculation from the center*************************************// 
			
		
	    for(i=1;i<=r;i++)
           {
              for(j=1;j<=c;j++)
              { 
//			  scar = 0.0;sum_avg=0.0;
                 for(d=0;d<k;d++)
                  {

					dist[i][j][d]=(float)(fabs(E[i][j]-center[d]));

					Xpsum=0.0;

                  }
                 
              }
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
						  temp2+=(pow(dist[i][j][d],2))/(pow(dist[i][j][x],2));
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
			      if(fabs(avg[d]-center[d])>0.5)
				  {   
					flag=1;
					break;
				  }
			
			   }
	}
	
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
