#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<limits.h>
#define fmp 2
#define Sr 9
#define max(a,b)    (((a>b))?(a):(b))
//#define k 4

//************WARNING !!!!!!!!!! (no. of clusters to be formed) k's value will be taken during runtime!!!!***********//

//*********P2 256 256 256 will be automatically generated on the output file************//

int main()
{
	
	int i,j,k,n,d,r,c,p;srand(time(0));char str[10];
	int m,pos,f=0; float *avg=NULL,*center=NULL;int fm;
	float fmvsum=0.0,sum =0.0,mini;	int flag = 1;int win[10]; float Xpsum=0.0;int q;int in,jn;
    int lx; int count = 0;float mu=0.0,tao=0.0;float **eta=NULL,**lamdaS=NULL,**lamdaG=NULL;int h;int maxi=INT_MIN;
	float **sigma=NULL;float temp5=0.0,deno=0.0;;float **C=NULL;int x;float arr[5];
//	int *hist = NULL;
	float **E=NULL; int s=0; int g;
	int **a=NULL,**b=NULL;
	float ***dist=NULL;float scar = 0.0;
	float ***temp=NULL;


	FILE *fp1 = fopen("m7_gray_fcm_nls_hist.txt","r");
	FILE *fp2 = fopen("m7_gray_cva_asifc.pgm","wb");
	
	fscanf(fp1,"%s %d %d %d",str,&c,&r,&p);
	fprintf(fp2,"%s %d %d %d\n",str,c,r,p);
	
	printf("\n\nEnter the value of k: ");
	scanf("%d",&k);                           //Input no. of clusters
	
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
	
	C=(float **)calloc(r,sizeof(float *));
    for(i=1;i<=r;i++) C[i]=(float *)calloc(c,sizeof(float));
	
	lamdaS=(float **)calloc(r,sizeof(float *));
    for(i=1;i<=r;i++) lamdaS[i]=(float *)calloc(c,sizeof(float));
	
	lamdaG=(float **)calloc(r,sizeof(float *));
    for(i=1;i<=r;i++) lamdaG[i]=(float *)calloc(c,sizeof(float));
	sigma=(float **)calloc(r,sizeof(float *));
    for(i=1;i<=r;i++) sigma[i]=(float *)calloc(c,sizeof(float));
	
	eta=(float **)calloc(r,sizeof(float *));
    for(i=1;i<=r;i++) eta[i]=(float *)calloc(c,sizeof(float));
	
	dist = (float***)calloc(r,sizeof(float**));
	for(i=1;i<=r;i++)
	{
		dist[i]=(float**)calloc(c,sizeof(float*));
		for(j=1;j<=c;j++)
		{
			dist[i][j]=(float*)calloc(k,sizeof(float));
		}
	}
	
	temp = (float***)calloc(r,sizeof(float**));
	for(i=1;i<=r;i++)
	{temp[i]=(float**)calloc(c,sizeof(float*));
		for(j=1;j<=c;j++)
		{
			temp[i][j]=(float*)calloc(k,sizeof(float));
		}
	}
	//*************************Input matrix of image************************************//	
    
   for(i=1;i<=r;i++) for(j=1;j<=c;j++) fscanf(fp1,"%d",&a[i][j]);
		//to read the input image

	//*************************matrix padding************************************************//
	
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
	//********************random center initialize first time*****************************************//
	
	center=(float *)calloc(k,sizeof(float));
    avg=(float *)calloc(k,sizeof(float));
    
    for(i=0;i<k;i++)
		{	
		    center[i]=rand()%(p+1);             //0-255
            printf("%2f ",center[i]);
		}
		
	//****************************************calculation of lamda smoothing factor****************************************************//

       for(i=1;i<=r;i++)
	     {   
		   for(j=1;j<=c;j++)
		    {	
			   E[i][j]=1+((float)1/(float)(r*c));
//			   printf("%.2f ",E[i][j]);
		    }
	     }

	float t1=0.0,t2=0.0,t3=0.0;
	for(i=1;i<=r;i++)
	     {   
		   for(j=1;j<=c;j++)
		    {
		    	for(in=-1;in<2;in++) 
				{
				  for(jn=-1;jn<2;jn++) 
				  {
				      t1+=fabs(a[i][j]-a[i+in][j+jn]);
				  }
			    }
		        mu=((float)1/(float)Sr)*t1;
		        t1=0.0;
		        
		        for(in=-1;in<2;in++)
                   for(jn=-1;jn<2;jn++)
           	          t2+=pow(((a[i][j]-a[i+in][j+jn])-mu),2);
           		eta[i][j]=(float)(0.5+pow((((float)1/(float)Sr)*t2),0.5));
           		t2=0.0;
//           		printf("%.2f ",eta[i][j]);
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
						   t3+=fabs(a[i+in][j+jn]); 
						}
           	    }
           	    tao=((float)1/(float)Sr)*t3;
//           	               	    printf("%f ",tao);
           	    t3=0.0;
           	    
           	    for(in=-1;in<2;in++)
				{
           	       for(jn=-1;jn<2;jn++)
					  {
           	           	  sigma[i][j]+=pow((((float)1/(float)Sr)*tao),0.5); 
//						  printf("%.2f ",sigma[i][j]);
					  }
				}
           	           
           	   	lamdaS[i][j]=(float)((eta[i][j])/(1.0+sigma[i][j]));
//           	   	printf("%.2f ",lamdaS[i][j]);
           	   	if(lamdaS[i][j]>maxi)
           	   	maxi=lamdaS[i][j];
           	}	
        }
        printf("\n%d",maxi);
        for(i=1;i<=r;i++)   
		   for(j=1;j<=c;j++)
		     lamdaG[i][j]=(float)(0.1+(lamdaS[i][j]/maxi));
//		    	printf("%.2f ",lamdaG[i][j]);
        
  	float temp3=0.0,temp4=0.0; 
	       
	//************************start of while****************************************//    	          
while(flag)
	{   
	  for(d=0;d<k;d++)
            {   
			    avg[d]=center[d];
            	center[d]=0.0;
			}
	  Xpsum=0.0;
	  
	 for(d=0;d<k;d++)
	  {
		 fmvsum=0.0;
		for(i=1;i<=r;i++)
	     {   
		   for(j=1;j<=c;j++)
		    {

              for(in=-1;in<2;in++)
              {
               	for(jn=-1;jn<2;jn++)
               	{
           	    	Xpsum+=0.5+(a[i+in][j+jn]);
           	    }
           	  }
//           		printf("%.2f ",Xpsum);
			      
			    center[d]+=(float)(pow(fmv[a[i][j]][d],fmp)*((float)E[i][j])*(lamdaG[i][j]*a[i][j] + (((float)fabs((1.0-lamdaG[i][j]))/(float)Sr)*Xpsum)));
			   	fmvsum+=(float)(pow(fmv[a[i][j]][d],fmp)*E[i][j]);    
				Xpsum=0.0;	                                 
            }
			 	Xpsum=0.0;                                            
	      }
	        center[d]=(center[d]/(fmvsum));
//	        printf("%.2f ",center[d]);
	    }
	    
		printf("\n\nCenters are: ");
			
		for(d=0;d<k;d++) 
		{
		   printf("%.2f ",center[d]);
	    } 

		temp3=0.0,temp4=0.0;deno=0.0;
	 	for(d=0;d<k;d++)
	 	{ 
			for(i=1;i<=r;i++)
        		{
	      			for(j=1;j<=c;j++)
	        		    { 
						   deno+=(float)(fmv[a[i][j]][d]*E[i][j]);	
					    }
	     		}
//	     		printf("%.2f ",deno);
	     		arr[d]=deno;
	     		deno=0.0;
	 	}
//	 	deno=0.0;
//	 	printf("%.2f ",temp4);

	 	for(i=1;i<=r;i++)
        {
	      for(j=1;j<=c;j++)
	        {     
	           	for(d=0;d<k;d++)
					{
					    temp3+=(fmv[a[i][j]][d]*(-1)*log((float)fmv[a[i][j]][d]/(float)arr[d]));
				    }
	        	
//	        	printf("%.2f ",temp3);
	        	C[i][j]=(float)(exp(temp3));
//				printf("%.2f ",C[i][j]);	
	        	temp3=0.0;
	        	temp5+=C[i][j];
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
	 
		for(i=1;i<=r;i++)
        {
	      for(j=1;j<=c;j++)
	        {
	        	E[i][j]=(float)(1+((E[i][j]*C[i][j])/(float)temp5));
//	        	printf("%.2f ",E[i][j]);
	        }
	    }
	 		temp5=0.0;
	 	//************************Distance calculation from the center*************************************// 
			
		
	    for(i=1;i<=r;i++) for(j=1;j<=c;j++) for(d=0;d<k;d++) dist[i][j][d]=(float)(fabs(a[i][j]-center[d]));
      
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
					flag=1; break;
				  }
			   }
	}
	
	//***********************************Partition Coefficient Calculation(Vpc)**************************//
	float Vpc=0.0,Vpe=0.0;
	for(d=0;d<k;d++) for(i=1;i<=r;i++) for(j=1;j<=c;j++) Vpc+=(float)pow(fmv[a[i][j]][d],2);
	
	Vpc=(float)((Vpc)/(r*c));
	printf("\n\nThe Partition Coefficient for ASIFC(Vpc): %.2f ",Vpc);

	//***********************************Partition Entropy Calculation(Vpe)**************************//
	for(d=0;d<k;d++) for(i=1;i<=r;i++) for(j=1;j<=c;j++) Vpe+=(float)(fmv[a[i][j]][d]*log(fmv[a[i][j]][d]));
	
	Vpe=-((Vpe)/(r*c));
	printf("\n\nThe Partition Entropy for ASIFC(Vpe): %.2f ",Vpe);
	
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
					if(center[pos]<0)
					b[i][j]=(int)((-1)*center[pos]);
					else
					b[i][j]=(int)(center[pos]);
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
