#include "windows.h"
#include "stdio.h"
#include "head.h"
int main()
{
	FILE *f,*save;
	int i,j,time=0,floodingtime=0,tmax=1000;
	int flag=0;
	float num;
	char buffer[20];
	char c1[10],c2[10];
	float A,l,Q,N,dx,dt,por,sw;
	float uo,uw;
	float k;
	float swc,kro,krw;
	static float Kr[3][13];
	double tr[2];
	float M[39][39],L[39][39],U[39][39];
	float p1[40],p2[40],reso[39];
	float la[40];
	float con[39];
	float s1[39],s2[39];
	p1[39]=0;
	p2[39]=0;

	if((f=fopen("sf.txt","r"))==NULL)
	{
		printf("cannot open this file\n");
		exit(0);
	}
	if((save=fopen("ps.txt","w"))==NULL)
	{
		printf("cannot open this file\n");
		exit(0);
	}
	//read the initial data

	while(fscanf(f,"%s %f\n",buffer,&num))
	{
		if(buffer[0]=='/')
		{
			break;
		}
		else
		{
			switch(buffer[0])
			{
			case 'k':
				k=num;
				break;
			case 'A':
				A=num;
				break;
			case 'l':
				l=num;
				break;
			case 'Q':
				Q=num;
				break;
			case 'N':
				N=num;
				break;
			case 's':
				sw=num;
			case 'd':
				if(buffer[1]=='x')
					dx=num;
				else if(buffer[1]=='t')
					dt=num;
				break;
			case 'u':
				if(buffer[1]=='o')
					uo=num;			
				else if(buffer[1]=='w')
					uw=num;			
				break;
			case 'p':
				por=num;
				break;
			default:
				break;
			}
		}		
	}

	i=0;
	fscanf(f,"%s %s %s\n",buffer,c1,c2);
	while(fscanf(f,"%f %f %f\n",&swc,&krw,&kro)!=EOF)
	{
		Kr[0][i]=swc;
		Kr[1][i]=kro;
		Kr[2][i]=krw;
		i++;
	}
		fclose(f);

		for(i=0;i<40;i++)
		{
			p1[i]=0;
			s1[i]=sw;
		}
		//repetion start
		while(time<=tmax)
		{
			for(i=0;i<40;i++)
			{
				intpoly(s1[i],Kr,tr);
				kro=tr[0];
				krw=tr[1];
				la[i]=k*kro/uo+k*krw/uw;
			}
			for(i=0;i<39;i++)
				for(j=0;j<39;j++)
					M[i][j]=0;
			M[0][0]=1;
			M[0][1]=-1;
			con[0]=dx/A*Q/la[0];
			for(i=1;i<38;i++)
			{
				M[i][i-1]=la[i-1];
				M[i][i]=-la[i-1]-la[i];
				M[i][i+1]=la[i];
				con[i]=0;
			}
			M[i][i-1]=la[i-1];
			M[i][i]=-la[i-1]-la[i];
			con[i]=0;
			LU(M,L,U);
			reso[0]=con[0]/L[0][0];
			for(i=1;i<39;i++)
			{
				reso[i]=(con[i]-reso[i-1]*L[i][i-1])/L[i][i];
			}
			p2[38]=reso[38];
			for(i=37;i>=0;i--)
			{
				p2[i]=reso[i]-p2[i+1]*U[i][i+1];
			}
			intpoly(s1[0],Kr,tr);
			krw=tr[1];
			s2[0]=s1[0]+dt/(por*dx)*(Q/A+k*krw/uw*(p2[1]-p2[0])/dx);
			for (i=1;i<39;i++)
			{
				intpoly(s1[i],Kr,tr);
				num=tr[1];//n time
				intpoly(s1[i-1],Kr,tr);
				krw=tr[1];//n-1 time
				s2[i]=s1[i]+dt/(por*dx)*(k*num/uw*(p2[i+1]-p2[i])/dx-k*krw/uw*(p2[i]-p2[i-1])/dx);				
			}
			for(i=0;i<40;i++)
			{
				s1[i]=s2[i];
				p1[i]=p2[i];
			}
			time=dt+time;
			if(s2[38]>sw&&floodingtime==0)
			{	
				floodingtime=time;				
			}
								
			if(time%100==0&&time>=100)
			{			
				fprintf(save,"print %ds p distribution \n",time);
				fprintf(save,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
				for(j=0;j<40;j++)
				{
					fprintf(save,"%f  ",p2[j]);
				}
				fprintf(save," \n");
				fprintf(save,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
				fprintf(save,"print %ds s distribution \n",time);
				fprintf(save,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
				for(j=0;j<39;j++)
				{
					fprintf(save,"%f   ",s2[j]);
				}
				fprintf(save," \n");			
				fprintf(save,"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
		  }	
			//go to repetion;			
		}
		fprintf(save,"%d s the rim of waterflooding arrive at the product well\n",floodingtime);			
		fclose(save);
}
