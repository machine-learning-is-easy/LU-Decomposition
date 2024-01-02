void intpoly(float x,float pt[3][13],double re[2])

{

	int i;
	for(i=0;i<13;i++)
	{
		if(x<pt[0][i])
			break;
	}
	if (i!=0&&i!=12)
	{
		re[0]=pt[1][i]*(x-pt[0][i-1])/0.05+pt[1][i-1]*(pt[0][i]-x)/0.05;
		re[1]=pt[2][i]*(x-pt[0][i-1])/0.05+pt[2][i-1]*(pt[0][i]-x)/0.05;
	}
	if(i==0)
	{
		re[0]=pt[1][0];
		re[1]=pt[2][0];
	}
	if(i==12)
	{
		re[0]=pt[1][12];
		re[1]=pt[2][12];
	}

	
}
