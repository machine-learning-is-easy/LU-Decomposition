void LU(float M[39][39],float L[39][39],float U[39][39])
//comput the matrix L and U L*U=M
{
	int i,j,k;
	int row=0,colum=0;
	int rows=39;
	for(i=0;i<rows;i++)
		for(j=0;j<rows;j++)
		{
			L[i][j]=0;
			U[i][j]=0;

		}
	for(i=0;i<rows;i++)
		U[i][i]=1;

	L[0][0]=M[0][0];
	U[0][1]=M[0][1]/L[0][0];

	for(i=1;i<rows-1;i++)
	{
		L[i][i-1]=M[i][i-1];
		L[i][i]=M[i][i]-L[i][i-1]*U[i-1][i];
		U[i][i+1]=M[i][i+1]/L[i][i];
	}
	L[rows-1][rows-2]=M[rows-1][rows-2];
	L[rows-1][rows-1]=M[rows-1][rows-1]-L[rows-1][rows-2]*U[rows-2][rows-1];
}