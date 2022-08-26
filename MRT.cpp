#include <iostream>
#include <cmath>
//1.transfer distribution function to momentum space;
//2.0.collision in momentum space;
//3.transfer back to velcocity space after collsion;
//4.streaming as normal in velocity space.
//**********************************************************************************************//
int const mx=100,my=100;
double f[9][mx][my],fmom[9][mx][my];
double fmeq[9][mx][my],rho[mx][my],u[mx][my],v[mx][my];
double stmiv[9][9];
//1.define M matrix, M inverse matrix
double tm[9][9]=
        {
            {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0},
            {-4.0,-1.0,-1.0,-1.0,-1.0,2.0,2.0,2.0,2.0},
            {4.0,-2.0,-2.0,-2.0,-2.0,1.0,1.0,1.0,1.0},
            {0.0,1.0,0.0,-1.0,0.0,1.0,-1.0,-1.0,1.0},
            {0.0,-2.0,0.0,2.0,0.0,1.0,-1.0,-1.0,1.0},
            {0.0,0.0,1.0,0.0,-1.0,1.0,1.0,-1.0,-1.0},
            {0.0,0.0,-2.0,0.0,2.0,1.0,1.0,-1.0,-1.0},
            {0.0,1.0,-1.0,1.0,-1.0,0.0,0.0,0.0,0.0},
            {0.0,0.0,0.0,0.0,0.0,1.0,-1.0,1.0,-1.0}
        };
double a=1.0/36.0;
double  tminv[9][9]=
        {
            {4.0*a, -4.0*a, 4.0*a, 0.0*a, 0.0*a,0.0*a,0.0*a,0.0*a,0.0*a},
            {4.0*a,-1.0*a,-2.0*a, 6.0*a, -6.0*a,0.0*a,0.0*a,9.0*a,0.0*a},
            {4.0*a,-1.0*a,-2.0*a, 0.0*a, 0.0*a,6.0*a,-6.0*a,-9.0*a,0.0*a},
            {4.0*a,-1.0*a,-2.0*a, -6.0*a,6.0*a,0.0*a,0.0*a,9.0*a,0.0*a},
            {4.0*a,-1.0*a,-2.0*a, 0.0*a, 0.0*a,-6.0*a,6.0*a,-9.0*a,0.0*a},
            {4.0*a,2.0*a,  1.0*a,  6.0*a,3.0*a,6.0*a,3*a,0.0*a,9.0*a},
            {4.0*a,2.0*a,  1.0*a,  -6.0*a,-3.0*a,6.0*a,3.0*a,0.0*a,-9.0*a},
            {4.0*a,2.0*a,  1.0*a,  -6.0*a,-3.0*a,-6.0*a,-3.0*a,0.0*a,9.0*a},
            {4.0*a,2.0*a,  1.0*a,  6.0*a, 3.0*a,-6.0*a,-3.0*a,0.0*a,-9.0*a}
        };
//2.define diagonal matrix S
double alpha=0.17;
double tau=(1+6*alpha)/2;
double s[9]={1.0,1.4,1.4,1.0,1.2,1.0,1.2,1.0/tau,1.0/tau};
void collision()
{
int i,j,k,l;
for(i=0;i<mx;i++)
{
    for(j=0;j<my;j++)
    {
        for(k=0;k<9;k++)
        {
            double sumb =0.0;
            for(l=0;l<9;l++)
            {
                sumb=sumb+stmiv[k][l]*(fmom[l][i][j]-fmeq[l][i][j]);
            }
            f[k][i][j]=f[k][i][j]-sumb;
        }
    }
}
}
int main()
{
//3.calculate M^-1*S
for (int i = 0; i < 9; i++)
{
    for (int j = 0; j < 9; j++)
    {
        stmiv[i][j]=tminv[i][j]* s[j];
    }
}
//4.caculate the equilibrium distribution function in momentum space, distribution fucntion M*f
int i,j,k,l;
for(i=0;i<mx;i++)
{
    for(j=0;j<my;j++)
    {
        for(k=0;k<9;k++)
        {
            double suma =0.0;
            for(l=0;l<9;l++)
            {
                suma=suma+tm[k][l]*f[l][i][j];
                fmom[k][i][j]=suma;
            }
        }
    }
}
for(i=0;i<mx;i++)
{
    for(j=0;j<my;j++)
    {
        fmeq[0][i][j]=rho[i][j];
        fmeq[1][i][j]=rho[i][j]*(-2.0+3.0*rho[i][j]*(u[i][j]*u[i][j]+v[i][j]*v[i][j]));
        fmeq[2][i][j]=rho[i][j]*(1.0-3.0*rho[i][j]*(u[i][j]*u[i][j]+v[i][j]*v[i][j]));
        fmeq[3][i][j]=rho[i][j]*u[i][j];
        fmeq[4][i][j]=-rho[i][j]*u[i][j];
        fmeq[5][i][j]=rho[i][j]*v[i][j];
        fmeq[6][i][j]=-rho[i][j]*v[i][j];
        fmeq[7][i][j]=rho[i][j]*(u[i][j]*u[i][j]-v[i][j]*v[i][j]);
        fmeq[8][i][j]=rho[i][j]*u[i][j]*v[i][j];
    }
}
//5.collision f*=f-M^-1*S(fmom-fmeq)
collision();
//identity matrix
double ev[9][9]={0};
for (int i = 0; i < 9; i++)
{
    for (int j = 0; j < 9; j++)
    {
        for (int k = 0; k < 9; k++)
        {
            ev[i][j]+=tm[i][k]*tminv[k][j];
        }
        
    }
}
for (int i = 0; i < 9; i++)
{
    for (int j = 0; j < 9; j++)
    {
        std::cout<<(int)(ev[i][j]+0.001)<<" ";
    }
    std::cout<<std::endl;
}
    return 0;
}