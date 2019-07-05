#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define MAXN 10   //x->t
#define g 9.81

void fun(double x,double y[],double RHS[],double K1,double K2,double b0,double m)
{
	RHS[0]=y[1];
	RHS[1]=((-2*K1*(1+K2*pow((-b0+sqrt(b0*b0+y[0]*y[0])),2.0))*y[0])/(m*sqrt(b0*b0+y[0]*y[0]))-g);
	//printf("\tRHS[0]=%lf\tRHS[1]=%lf\n",RHS[0],RHS[1]);
}

void vrk4(double x0,double y0[],double h,int n,void (*fun)(double,double*,double*,double,double,double,double),double y1[],double K1,double K2,double b0,double m)
{
	int		i;
	double	k1[MAXN],k2[MAXN],k3[MAXN],k4[MAXN];
	double	ytmp[MAXN];

	fun(x0,y0,k1,K1,K2,b0,m);
	for(i=0; i<n; ++i)
	{
		k1[i] *= h;
		ytmp[i] = y0[i] + k1[i]/2.0;
	}

	fun(x0+h/2.0,ytmp,k2,K1,K2,b0,m);
	for(i=0; i<n; ++i)
	{
		k2[i] *= h;
		ytmp[i] = y0[i] + k2[i]/2.0;
	}

	fun(x0+h/2.0,ytmp,k3,K1,K2,b0,m);
	for(i=0; i<n; ++i)
	{
		k3[i] *= h;
		ytmp[i] = y0[i] + k3[i];
	}

	fun(x0+h,ytmp,k4,K1,K2,b0,m);
	for(i=0; i<n; ++i)
		k4[i] *= h;

	for(i=0; i<n; ++i)
		y1[i] = y0[i] + (k1[i] + 2.*k2[i] + 2.*k3[i] + k4[i])/6.;

}



int main(void)
{
	double y0[2]={0.1,0.1}; //z0,v0
	double y1[2];           //wyniki z,v w danej chwili (po k-tym kroku calkowania)
	double t,t0,t1;         //czas
	double h;
	double E,Eps,Epc,Ek,A;   // krok calkowania,energie,pewien wspolczynnik
	double K1,K2,b0,m;      //dane do zadania
	double max,min,zs;  //skrajne wartosci z, aby ustalic punkt rownowagi
	int n,N;					//liczba rownan
	FILE* f=NULL;

	t0=0.0;
	t1=60.0;
	t=t0;
	h=0.01;                //60/0.125=480=N cos sensownego
	n=2;                    //2gi rzad RR
	N=(int)(fabs(t1-t0)/h);

	printf("Wprowadz dane do zadania:\n[N] k1=");
	scanf_s("%lf",&K1);
	printf("[1/m^2] k2=");
	scanf_s("%lf",&K2);
	printf("[m] b0="); //means 'a' in the picture
	scanf_s("%lf",&b0);
	printf("[kg] m=");//means 'M' in the picture
	scanf_s("%lf",&m);

	/*b0=0.25;  examples
	K1=1;
	K2=0.8;
	m=1;*/

	max=y0[0];
	min=y0[0];


	errno_t err_f =fopen_s(&f,"data.txt","w");
	if(err_f!=0)
	{
		printf("BLAD ODCZYTU");
		return 1;
	}

	fprintf(f,"t\tz\tvz\tm\t\t\t\t\t\tK1=%lf\tK2=%lf\tb0=%lf\tm=%lf\n\n",K1,K2,b0,m);

	for(int j=0;j<N+1;j++)
	{
		vrk4(t,y0,h,n,fun,y1,K1,K2,b0,m);
		A=-b0+sqrt(b0*b0+y1[0]*y1[0]);
		Eps=2*(K1*A+(1.0/3.0)*K1*K2*A*A*A); //E potencjalna sprezystosci
		Epc=m*g*y1[0];                      //E potencjalna ciezkosci
		Ek=0.5*m*y1[1]*y1[1];
		E=Eps+Epc+Ek;

		if(y1[0]>max)
			max=y1[0];
		if(y1[0]<min)
			min=y1[0];


		printf("Dla t=%lf\tz=%lf\tv=%lf\tE=%lf\n",t,y1[0],y1[1],E);
		fprintf(f,"%lf\t%lf\t%lf\t%lf\n",t,y1[0],y1[1],E);
		y0[0] = y1[0];
		y0[1] = y1[1];
		t+=h;
	}
	zs=(max+min)/2.0;  //mozna by pododawac zs do ztow zeby przesunac wykres w gore, ale zrobienie tego znacznie spowolniloby dzialanie programu
	printf ("Srodek rownowagi wynosi z=%lf\n",zs);
	fclose(f);
	return 0;
}




