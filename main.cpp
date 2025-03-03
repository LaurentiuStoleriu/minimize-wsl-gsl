#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <random>
//#include <windows.h>
#include <sys/time.h>

unsigned int timeGetTime()
{
	struct timeval now;
	gettimeofday(&now, NULL);
	return now.tv_usec/1000;
}


#define metode_fdf 1
//#undef metode_fdf

constexpr int Npart = 6;				// numar de particule
constexpr int Npasi = 1000;				// numar de pasi (de pozitii de echilibru)
constexpr double l0 = 0.4;				// lungimea resortului
constexpr double r_mic = 0.25;			// raza mica
constexpr double R_mare = 0.45;			// raza mare
constexpr double A = 0.001;				// "adancimea" gropii de potential
constexpr double perioada = 5.0*M_PI;	// perioada potentialului de suprafata
constexpr double k_el = 1.0;			// constanta elastica

double p[Npart], r[Npart];	// sirurile de pozitii si de raze

double v(double loco_p);
double fn1(const gsl_vector q[], void *params);
void dfn1(const gsl_vector q[], void *params, gsl_vector *df);
void fdfn1(const gsl_vector q[], void *params, double *f, gsl_vector *df);

int main(void)
{
	int j;
	//DWORD starttime, elapsedtime;     //windows
	unsigned int starttime, elapsedtime;     //windows

	FILE *fp;
	char numefis[200];

	// varianta mea favorita de generat numere aleatorii conform C++11 e cu std::random_device
	// (varianta C-clasic, cu functia rand() e invechita si n-ar mai trebui folosita)
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(-0.01, 0.01);
	std::uniform_int_distribution<int> decizie(0, 1);
	// la fiecare apel de tipul dist(mt) vom primi un double intre -0.01 si +0.01 (distrib. uniforma)


	p[0] = 0.4 /*+ dist(mt)*/;
	r[0] = r_mic;
	for (j = 1; j < Npart; j++)
	{
		r[j] = r_mic;
		p[j] = p[j - 1] + r[j - 1] + l0 + r[j] /*+ dist(mt)*/;
	}

	const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_fr;
	gsl_multimin_fdfminimizer* s;
	gsl_vector* x;
	gsl_multimin_function_fdf minex_func;

	size_t iter = 0;
	int status;
	double size;

	/* Starting point */
	x = gsl_vector_alloc(Npart);
	for (int i = 0; i < Npart; i++)
	{
		gsl_vector_set(x, i, p[i]);
	}

	double *par = NULL;
	/* Initialize method and iterate */
	minex_func.n = Npart;
	minex_func.f = fn1;
	minex_func.df = dfn1;
	minex_func.fdf = fdfn1;
	minex_func.params = par;

	s = gsl_multimin_fdfminimizer_alloc(T, Npart);
	gsl_multimin_fdfminimizer_set(s, &minex_func, x, 0.01, 1.0e-4);

	//////////////////////////////////////////////////////////////////////////
	//echilibru initial
	//////////////////////////////////////////////////////////////////////////
	starttime = timeGetTime();
	do
	{
		iter++;

		status = gsl_multimin_fdfminimizer_iterate(s);

		if (status)
			break;

		status = gsl_multimin_test_gradient(s->gradient, 1.0e-3);
		printf("%5zd %10.3e %10.3e f() = %7.3f\n", iter, gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1), s->f);

	} while (status == GSL_CONTINUE);

	for (int i = 0; i < Npart; i++)
	{
		p[i] = gsl_vector_get(s->x, i);				// salvare pozitii in p[]
		gsl_vector_set(x, i, p[i]);					// si in x pentru ca s->x != x
	}

	elapsedtime = timeGetTime() - starttime;
	printf("\n echilibru initial - DONE IN %d milliseconds\n", elapsedtime);

	sprintf(numefis, "/home/lali/TITAN-ROG-sync/c/rez/minim_%d.dat", Npart);
	fp = fopen(numefis, "w");
	for (int i = 0; i < Npart; i++)
	{
		fprintf(fp, "%d %20.16lf %4.2lf\n", i, p[i], r[i]);
	}
	fclose(fp);

	//////////////////////////////////////////////////////////////////////////
	//comutari
	//////////////////////////////////////////////////////////////////////////
	int contor_particule_comutate = 0;
// 	for (int i=0; i<Npart; i++)
// 	{
// 		if (decizie(mt))
// 		{
// 			r[i] = R_mare;
// 			contor_particule_comutate++;
// 		}
// 	}
	r[2] = R_mare;
	r[3] = R_mare;
	contor_particule_comutate += 2;
	printf("am comutat %d\n", contor_particule_comutate);

	status = gsl_multimin_fdfminimizer_restart(s);					// !!! neaparat restart!
	gsl_multimin_fdfminimizer_set(s, &minex_func, x, 0.01, 1.0e-8);

	starttime = timeGetTime();
	do
	{
		iter++;
		status = gsl_multimin_fdfminimizer_iterate(s);
	
		if (status)
			break;

		status = gsl_multimin_test_gradient(s->gradient, 1.0e-3);
		printf("%5zd %10.3e %10.3e f() = %7.3f\n", iter, gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1), s->f);

	} while (status == GSL_CONTINUE);

	for (int i = 0; i < Npart; i++) // salvare pozitii in p[]
	{
		p[i] = gsl_vector_get(s->x, i);
	}

	elapsedtime = timeGetTime() - starttime;
	printf("\n comutari - DONE IN %d milliseconds\n", elapsedtime);

	sprintf(numefis, "/home/lali/TITAN-ROG-sync/c/rez/minim_%d_sw.dat", Npart);
	fp = fopen(numefis, "w");
	for (int i = 0; i < Npart; i++)
	{
		fprintf(fp, "%d %20.16lf %4.2lf\n", i, p[i], r[i]);
	}
	fclose(fp);

	gsl_vector_free(x);
	gsl_multimin_fdfminimizer_free(s);

	return status;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

double v(double loco_p)
{
	return -A * sin(perioada * loco_p);
}

double fn1(const gsl_vector q[], void *params)
{
	(void)(params); /* avoid unused parameter warning */
	double sv = 0, sx = 0, f;
	double xi, xip1, delta;
	
	xi = gsl_vector_get(q, 0);
	for (int i = 0; i < (Npart - 1); i++)
	{
		//xi = gsl_vector_get(q, i);
		xip1 = gsl_vector_get(q, (i + 1));
		sv = sv + v(xi);
		delta = (xip1 - xi - r[i] - r[i + 1] - l0);
  		sx = sx + delta * delta;
		xi = xip1;
	}

	sv += v(gsl_vector_get(q, (Npart-1)));

	f = sv + (k_el / 2.0) * sx;
	return f;
}

void dfn1(const gsl_vector q[], void *params, gsl_vector *df)
{
	(void)(params); /* avoid unused parameter warning */
	double xim1, xi, xip1, deriv_xi;

	xi   = gsl_vector_get(q, 0);
	xip1 = gsl_vector_get(q, 1);
	deriv_xi = -perioada * A * cos(perioada * xi) - k_el * (xip1 - xi - r[0] - r[1] - l0);
	gsl_vector_set(df, 0, deriv_xi);

	for (int i = 1; i < (Npart - 1); i++)
	{
		xim1 = xi;
		xi   = xip1;
		xip1 = gsl_vector_get(q, (i + 1));
		deriv_xi = -perioada * A * cos(perioada * xi) - k_el * (r[i-1] - r[i+1] + xim1 - 2.0 * xi + xip1);
		gsl_vector_set(df, i, deriv_xi);
	}

	xim1 = xi;
	xi = xip1;
	deriv_xi = -perioada * A * cos(perioada * xi) + k_el * (xi - xim1 - r[Npart-2] - r[Npart-1] - l0);
	gsl_vector_set(df, (Npart-1), deriv_xi);
}

void fdfn1(const gsl_vector q[], void *params, double *f, gsl_vector *df)
{
	*f = fn1(q, params);
	dfn1(q, params, df);
}