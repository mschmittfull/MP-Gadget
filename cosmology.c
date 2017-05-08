#include <math.h>
#include "allvars.h"
#include "cosmology.h"
#include "utils-string.h"

/*Hubble function at scale factor a, in dimensions of All.Hubble*/
double hubble_function(double a)
{
    /* for compatibility only */
    return All.Hubble * HubbleEa(a);
}

static double growth(double a);
static double growth_int(double a, void * params);

static double growth_int(double a, void * params)
{
    /* FIXME: taylor expand for small a */
    if(a == 0) return 0;
    return pow(1 / (a * HubbleEa(a)), 3);
}

double F_Omega(double a)
{
    return pow(OmegaA(a), 0.6);
}

static double growth(double a)
{
    /* NOTE that the analytic COLA growthDtemp() is 6 * pow(1 - c.OmegaM, 1.5) times growth() */

    int WORKSIZE = 100000;

    double result, abserr;
    gsl_integration_workspace *workspace;
    gsl_function F;


    workspace = gsl_integration_workspace_alloc(WORKSIZE);

    F.function = &growth_int;
    F.params = (double[]) {All.CP.Omega0, All.CP.OmegaLambda};

    gsl_integration_qag(&F, 0, a, 0, 1.0e-9, WORKSIZE, GSL_INTEG_GAUSS41, 
            workspace, &result, &abserr);

    gsl_integration_workspace_free(workspace);

    return HubbleEa(a) * result;
}

double OmegaA(double a) {

    /* FIXME: radiation is not there! */

    return All.CP.Omega0 / (All.CP.Omega0 + a * (1 - All.CP.Omega0 - All.CP.OmegaLambda) + a * a * a * All.CP.OmegaLambda);
}

double GrowthFactor(double a) { // growth factor for LCDM
    return growth(a) / growth(1.0);
}

double DLogGrowthFactor(double a) {
    /* Or OmegaA^(5/9) */
    return pow(OmegaA(a), 5.0 / 9);
}

double GrowthFactor2(double a) {
    /* Second order growth factor */
    /* 7 / 3. is absorbed into dx2 */
    double d = GrowthFactor(a);
    return d * d * pow(OmegaA(a) / OmegaA(1.0), -1.0/143.);
}

double DLogGrowthFactor2(double a) {
    return 2 * pow(OmegaA(a), 6.0/11.);
}

double HubbleEa(double a)
{
    /* H(a) / H0 */
    double hubble_a;
    /* first do the terms in SQRT */
    hubble_a = All.CP.OmegaLambda;

    hubble_a += All.CP.OmegaK / (a * a);
    hubble_a += All.CP.Omega0 / (a * a * a);

    if(All.CP.RadiationOn) {
        hubble_a += All.CP.OmegaG / (a * a * a * a);
        /* massless neutrinos are added only if there is no (massive) neutrino particle.*/
        if(!NTotal[2])
            hubble_a += All.CP.OmegaNu0 / (a * a * a * a);
    }
    return hubble_a;
}

double DHubbleEaDa(double a) {
    /* FIXME: add radiation ! */
    /* d E / d a*/
    double E = HubbleEa(a);
    return 0.5 / E * (-3 * All.CP.Omega0 / (a * a * a * a));
}
double D2HubbleEaDa2(double a) {
    /* FIXME: add radiation ! */
    double E = HubbleEa(a);
    double dEda = DHubbleEaDa(a);
    return - dEda * dEda / E + dEda * (-4 / a);
}
double DGrowthFactorDa(double a) {
    /* FIXME: add radiation ! */
    double E = HubbleEa(a);

    double EI = growth(1.0);

    double t1 = DHubbleEaDa(a) * GrowthFactor(a) / E;
    double t2 = E * pow(a * E, -3) / EI;
    return t1 + t2;
}

double D2GrowthFactorDa2(double a) {
    /* FIXME: add radiation ! */
    double d2Eda2 = D2HubbleEaDa2(a);
    double dEda = DHubbleEaDa(a);
    double E = HubbleEa(a);
    double EI = growth(1.0);
    double t1 = d2Eda2 * GrowthFactor(a) / E;
    double t2 = (dEda + 3 / a * E) * pow(a * E, -3) / EI;
    return t1 - t2;
}

static double sigma2_int(double k, void * p)
{
    void ** params = p;
    FunctionOfK * fk = params[0];
    double * R = params[1];
    double kr, kr3, kr2, w, x;

    kr = *R * k;
    kr2 = kr * kr;
    kr3 = kr2 * kr;

    if(kr < 1e-8)
        return 0;

    w = 3 * (sin(kr) / kr3 - cos(kr) / kr2);
    x = 4 * M_PI * k * k * w * w * function_of_k_eval(fk, k);

    return x;
}

FunctionOfK * function_of_k_new_from_string(const char * string, int logscale)
{
    FunctionOfK * fk = NULL;
    char ** list = fastpm_strsplit(string, "\n");
    char ** line;
    int i;
    int pass = 0;
    /* two pass parsing, first pass for counting */
    /* second pass for assignment */
    while(pass < 2) {
        i = 0;
        for (line = list; *line; line++) {
            double k, p;
            if(2 == sscanf(*line, "%lg %lg", &k, &p)) {
                if(logscale) {
                    k = pow(10, k);
                    p = pow(10, p);
                }
                if(pass == 1) {
                    fk->table[i].k = k;
                    fk->table[i].P = p;
                }
                i ++;
            }
        }

        if(pass == 0) {
            fk = malloc(sizeof(*fk) + sizeof(fk->table[0]) * i);
            fk->bytesize = sizeof(*fk) + sizeof(fk->table[0]) * i;
            fk->size = i;
            fk->normfactor = 1;
        }
        pass ++;
    }

    free(list);

    return 0;
}

double function_of_k_eval(FunctionOfK * fk, double k)
{
    /* ignore the 0 mode */

    if(k == 0) return 1;

    int l = 0;
    int r = fk->size - 1;

    while(r - l > 1) {
        int m = (r + l) / 2;
        if(k < fk->table[m].k)
            r = m;
        else
            l = m;
    }
    double k2 = fk->table[r].k,
           k1 = fk->table[l].k;
    double p2 = fk->table[r].P,
           p1 = fk->table[l].P;

    if(l == r) {
        return fk->table[l].P;
    }

    if(p1 == 0 || p2 == 0 || k1 == 0 || k2 == 0) {
        /* if any of the p is zero, use linear interpolation */
        double p = (k - k1) * p2 + (k2 - k) * p1;
        p /= (k2 - k1);
        return p;
    } else {
        k = log(k);
        p1 = log(p1);
        p2 = log(p2);
        k1 = log(k1);
        k2 = log(k2);
        double p = (k - k1) * p2 + (k2 - k) * p1;
        p /= (k2 - k1);
        return exp(p);
    }
}

double function_of_k_tophat_sigma(FunctionOfK * fk, double R)
{
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    void * params[] = {fk, &R};
    double result,abserr;
    gsl_function F;
    F.function = &sigma2_int;
    F.params = params;

    /* note: 500/R is here chosen as integration boundary (infinity) */
    gsl_integration_qags (&F, 0, 500. / R, 0, 1e-4,1000,w,&result, &abserr);
    //   printf("gsl_integration_qng in TopHatSigma2. Result %g, error: %g, intervals: %lu\n",result, abserr,w->size);
    gsl_integration_workspace_free (w);
    return sqrt(result);
}

void function_of_k_normalize_sigma(FunctionOfK * fk, double R, double sigma) {
    double old = function_of_k_tophat_sigma(fk, R);
    int i;
    for(i = 0; i < fk->size; i ++) {
        fk->table[i].P *= sigma / old;
    };
}


