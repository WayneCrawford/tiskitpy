//#include <malloc.h>

void  decim(tr1,n,dec_ratio,pos)
double  *tr1;
long *n;
long dec_ratio,pos;
{
    long j;
    long max_pos;
    double max;
    long ndat,ndat2;
    double *x;
    ndat = *n;
    max = fabs(tr1[0]);
    max_pos = 0;
    for (j = 0; j < ndat;j++)
    {
        if(fabs(tr1[j]) > max){
            max = fabs(tr1[j]);
            max_pos = j;
        }
    }
    /* only for negative start positions take the 
       position of the maximum */
    if (pos >= 0) max_pos = pos;  
  
    ndat2 = ndat/dec_ratio;
    x  = (double *)calloc(ndat2+1,sizeof(double));
    for (j = max_pos; j < ndat;j = j + dec_ratio)
    {
        if((j/dec_ratio < ndat2) && (j/dec_ratio>=0))
            x[j/dec_ratio] = tr1[j];
    }
    for (j = max_pos - dec_ratio; j >= 0;j = j - dec_ratio)
    {
        if((j/dec_ratio < ndat2) && (j/dec_ratio>=0))
            x[j/dec_ratio] = tr1[j];
    }
    *n = ndat2;
    /* copy data back to tr[] */
    for (j = 0; j < ndat;j++)
        tr1[j] = 0.0;
    for (j = 0; j < ndat2;j++)
        tr1[j] = x[j];
    free((char *)x);
   
}
