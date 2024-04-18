//#include <malloc.h>

/* E. Wielandt's filter coefficients for interpolation ratios 1 - 6 */
/* Only the ratios 2, 3, and 5 are used presently.E. g. an interpolation */
/* by a factor of 4  is performed as 2 x 2 interpolations */

/* g2_1 RMS rel. error for k_rel <0.8: 0.009% */
static double g2_1[] = { -0.0002616,  0.0009302, -0.0023258,  0.0049142, -0.0092537,  0.0161355, -0.0266150,  0.0424554, -0.0670612,  0.1091686, -0.2008408,  0.6327541,  0.6327542, -0.2008408,  0.1091686, -0.0670612,  0.0424554, -0.0266150,  0.0161355, -0.0092537,  0.0049142, -0.0023258,  0.0009302, -0.0002616};

/* g3_1 RMS rel. error for k_rel <0.8: 0.008% */
static double g3_1[] = { -0.0002324,  0.0008256, -0.0020654,  0.0043684, -0.0082410,  0.0144087, -0.0238643,  0.0383080, -0.0611438,  0.1015046, -0.1959616,  0.8226883,  0.4111037, -0.1564917,  0.0885536, -0.0553522,  0.0353777, -0.0223077,  0.0135759, -0.0078056,  0.0041522, -0.0019670,  0.0007871, -0.0002211};

/* g4_1 RMS rel. error for k_rel <0.8: 0.007% */
static double g4_1[] = { -0.0001922,  0.0006827, -0.0017085,  0.0036157, -0.0068275,  0.0119544, -0.0198421,  0.0319579, -0.0512943,  0.0860733, -0.1708677,  0.8964090,  0.2985412, -0.1217245,  0.0701113, -0.0441759,  0.0283604, -0.0179323,  0.0109329, -0.0062935,  0.0033505, -0.0015879,  0.0006355, -0.0001784};

/* g5_1 RMS rel. error for k_rel <0.8: 0.006% */
static double g5_1[] = { -0.0001611,  0.0005720, -0.0014316,  0.0030307, -0.0057263,  0.0100352, -0.0166788,  0.0269186, -0.0433567,  0.0732492, -0.1480766,  0.9320452,  0.2327664, -0.0984035,  0.0572469, -0.0362360,  0.0233232, -0.0147709,  0.0090151, -0.0051932,  0.0027660, -0.0013112,  0.0005248, -0.0001473};


/* g5_2 RMS rel. error for k_rel <0.8: 0.009% */
static double g5_2[] = { -0.0002526,  0.0008977, -0.0022452,  0.0047467, -0.0089479,  0.0156272, -0.0258392,  0.0413715, -0.0657512,  0.1082703, -0.2048060,  0.7525200,  0.5015040, -0.1790148,  0.0997642, -0.0619420,  0.0394431, -0.0248145,  0.0150789, -0.0086611,  0.0046043, -0.0021804,  0.0008723, -0.0002452};

/* g6_1 RMS rel. error for k_rel <0.8: 0.005% */
static double g6_1[] = { -0.0001378,  0.0004891, -0.0012243,  0.0025926, -0.0049006,  0.0085933, -0.0142952,  0.0231041, -0.0373015,  0.0633145, -0.1296402,  0.9518885,  0.1901555, -0.0822067,  0.0481311, -0.0305567,  0.0197009, -0.0124900,  0.0076282, -0.0043963,  0.0023423, -0.0011105,  0.0004445, -0.0001247};

static int g2_1_len=24;
static int g3_1_len=24;
static int g4_1_len=24;
static int g5_1_len=24;
static int g5_2_len=24;
static int g6_1_len=24;

/* interpolation by factor 2 */
void  ipol2(xin,xout,n)
double  *xin, *xout;
long n;
{
    double xx;
    long i,j, k;
    /* xx = interpolated sample at center between index i and i+1 */
    for(i = 0; i < n-1; i++)
    {
        xx =0.0;
        for(j = 0; j < g2_1_len/2; j++)
        {
            if((i-j) >= 0) /* left side of filter */
               xx+=g2_1[g2_1_len/2-(j+1)]*xin[i-j];
            if((i+1+j) < n) /* right side of filter */
               xx+=g2_1[g2_1_len/2+j]*xin[i+1+j];
        }
        xout[2*i] = xin[i];
        xout[2*i+1] = xx;
        xout[2*i+2] = xin[i+1];
    }
}

/* interpolation by factor 3 */
void  ipol3(xin,xout,n)
double  *xin, *xout;
long n;
{
    double xx;
    long i,j, k;
    long ifac = 3; /* interpolation factor */
    
    for(i = 0; i < n-1; i++)
    {
        xout[ifac*i] = xin[i];
        xout[ifac*i+ifac] = xin[i+1];
    }

    /* here xx = interpolated sample 1/3 between index i and i+1 */
    for(i = 0; i < n-1; i++)
    {
        xx =0.0;
        for(j = 0; j < g3_1_len/2; j++)
        {
            if((i-j) >= 0) /* left side of filter */
               xx+=g3_1[g3_1_len/2-(j+1)]*xin[i-j];
            if((i+1+j) < n) /* right side of filter */
               xx+=g3_1[g3_1_len/2+j]*xin[i+1+j];
        }
        xout[ifac*i+1] = xx;
    } 

    /* here xx = interpolated sample 2/3 between index i and i+1 */
    /*       use g3_1[] in reverse */ 
    for(i = 0; i < n-1; i++)
    {
        xx =0.0;
        for(j = 0; j < g3_1_len/2; j++)
        {
            if((i-j) >= 0) /* left side of filter */
               xx+=g3_1[g3_1_len/2+j]*xin[i-j];
            if((i+1+j) < n) /* right side of filter */
               xx+=g3_1[g3_1_len/2-(j+1)]*xin[i+1+j];
        }
        xout[ifac*i+2] = xx;
     }

}

/* interpolation by factor 5 */
void  ipol5(xin,xout,n)
double  *xin, *xout;
long n;
{
    double xx;
    long i,j, k;
    long ifac = 5; /* interpolation factor */
    
    for(i = 0; i < n-1; i++)
    {
        xout[ifac*i] = xin[i];
        xout[ifac*i+ifac] = xin[i+1];
    }

    /* here xx = interpolated sample 1/5 between index i and i+1 */
    for(i = 0; i < n-1; i++)
    {
        xx =0.0;
        for(j = 0; j < g5_1_len/2; j++)
        {
            if((i-j) >= 0) /* left side of filter */
               xx+=g5_1[g5_1_len/2-(j+1)]*xin[i-j];
            if((i+1+j) < n) /* right side of filter */
               xx+=g5_1[g5_1_len/2+j]*xin[i+1+j];
        }
        xout[ifac*i+1] = xx;
    } 

    /* here xx = interpolated sample 2/5 between index i and i+1 */
    for(i = 0; i < n-1; i++)
    {
        xx =0.0;
        for(j = 0; j < g5_2_len/2; j++)
        {
            if((i-j) >= 0) /* left side of filter */
               xx+=g5_2[g5_2_len/2-(j+1)]*xin[i-j];
            if((i+1+j) < n) /* right side of filter */
               xx+=g5_2[g5_2_len/2+j]*xin[i+1+j];
        }
        xout[ifac*i+2] = xx;
    } 

    /* here xx = interpolated sample 3/5 between index i and i+1 */
    /*       use g5_2[] in reverse */ 
    for(i = 0; i < n-1; i++)
    {
        xx =0.0;
        for(j = 0; j < g5_2_len/2; j++)
        {
            if((i-j) >= 0) /* left side of filter */
               xx+=g5_2[g5_2_len/2+j]*xin[i-j];
            if((i+1+j) < n) /* right side of filter */
               xx+=g5_2[g5_2_len/2-(j+1)]*xin[i+1+j];
        }
        xout[ifac*i+3] = xx;
     }

    /* here xx = interpolated sample 4/5 between index i and i+1 */
    /*       use g5_1[] in reverse */ 
    for(i = 0; i < n-1; i++)
    {
        xx =0.0;
        for(j = 0; j < g5_1_len/2; j++)
        {
            if((i-j) >= 0) /* left side of filter */
               xx+=g5_1[g5_1_len/2+j]*xin[i-j];
            if((i+1+j) < n) /* right side of filter */
               xx+=g5_1[g5_1_len/2-(j+1)]*xin[i+1+j];
        }
        xout[ifac*i+4] = xx;
     }

}



