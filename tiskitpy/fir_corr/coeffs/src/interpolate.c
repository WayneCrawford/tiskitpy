#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
//#include <malloc.h>
#include <memory.h>

/**
   NAME: interpolate
   SYNOPSIS:

       interpolate

   VERSION: 1.0
   DATE: 1997-01-31
   DESCRIPTION: Change sampling rate of an ASCII input file by interpolation.
   Interpolation operators by E. Wielandt. 

**/
main(argc,argv)
        int argc;
        char    *argv[];

{
    long i,k;

    FILE *fo;
    FILE *fi;
    char in_name[80];
    char out_name[80];
    long nskip;
    int conv_fac,ipol_fac;
    int no_it2,no_it3,no_it5;     /* no of iterations */
    int pow2,pow3,pow5;          /* powers of 2,3,5 in interpolation factor */
    char ch;
    long ntot;
    long ndat;
    double *tr_in, *tr_out;
    double x;
    char text_buffer[255];
    char *string_end;    /* pointer to possible CR */

    void ipol2();
    void ipol3();   
    void ipol5();

    nskip = 0L;
    if (argc < 5)
    {
       printf("Interpolation based on E. Wielandt's optimum FIR filters (Wielandt, 1997)\n");
       printf("USAGE: interpolate  -i <input file> -o <output file> -f <integer interpolation factor>\n");
       printf("                   [-s <no. pts to skip (0)>]\n");
       exit(1);
    }
    conv_fac = 2;
    for (i=1; i<argc;i++)
    {
        if (argv[i][0] == '-' && strlen(argv[i]) >= 2 )
        {
            switch( argv[i][1])
            {
                case 'i':  /* input file name */
                case 'I':  /* input file name */
                {
                    i++;
                    strcpy(in_name,&argv[i][0]);
                    break;
                }
                case 'o':  /* output file name */
                case 'O':  /* output file name */
                {
                    i++;
                    strcpy(out_name,&argv[i][0]);
                    break;
                }
                case 's':  /* number of lines to skip on input */
                case 'S':
                {
                    i++;
                    nskip = atol(&argv[i][0]);
                    break;
                }
                case 'f':   /* interpolation factor */
                case 'F':
                {
                    i++;
                    ipol_fac = atoi(&argv[i][0]);
                    break;
                }
            }
        }
    }

    /* estimate no. of lines in file */
    ntot = 0L;
    if ((fi = fopen(in_name,"rt")) != NULL)
    {
        while((ch=fgetc(fi)) != EOF)
        {
            if (ch == '\n')
                ntot++;
        }
        fclose(fi);
    }
    ndat = ntot - nskip;
    if(ndat <2)
    {
        printf("ERROR: trace too short\n");
        exit(1);
    }

    if( ipol_fac == 1)
    {
        printf("Nothing to interpolate, exit...\n");
        exit(1);
    }

    conv_fac = ipol_fac;
    /*  number of divisions by 2 */
    pow2 = 2;
    no_it2 = 0;
    while(((float)conv_fac/(float)pow2) == (float)(conv_fac/pow2))
    { 
        pow2 *=2;
        no_it2 +=1;
    }
    if(no_it2 >0)
        conv_fac = 2*conv_fac/pow2;
    /*  number of divisions by 3 */
    pow3 = 3;
    no_it3 = 0;
    while(((float)conv_fac/(float)pow3) == (float)(conv_fac/pow3))
    {
        pow3 *=3;
        no_it3 +=1;
    }
    if(no_it3 > 0)
        conv_fac = 3*conv_fac/pow3;
    /*  number of divisions by 5 */
    pow5 = 5;
    no_it5 = 0;
    while(((float)conv_fac/(float)pow5) == (float)(conv_fac/pow5))
    {
        pow5 *=5;
        no_it5 +=1;
    }
    if (no_it5 >0)
        conv_fac = 5*conv_fac/pow5; 
   
    if( conv_fac != 1)
    {
        printf("Factor %d cannot be separated in factors 2, 3, and 5\n",ipol_fac);
        exit(1);
    }

    tr_in = (double *)calloc(ndat,sizeof(double));
    /* INPUT  */
    if((fi = fopen(in_name,"rt")) == NULL)
    {
            printf("input file %s cannot be opened!\n",in_name);
            exit(1);
    }

    for(k=0;k<ndat+nskip;k++)
    {
        if(fgets(text_buffer,255,fi) != NULL)
        {
            if ((string_end = memchr(&text_buffer[0],'\n',255)) != NULL)
                    *string_end = '\0';
            sscanf(text_buffer,"%lf",&x);
        }
        else
            x = 0.0;
        if(k >= nskip)
            *(tr_in+k-nskip) = x;
    }
    fclose(fi);

    /* interpolations by factors of  2 */
    for(k=0;k<no_it2;k++)
    {
       tr_out = (double *)calloc(ndat*2,sizeof(double));      
       ipol2(tr_in,tr_out,ndat);
       ndat *= 2;
       free((char *) tr_in);
       tr_in = tr_out;
    }
    /* interpolations by factors of 3 */
    for(k=0;k<no_it3;k++)
    {
       tr_out = (double *)calloc(ndat*3,sizeof(double));      
       ipol3(tr_in,tr_out,ndat);
       ndat *= 3;
       free((char *) tr_in);
       tr_in = tr_out;
    }
    /* interpolations by factors of  5 */
    for(k=0;k<no_it5;k++)
    {
       tr_out = (double *)calloc(ndat*5,sizeof(double));      
       ipol5(tr_in,tr_out,ndat);
       ndat *= 5;
       free((char *) tr_in);
       tr_in = tr_out;
    }

    /* OUTPUT  */
    if((fo = fopen(out_name,"wt")) == NULL)
    {
            printf("output file %s cannot be opened!\n",out_name);
            exit(1);
    } 
    for(k=0;k<ndat;k++)
    {
        fprintf(fo,"%g\n",*(tr_in+k));

    }
    fclose(fo);
    free((char *)tr_in);

}

