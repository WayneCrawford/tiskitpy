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
   NAME: decimate
   SYNOPSIS:

       decimate

   VERSION: 1.0
   DATE: 1996-10-25
   DESCRIPTION: Decimate the  ASCII input file in the time domain.


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
    long conv_fac;   /* decimation ratio */
    long start_dec;  /* start index for decimation */

    char ch;
    long ntot;
    long ndat,ndat2;
    double *tr;
    double x;
    char text_buffer[255];
    char *string_end;    /* pointer to possible CR */
    void decim();

    nskip = 0L;
    if (argc < 5)
    {
       printf("USAGE: interpol  -i <input file> -o <output file> -d <decimation factor>\n");
       printf("                 -f <force decimation to include this sample index; first index := 0>\n ");
       printf("                [-s <no. pts to skip (0)>]\n ");
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
                case 'd':   /* decimation factor */
                case 'D':
                {
                    i++;
                    conv_fac = atol(&argv[i][0]);
                    break;
                }
                case 'f':   /* force decimation to include this index */
                case 'F':
                {
                    i++;
                    start_dec = atol(&argv[i][0]);
                    break;
                }
            }
        }
    }

    /* estimate no. of lines in file */
    ntot = 0L;
    if ((fo = fopen(in_name,"rt")) != NULL)
    {
            while((ch=fgetc(fo)) != EOF)
            {
                if (ch == '\n')
                    ntot++;
            }
        fclose(fo);
    }
    ndat = ntot - nskip;

    tr = (double *)calloc(ndat,sizeof(double));

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
            *(tr+k-nskip) = x;
    }
    fclose(fi);
    
    decim(tr,&ndat,conv_fac,start_dec);

    /* OUTPUT  */
    if((fo = fopen(out_name,"wt")) == NULL)
    {
            printf("output file %s cannot be opened!\n",out_name);
            exit(1);
    }
  
    for(k=0;k<ndat;k++)
    {
        fprintf(fo,"%g\n",*(tr+k));

    }
    fclose(fo);
    free((char *)tr);

}

