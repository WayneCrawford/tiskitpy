/**
   NAME: fir2caus
   SYNOPSIS:
       fir2caus <flags>
       flags: -c <filter coefficient file for correction filter>
              -i <data input file>
              -o <data output file>
              -z <no. of zeros to append to input file on both ends> (Optional!)
              -t (shift output trace back by time tag difference. Optional!)
              -f fdig (overwrite fdig [Hz] from protocol file by this value)
                 Use in case same filter stage is used for different data streams!

   VERSION: 2.0
   DATE: 1996-10-18 (Frank Scherbaum)
   DESCRIPTION: This program implements the FIR filter correction described in
   chapter 8 of Scherbaum, F., Of poles and zeros: fundamentals of digital
   seismology, Kluwer Academic Publishers, 1996. The key equation
   is (8.15) on page 120. The ARMA coefficients are assumed to be 
   stored in a '*.prt' file produced by Frank Scherbaum's analyse.m
   Mathematica program. This is the  <filter coefficient file
   for correction filter> above. 

   The seismic recording is assumed to be in the
   <data input file>. The corrected trace is written into the <output file>.
   The data  files are 1 value per row ASCII text files.

   Using a -z <n> option, <n> number of zeros are appended to the data from the
   input file on input on both ends. This option prevents the supression of cutting off data from
   the output file in case the response of the correction filter is longer than the
   input trace (e.g. if a FIR filter response itself is being corrected).
   
   Using the -t option, the output data are time shifted by
   that number of samples corresponding to the time tag found in the
   correction file. Notice!!!! This option may cause 'acausal' oscillations
   with the Nyquist frequency due to the fact that non-integer multiple shifts
   of one sample are involved. For a discussion of this effect see Scherbaum,F.,
   Of poles and zeros, page 189.

**/
/* begin Makefile defines     */
/* define DBL                 for double  precision calculation */
/* define SNG                 for single precision calculation */
/* end  Makefile defines      */
#ifndef DBL
#    define SNG               /* default precision */
#endif

#define PI 3.141592653589793

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
//#include <malloc.h>
#include <memory.h>
main(argc,argv)
        int argc;
        char    *argv[];
{
    FILE *fo;                   /* FILE pointer output file            */
    FILE *fi;                   /* FILE pointer input file             */
    char fir_name[80];          /* input file name FIR filter          */
    char long_fir_name[300];    /* input file name FIR filter including path  */
    char in_name[80];           /* input file name trace file          */
    char out_name[80];          /* output file name                    */
    char *buffer;               /* input buffer                        */
    char *getline_fs();
    int i,ii,j,k,l;             /* indeces                             */
    int nskip = 0;              /* no. of lines to skip on input       */
    char ch;                    /* char. buffer                        */
    int ntot;                   /* total no. of lines in file          */
#ifdef SNG
    float x_inp;                /* input float buffer                  */
#endif
#ifdef DBL
    double x_inp;               /* input float buffer                  */
#endif
    int ndat;                   /* length of wavelet                   */
    int wv_m;                   /* degree of wavelet polynomial        */
    char text_buffer[255];      /* input text buffer                   */
    char *string_end;           /* pointer to possible CR in buffer    */
    float freq;                 /* frequency in Hz                     */
    float fdig;                 /* digitization frequency in Hz        */
    float fdig2;                /* fdig from command line overwites fdig */

    int pad_zeros;              /* no. of padded zeroes for filtering  */
    int lead_zeros;             /* no. of leading zeroes automatically */
                                /* put before input trace              */
    int n,m;                    /* degrees of ARMA coeeficienst        */
    int ndat2;                  /* length of output file               */

#ifdef SNG
    float *x2;                  /* input trace                         */
    float *x;                   /* time reversed input trace           */
    float *y;                   /* filtered output trace               */
    float *a,*b;                /* ARMA coefficients                   */
#endif
#ifdef DBL
    double *x2;                 /* input trace                         */
    double *x;                  /* time reversed input trace           */
    double *y;                  /* filtered output trace               */
    double *a,*b;               /* ARMA coefficients                   */
#endif
    float *temp;                /* buffer for time shift               */
    int max_pos;                /* index of maximum point of FIR filter*/
    float shift_samples;        /* time shift in samples caused by the */
                                /* causal linear phase FIR filter:     */
                                /* = (float)(ndat-1)/2.0               */
    int correct_time;           /* if 1 -> time correction for linear  */
                                /* phase performed, else not           */
    char *envptr;


    int no_ar,no_ma;
    int mx;
    char decbuff[50];
    int cnt;
    float time_tag;

    pad_zeros =  0;
    fdig2 = -999; /* is only used if -f flag given */
    correct_time = 0;
    for (i=1; i<argc;i++)
    {
        if (argv[i][0] == '-' && strlen(argv[i]) >= 2 )
        {
            switch( argv[i][1])
            {
                case 'c':  /* correction filter coefficient */
                case 'C':
                {
                    i++;
                    strcpy(fir_name,&argv[i][0]);
                    break;
                }
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
                case 'z':  /* no. of zeros to append on input of data file */
                case 'Z':
                {
                    i++;
                    pad_zeros = atoi(&argv[i][0]);
                    break;
                }
                case 't':  /* correct for linear phase shift */
                case 'T':
                {
                    correct_time = 1;
                    break;
                }
                case 'f':  /* overwrite fdig by command line value */
                case 'F':
                {
                    i++;
                    fdig2 = atof(&argv[i][0]);
                    break;
                }
            }
        }
    }

    if ((argc < 3) ||
        (strlen(fir_name) == 0)||
        (strlen(in_name) == 0)||
        (strlen(out_name) == 0))
    {
       printf("fir2caus Version 2.0 (1996/10/19)\n");
       printf("USAGE: fir2caus <flags>\n");
       printf("flags: -c <filter coefficients correction filter>\n");
       printf("       -i <data input file>\n");
       printf("       -o <output file>\n");
       printf("       -f fdig (overwrite fdig [Hz] from protocol file by this value)\n");
       printf("          Use in case same filter stage is used for different data streams!\n");
       printf("       -z <no. of zeros to append to input file on both ends> (Optional!)\n");
       printf("       -t (shift output trace back by time tag difference (Optional!)\n");
       exit(1);
    }

    if ((envptr=getenv("FIR_CORR_COEFF_PATH")) == NULL)
    {
        fprintf(stderr,"fir2caus: env.var. FIR_CORR_COEFF_PATH not set: can't continue\n");
        exit(1);
    }
    else
    {
        strcpy(long_fir_name,envptr);
    }
    strcat(long_fir_name,fir_name);

    printf("Correction coefficient file: %s\n",long_fir_name);

    /* read correction filter coefficients */
    if ((fi = fopen(long_fir_name,"rt")) != NULL)
    {
      /* find FDIG */
      do
      {
          if ((buffer = getline_fs(fi)) == NULL)
          {
             printf("no FDIG information found!\n");
             fclose(fi);
             exit(1);
          }
          free((char *)buffer);
      } while (strncmp(buffer,"#FDIG_EFF",9) != 0);

      if ((buffer = getline_fs(fi)) == NULL)
      {
         printf("Error on reading FDIG!\n");
         fclose(fi);
         exit(1);
      } else {
        sscanf(buffer," %f ",&fdig);
      }
      free((char *)buffer);
      /* Overwrite fdig in case fdig2 has been defined */
      if (fdig2 > 0) fdig = fdig2;

      /* find AR coefficients */
      do
      {
          if ((buffer = getline_fs(fi)) == NULL)
          {
             printf("no AR coefficients found!\n");
             fclose(fi);
             exit(1);
          }
          free((char *)buffer);
      } while (strncmp(buffer,"#CORR_AR",8) != 0);

      if ((buffer = getline_fs(fi)) == NULL)
      {
         printf("Error on reading no of AR coefficients!\n");
         fclose(fi);
         exit(1);
      } else {
        sscanf(buffer," %d ",&no_ar);
      }
      free((char *)buffer);

#ifdef SNG
      a = (float *)calloc(no_ar+1,sizeof(float));
#endif
#ifdef DBL
      a = (double *)calloc(no_ar+1,sizeof(double));
#endif

      do
      {
          if ((buffer = getline_fs(fi)) == NULL)
          {
             printf("no AR coefficients found!\n");
             fclose(fi);
             exit(1);
          }
          free((char *)buffer);
      } while (strncmp(buffer,"#CORR_AR_DATA",13) != 0);

      ii=1; /* the first AR coefficient in the coefficient file is a[1] */
            /* a[1] = 1 anyway and is not used. See equ. 8.15 in */
            /* Of poles and zeros. */
      while(ii <= no_ar )
      {
          if ((buffer = getline_fs(fi)) == NULL)
          {
             printf("Error on reading AR coefficients!\n");
             fclose(fi);
             exit(1);
          } else {  /* read this line */
              cnt = 0;
              while(cnt < strlen(buffer)) {
                  /* skip initial whitespaces */
                  while (isspace(buffer[cnt]) && ( cnt < strlen(buffer) ))
                      cnt++;
                  i=0;
                  while(!isspace(buffer[cnt]) && ( cnt < strlen(buffer))){
                      decbuff[i++] = buffer[cnt++];
                  }
                  decbuff[i] = '\0';
#ifdef SNG
                  sscanf(decbuff," %f ",&a[ii++]);
#endif
#ifdef DBL
                  sscanf(decbuff," %lf ",&a[ii++]);
#endif
              }
          }
          free((char *)buffer);
      }
      /* find MA coefficients */
      do
      {
          if ((buffer = getline_fs(fi)) == NULL)
          {
             printf("no MA coefficients found!\n");
             fclose(fi);
             exit(1);
          }
          free((char *)buffer);
      } while (strncmp(buffer,"#CORR_MA",8) != 0);
      if ((buffer = getline_fs(fi)) == NULL)
      {
         printf("Error on reading no of MA coefficients!\n");
         fclose(fi);
         exit(1);
      } else {
        sscanf(buffer," %d ",&no_ma);
      }
      free((char *)buffer);
#ifdef SNG
      b = (float *)calloc(no_ma+1,sizeof(float));
#endif
#ifdef DBL
      b = (double *)calloc(no_ma+1,sizeof(double));
#endif

      do
      {
          if ((buffer = getline_fs(fi)) == NULL)
          {
             printf("no AR coefficients found!\n");
             fclose(fi);
             exit(1);
          }
          free((char *)buffer);
      } while (strncmp(buffer,"#CORR_MA_DATA",13) != 0);

      ii=0; /* the first MA coefficient in coefficient file is b[0] */
      while(ii <= no_ma )
      {
          if ((buffer = getline_fs(fi)) == NULL)
          {
             printf("Error on reading AR coefficients!\n");
             fclose(fi);
             exit(1);
          } else {  /* read this line */
              cnt = 0;
              while(cnt < strlen(buffer)) {
                  /* skip initial whitespaces */
                  while(isspace(buffer[cnt]) && ( cnt < strlen(buffer)))
                      cnt++;
                  i=0;
                  while(!isspace(buffer[cnt]) && ( cnt < strlen(buffer))){
                      decbuff[i++] = buffer[cnt++];
                  }
                  decbuff[i] = '\0';
#ifdef SNG
                  sscanf(decbuff," %f ",&b[ii++]);
#endif
#ifdef DBL
                  sscanf(decbuff," %lf ",&b[ii++]);
#endif
              }
          }
          free((char *)buffer);
      }

      /* find TIMETAG correction */
      do
      {
          if ((buffer = getline_fs(fi)) == NULL)
          {
             printf("no time tage correction found!\n");
             fclose(fi);
             exit(1);
          }
          free((char *)buffer);
      } while (strncmp(buffer,"#TIMETAG",8) != 0);

      if ((buffer = getline_fs(fi)) == NULL)
      {
         printf("Error on reading no of AR coefficients!\n");
         fclose(fi);
         exit(1);
      } else {
        sscanf(buffer," %f ",&time_tag);
        shift_samples = -time_tag;
      }
      free((char *)buffer);

    }  else {
       printf("Error on opening file %s\n",fir_name);
       exit(1);
    }
    fclose(fi);

    /* INPUT: trace to be filtered   */
    /* estimate no. of lines in file */
    ntot  = 0;
    if ((fo = fopen(in_name,"rt")) != NULL)
    {
            while((ch=fgetc(fo)) != EOF)
            {
                if (ch == '\n')
                    ntot++;
            }
        fclose(fo);
    }

    lead_zeros = 0; /* use as default */
    if(pad_zeros > 0) lead_zeros = pad_zeros;
    ndat2 = ntot - nskip + lead_zeros + pad_zeros;

#ifdef SNG
    x  = (float *)calloc(ndat2+1,sizeof(float));
    x2 = (float *)calloc(ndat2+1,sizeof(float));
#endif
#ifdef DBL
    x  = (double *)calloc(ndat2+1,sizeof(double));
    x2 = (double *)calloc(ndat2+1,sizeof(double));
#endif

    /* INPUT data file */
    if((fi = fopen(in_name,"rt")) == NULL)
    {
            printf("input file %s cannot be opened!\n",in_name);
            exit(1);
    }
    for(k=0;k < ndat2-nskip-pad_zeros-lead_zeros;k++)
    {
        if(fgets(text_buffer,255,fi) != NULL)
        {
            if ((string_end = memchr(&text_buffer[0],'\n',255)) != NULL)
                    *string_end = '\0';
#ifdef SNG
            sscanf(text_buffer,"%f",&x_inp);
#endif
#ifdef DBL
            sscanf(text_buffer,"%lf",&x_inp);
#endif
        }
        else
            x_inp = 0.0;

        if(k >= nskip)
            x2[k-nskip+lead_zeros] = x_inp;
    }
    fclose(fi);

    /* flip in time */
    for(k=0;k<ndat2;k++)
    {
        x[k]=x2[ndat2-1-k];
    }

    /* output trace                        */
#ifdef SNG
    y = (float *)calloc(ndat2+1,sizeof(float));
#endif
#ifdef DBL
    y = (double *)calloc(ndat2+1,sizeof(double));
#endif


/*

    Filter x[] which is the reversed input sequence x2[] 
    using the difference equation:

            mx                mx
            --                -----
    y'[i] =  > a[k]*y'[i-k]   + > b[l] x[i-k]
            __                __
            k=1               l=0

    This corresponds to equ. (8.15) in 'Scherbaum, F: Of poles and zeros,
    Fundamentals of Digital Seismology, Kluwer Academic Publ., 1996' 
    mx = number of AR coefficients
    b[l] = MA coefficients for l = 0, mx
    a[k] = AR coefficients for k = 1, mx

    x[i] = reversed input sequence for i = 0 ..... 
    y'[] = output sequence 

    Reverse the output sequence y'[] in time again to obtain the 
    corrected sequence y[n]!

*/

    mx = no_ar;

    /* filter                              */
    for (i=0; i< ndat2; i++) {
        y[i] = 0.0;
        /* MA */
        for (l=0; l <= mx; l++) {
            if ((i-l) >= 0) {
                y[i] += x[i-l]*b[l];
            }
        }
        /* AR */
        for (k=1; k<= mx; k++) {
            if ((i-k) >= 0) {
                y[i] += y[i-k]*a[k];
            }
        }
    }

    temp = (float *)calloc(ndat2+1,sizeof(float));
    /* flip back in time */

    for(k=0;k<ndat2;k++)
    {
        temp[k]=(float)y[ndat2-1-k];
    }

    if (correct_time == 1)
    {
        time_shift(temp,ndat2,-shift_samples);

        printf("Trace physically delayed by: %f [samples] = %f [sec]\nNo signal front advance in output trace.",-shift_samples,-shift_samples/fdig);

        /*
        printf("Trace physically delayed by: %f [samples]\nNo signal front advance in output trace.",-shift_samples);
        */
    }
    else

        printf("Signal front advance in output trace: %f [samples] at %f [Hz] += %f [sec]\n",-shift_samples,fdig,-shift_samples/fdig);
        /*
        printf("Signal front advance in output trace: %f [samples]\n",-shift_samples);
        */

    /* output                              */
    if((fo = fopen(out_name,"wt")) == NULL)
    {
            printf("output file %s cannot be opened!\n",out_name);
            exit(1);
    }
    if(pad_zeros == 0)
    {
        /* write out trace without padding zeros */
        for(j=0;j<ndat2-lead_zeros;j++)
        {
                fprintf(fo,"%g\n",temp[j+lead_zeros]);

        }
    } else
    {
        /* write out trace with padding zeros */
        for(j=0;j<ndat2;j++)
        {
                fprintf(fo,"%g\n",temp[j]);

        }
    }

    fclose(fo);

    /* free allocated memory               */
    free((char *)a);
    free((char *)b);
    free((char *)x);
    free((char *)x2);
    free((char *)y);
    free((char *)temp);

}
/**
   NAME: time_shift
   SYNOPSIS:
   float *y;
   int ndat;
   float shift_samples;
   time_shift(y,ndat,extra_samples);
   DESCRIPTION: Performs a time shift in the frequency domain by
   multiplication with the corresponding phase shift operator.
   DATE: July 2, 1993 (Frank Scherbaum)
**/
int time_shift(y,ndat,shift_samples)
float *y;
int ndat;
float shift_samples;
{
   float *shift;           /* spectrum corresponding to the desired time shift */
   int i,j;
   float amp,phase;        /* amplitude and phase of shifting spectrum */
   float real, imag;       /* real, imaginary of shifting spectrum */
   float t_samp;           /* sampling interval */
   int ia;                 /* integer part of shift */
   float extra_samples;    /* fraction of samples to shift */
   float *b1;              /* trace buffers */
   float *bb1, *bb2;       /* trace buffers */
   float x1;               /* dummy variable */
   int nfft;
   /*  allocate buffer */
   b1 = (float *)calloc(ndat,sizeof(float));
   ia = 0;
   if(shift_samples == (int)shift_samples)
   { /*integer multiple of 1 sample */
       ia = (int)shift_samples;
       extra_samples = 0.0;
   }
   else
   {
       ia = (int)shift_samples;
       extra_samples = shift_samples - ia;
   }
   if(ia >0)
   {
       for (j = ndat-1 ; j >= ia; j--)
       {
           x1 = *(y + j - ia);
           *(b1+j) = x1;
       }
       for (j = 0; j < ia; j++)
           *(b1+j) = 0.0;
   }
   else if(ia < 0)
   {
       for (j = 0 ; j < ndat +ia -1 ; j++)
       {
           x1 = *(y + j - ia);
           *(b1+j) = x1;
       }
       for (j = ndat-ia ; j <ndat; j++)
       {
           *(b1+j) = 0.0;
       }
   }
   else if (ia == 0)
   {
       for (j = 0; j < ndat; j++)
           *(b1+j) = *(y + j);
   }
   /*
      every shift less than a sample is done in the
      frequency domain
   */
   if (extra_samples != 0.0) /* non-integer part of shift */
   {
       nfft = 1;
       while(nfft < ndat)
           nfft *= 2;
       nfft *= 2; /* to avoid wrap around */
       bb1 = (float *)calloc(nfft,sizeof(float));
       bb2 = (float *)calloc(nfft,sizeof(float));
       for(j=0;j<ndat;j++)
       {
           bb1[j] = b1[j];
       }
       realft(bb1-1,(long)nfft,1);
       amp = 1.0;
       /*
       1 sample shift ==
       phase of PI at f Nyquist (linear in between)
       */
       bb2[0] = amp;
       bb2[1] = amp;
       for (j=1; j < nfft/2; j++)
       {
           phase = extra_samples*
                   ((float)j/(float)nfft)*PI;
           real = amp*cos(phase);
           imag = amp*sin(phase);
           bb2[2*j] = real;
           bb2[2*j+1] = imag;
       }
       /* spectral multiplication */
       bb1[0] *= bb2[0];
       bb1[1] *= bb2[1];
       for (j=1; j < nfft/2; j++)
       {
           real = bb1[2*j]*bb2[2*j] - bb1[2*j+1]*bb2[2*j+1];
           imag = bb1[2*j]*bb2[2*j+1] + bb1[2*j+1]*bb2[2*j];
           bb1[2*j] = real;
           bb1[2*j+1] = imag;
       }
       /* inverse FFT */
       realft(bb1-1,(long)nfft,-1);
       /* scale amplitudes back */
       for(j=0;j<ndat;j++)
       {
           y[j] = 2*bb1[j]/nfft;
       }
       free((char *)b1);
       free((char *)bb1);
       free((char *)bb2);
   } else {
       for(j=0;j<ndat;j++)
       {
           y[j] = b1[j];
       }
   };
}
/* getline_fs.c -- Read a line */
/* and dynamically allocate memory for it */
/* Taken from BYTE*June 1988 p.313-318:*/
/* Dynamic memory management in C, by */
/* David L. Fox */
/* CHANGED: NOV 2, 1992 since it did not replace \n by \0 properly */
/* Frank Scherbaum */

/* Return pointer to string or NULL on error */
#define MEMINCR 256 /* increment in size of memory block */

char *getline_fs(infile)
FILE *infile;
{
    char *r;
    size_t n,m;
    int c;

    n = 0;  /* # of bytes read */
    m = MEMINCR;    /* available space */
    r = malloc(m+1); /* allow room for \0 */
    do {
        if (--m == 0) {
            if((r = realloc(r, n+MEMINCR+1)) == NULL) {
                return NULL;
            }
            m = MEMINCR;
        }
        if ((c = getc(infile)) == EOF) {
            free(r);
            return NULL;
        }
        r[n++] = c;
    } while (c != '\n');
    r[n-1] = '\0';    /* terminate the string */
    if ((r = realloc(r,n)) == NULL) {
        return NULL;
    }
    return r;
}
