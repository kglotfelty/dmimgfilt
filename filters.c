/*                                                                
**  Copyright (C) 2004-2008  Smithsonian Astrophysical Observatory 
*/                                                                

/*                                                                          */
/*  This program is free software; you can redistribute it and/or modify    */
/*  it under the terms of the GNU General Public License as published by    */
/*  the Free Software Foundation; either version 3 of the License, or       */
/*  (at your option) any later version.                                     */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful,         */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*  GNU General Public License for more details.                            */
/*                                                                          */
/*  You should have received a copy of the GNU General Public License along */
/*  with this program; if not, write to the Free Software Foundation, Inc., */
/*  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.             */
/*                                                                          */

#include <dsnan.h>
#include <dmimgfilt.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>  /* for memset */
#include <time.h>

#ifndef MIN
#define MIN(a,b) (((a)<=(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b) (((a)>=(b))?(a):(b))
#endif


double minpixel;
double maxpixel;
double pixelrange;

double *kuw_vals[4];
long kuw_num[4];

long *histogram;
long num_hist_bins;

double centerval;



#define SIGMA_ITER 5



double _filtCOUNT( double *vals, long nvals ) {
  return(nvals);
}


double _filtSUM(double *vals, long nvals ) { 
  double retval;
  long ii;

  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }

  retval = vals[0];
  for(ii=nvals;--ii;)
    retval += vals[ii];
  return(retval); 
}


double _filtMIN(double *vals, long nvals ) { 
  double retval;
  long ii;


  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }


  retval = vals[0];
  for(ii=nvals;--ii;)
    retval = MIN(retval, vals[ii] );
  return(retval); 
}


double _filtMAX(double *vals, long nvals ) {
  double retval;
  long ii;
  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }

  retval = vals[0];
  for(ii=nvals;--ii;)
    retval = MAX(retval, vals[ii] );
  return(retval); 
}

double _filtEXTREME(double *vals, long nvals ) {
  double retval;
  long ii;
  double min;
  double max;
  double mean;
  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }

  min = vals[0];
  max = vals[0];
  mean = vals[0];
  for(ii=nvals;--ii;) {
    mean += vals[ii];
    min = MIN( min, vals[ii]);
    max = MAX( max, vals[ii]);
  }
  mean /= nvals;
  retval = ( (mean-min) < (max-mean) ) ? min : max;
  return(retval);
}

double _filtRANGE(double *vals, long nvals ) {
  double retval;
  long ii;
  double min;
  double max;

  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }
  min = vals[0];
  max = vals[0];

  for(ii=nvals;--ii;) {
    min = MIN( min, vals[ii]);
    max = MAX( max, vals[ii]);
  }
  retval = (max - min );
  return(retval);
}

double _filtMID(double *vals, long nvals ) {
  double retval;
  long ii;
  double min;
  double max;

  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }
  min = vals[0];
  max = vals[0];

  for(ii=nvals;--ii;) {
    min = MIN( min, vals[ii]);
    max = MAX( max, vals[ii]);
  }
  retval = (max + min )/2.0;
  return(retval);
}

double _filtRCLIP(double *vals, long nvals ) {
  double retval;
  long ii;
  double min;
  double max;
  short center=0;

  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }
  min = vals[0];  /* Assume center value not 1st value */
  max = vals[0];

  for(ii=nvals;ii--;) {
    if (( vals[ii] == centerval ) && (center == 0 )) {
      center=1; /* want the min/max excludeding the center value */
      continue;
    }
    min = MIN( min, vals[ii]);
    max = MAX( max, vals[ii]);
  }
  if ( centerval < min ) 
    retval = min;
  else if ( centerval > max )
    retval = max;
  else
    retval = centerval;

  return(retval);
}


double _filtLOCEQ(double *vals, long nvals ) {
  double retval;
  long ii;
  double min;
  double max;
  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }
  min = vals[0];
  max = vals[0];


  for(ii=nvals;--ii;) {
    min = MIN( min, vals[ii]);
    max = MAX( max, vals[ii]);
  }
  if ( max == min ) 
    retval = (min-minpixel)*pixelrange + minpixel;
  else
    retval = ((centerval - min)*pixelrange/(max - min)) + minpixel;
  return(retval);
}

double _filtMEAN(double *vals, long nvals ) {
  double retval;
  long ii;
  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }
  retval=vals[0];
  for (ii=nvals;--ii;)
    retval += vals[ii];
  retval /= nvals; 
  return(retval); 
}


double kth_smallest(double *a, long n, long k);


double _filtMEDIAN(double *vals, long nvals ) {
  double retval;
  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }

  if ( 1==nvals )
    retval = vals[0];
  else
    retval = kth_smallest(vals,nvals,(((nvals)&1)?((nvals)/2):(((nvals)/2)-1))); 
  return(retval); 
}



/* http://ndevilla.free.fr/median/median/node20.html
 *
 * Algorithm from N. Wirth's book, implementation by N. Devillard.
 * This code in public domain.
 */


#define ELEM_SWAP(a,b) { register double t=(a);(a)=(b);(b)=t; }


/*---------------------------------------------------------------------------
   Function :   kth_smallest()
   In       :   array of elements, # of elements in the array, rank k
   Out      :   one element
   Job      :   find the kth smallest element in the array
   Notice   :   use the median() macro defined below to get the median. 

                Reference:

                  Author: Wirth, Niklaus 
                   Title: Algorithms + data structures = programs 
               Publisher: Englewood Cliffs: Prentice-Hall, 1976 
    Physical description: 366 p. 
                  Series: Prentice-Hall Series in Automatic Computation 

 ---------------------------------------------------------------------------*/


double kth_smallest(double *a, long n, long k)
{
    register long i,j,l,m ;
    register double x ;

    l=0 ; m=n-1 ;
    while (l<m) {
        x=a[k] ;
        i=l ;
        j=m ;
        do {
            while (a[i]<x) i++ ;
            while (x<a[j]) j-- ;
            if (i<=j) {
                ELEM_SWAP(a[i],a[j]) ;
                i++ ; j-- ;
            }
        } while (i<=j) ;
        if (j<k) l=i ;
        if (k<i) m=j ;
    }
    return a[k] ;
}



double _filtMODE(double *vals, long nvals ) {
  double retval;
  long ii;
  double mean;
  double median;
  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }
  mean = vals[0];
  for (ii=nvals;--ii;)
    mean += vals[ii];
  mean /= nvals; 
  median = _filtMEDIAN( vals, nvals );
  retval = 3 * median - 2 * mean;
  return(retval);
}

double _filtNORM_MODE(double *vals, long nvals ) {
  double retval;
  long ii;
  double mean;
  double median;
  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }
  mean = vals[0];
  for (ii=nvals;--ii;)
    mean += vals[ii];
  mean /= nvals; 
  median = _filtMEDIAN( vals, nvals );
  retval = 3 * median - 2 * mean;
  retval /= mean;
  return(retval);
}

double _filtSIG(double *vals, long nvals ) {
  double retval;
  long ii;
  double mean;
  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }
  mean=vals[0];
  for (ii=nvals;--ii;)
    mean += vals[ii];
  mean /= nvals; 
  retval = 0;
  for (ii=nvals;ii--;) {
    double diff = (vals[ii]-mean);
    retval += (diff*diff);
  }
  retval = sqrt( retval / (nvals-1) ) ;
  return(retval);
}

double _filtKUWAHARA(double *vals, long nvals ) {
  double retval;
  long ii;
  long jj;
  double mean[4];
  double sig[4];

  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }
  if ( nvals > 1 ) {
    memset( mean, 0, 4*sizeof(double));
    memset( sig, 0, 4*sizeof(double));
    
    for (ii=4;ii--;) {
      for (jj=kuw_num[ii];jj--; ) {
	mean[ii] += kuw_vals[ii][jj];
      }
      mean[ii] /= kuw_num[ii];
    }
    for (ii=4;ii--;) {
      for (jj=kuw_num[ii];jj--; ) {
	double diff = (kuw_vals[ii][jj]-mean[ii]);
	sig[ii] += (diff*diff);
      }
      sig[ii] = sqrt( sig[ii] / ( kuw_num[ii] -1 ));
    }
    
    jj = 0;
    retval = sig[0];
    if (( sig[1] < retval ) && ((sig[1]-retval)>2*DBL_EPSILON*retval) ) {
      jj = 1; 
      retval=sig[1]; 
    }
    if (( sig[2] < retval ) && ((sig[2]-retval)>2*DBL_EPSILON*retval) ) {
      jj = 2; 
      retval=sig[2]; 
    }
    if (( sig[3] < retval ) && ((sig[3]-retval)>2*DBL_EPSILON*retval) ) {
      jj = 3; 
      retval=sig[3]; 
    }

    retval=mean[jj];
  } else {
    retval=vals[0];
  }
  return(retval);
}

double _filtUNSHARP(double *vals, long nvals ) {
  double retval;
  long ii;
  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }
  retval=vals[0];
  for (ii=nvals;--ii;)
    retval += vals[ii];
  retval /= nvals; 
  retval = centerval - retval;
  return(retval); 
}

double _filtVARIANCE(double *vals, long nvals ) {
  double retval;
  long ii;
  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }
  retval=0;
  for (ii=nvals;ii--;) {
    double diff=(vals[ii]-centerval);
    retval += (diff*diff);
  }
  retval=sqrt(retval);
  return(retval); 
}


double _filtQUANTILE_XX(double *vals, long nvals ) {
  double retval;
  long idx;
  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }
  
  char *xx_str = getenv("QXX");
  if (strlen( xx_str ) == 0 ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }
  

   
  char *end = NULL;
  float ff = strtof( xx_str, &end );

  
  if (( ff > 1) || (ff < 0 ) ) {
    ds_MAKE_DNAN(retval);
    return(retval);      
  }
  
  retval = nvals * ff ;
  idx = retval; /* integer case */
  
  if ( 1 == nvals ) 
    retval = vals[0];
  else
    retval = kth_smallest( vals, nvals, idx );
  return(retval); 
}



double _filtQUANTILE_25(double *vals, long nvals ) {
  double retval;
  long idx;
  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }
  
  retval = nvals / 4.0 ;
  idx = retval; /* integer case */
  
  if ( 1 == nvals ) 
    retval = vals[0];
  else
    retval = kth_smallest( vals, nvals, idx );
  return(retval); 
}

double _filtQUANTILE_33(double *vals, long nvals ) {
  double retval;
  long idx;
  
  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }
  retval = nvals / 3.0 ;
  idx = retval; /* integer case */
  
  if ( 1 == nvals ) 
    retval = vals[0];
  else
    retval = kth_smallest( vals, nvals, idx );
  return(retval); 
}
double _filtQUANTILE_67(double *vals, long nvals ) {
  double retval;
  long idx;

  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }
  retval = 2.0 * nvals / 3.0 ;
  idx = retval; /* integer case */

  if ( 1 == nvals ) 
    retval = vals[0];
  else
    retval = kth_smallest( vals, nvals, idx );
  
  return(retval);
}


double _filtQUANTILE_75(double *vals, long nvals ) {
  double retval;
  long idx;
  
  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }
  retval = 3.0 * nvals / 4.0 ;
  idx = retval; /* integer case */
  if ( 1 == nvals ) 
    retval = vals[0];
  else
    retval = kth_smallest( vals, nvals, idx );

  return(retval);
}

double _filtMOST_COMMON(double *vals, long nvals ) {
  double retval;
  long ii;
  double delta;
  double min, max;

  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }

  min = vals[0];
  max = vals[0];

  if ( nvals>1 ) {
    for(ii=nvals;--ii;) {
      min = MIN( min, vals[ii]);
      max = MAX( max, vals[ii]);
    }
    
    delta = (max - min ) / num_hist_bins;
    
    if ( delta == 0 ) {
      retval=vals[0];
      return(retval);
    }
    
    memset( histogram, 0, num_hist_bins * sizeof(long));
    for(ii=nvals;ii--;) {
      long bin;
      bin = (vals[ii] - min ) / delta;
      if ( bin >= num_hist_bins ) bin = (num_hist_bins -1);
      histogram[bin] +=1 ;
    }
    max = 0;
    retval = 0;
    for (ii=num_hist_bins;ii--; ) {
      if ( histogram[ii] > max ) {
	max = histogram[ii];
	retval = ii;
      }
    }
    /*     printf("num_hist_bin=%d\tretval=%f\n", num_hist_bins, retval); */
    retval = (retval+0.5)*delta + min;
  } else {
    retval = vals[0];
  }
  return(retval);
}


double _filtPEAK(double *vals, long nvals ) { 

  double retval;
  short center;
  long ii;

  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }


  retval = centerval; /* assum peak */
  center = 0;
  for(ii=nvals;ii--;) {
    if (( vals[ii] == centerval ) && (center == 0 )) {
      center=1; /* want the min/max excludeding the center value */
      continue;
    }
    if ( vals[ii] >= centerval ) {
      ds_MAKE_DNAN(retval);
      break;
    }
  }

  return(retval);

}


double _filtRIDGE(double *vals, long nvals ) { 
  /* Same as above but now > instead of >= */
  double retval;
  short center;
  long ii;

  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }


  retval = centerval; /* assum peak */
  center = 0;
  for(ii=nvals;ii--;) {
    if (( vals[ii] == centerval ) && (center == 0 )) {
      center=1; /* want the min/max excludeding the center value */
      continue;
    }
    if ( vals[ii] > centerval ) {
      ds_MAKE_DNAN(retval);
      break;
    }
  }

  return(retval);

}




double _filtVALLEY(double *vals, long nvals ) { 

  double retval;
  short center;
  long ii;

  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }


  retval = centerval; /* assume valley */
  center = 0;
  for(ii=nvals;ii--;) {
    if (( vals[ii] == centerval ) && (center == 0 )) {
      center=1; /* want the min/max excludeding the center value */
      continue;
    }
    if ( vals[ii] <= centerval ) {
      ds_MAKE_DNAN(retval);
      break;
    }
  }

  return(retval);

}

double _filtPLAIN(double *vals, long nvals ) { 
  /* Same as valley but < instead of <= */

  double retval;
  short center;
  long ii;

  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }


  retval = centerval; /* assume valley */
  center = 0;
  for(ii=nvals;ii--;) {
    if (( vals[ii] == centerval ) && (center == 0 )) {
      center=1; /* want the min/max excludeding the center value */
      continue;
    }
    if ( vals[ii] < centerval ) {
      ds_MAKE_DNAN(retval);
      break;
    }
  }

  return(retval);

}



double _filtOLYMPIC(double *vals, long nvals ) {
  double retval;
  double minval;
  double maxval;
  long ii;
  if (  3 > nvals ) {  /* need to have at least 3 points */
    ds_MAKE_DNAN(retval);
    return(retval);
  }

  retval = vals[0];
  minval = retval;
  maxval = retval;
  for(ii=nvals;--ii;) {
    maxval = MAX(maxval, vals[ii] );
    minval = MIN(minval, vals[ii] );
    retval += vals[ii];
  }

  retval -= minval;
  retval -= maxval;

  retval /= ( nvals-2.0);


  return(retval); 
}


double _filtPMEAN(double *vals, long nvals ) {
  double retval;
  long ii;
  long n0, nn;

  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }
  retval=0;
  nn = 0;
  n0 = 0;
  for (ii=nvals;ii--;) {
    long ival = vals[ii];
    
    if ( (ival < 0) || (ival > 1 )) continue;
    nn++;
    if ( ival == 0 ) n0++;
  }
  
  if ( 0 == n0 ) 
    retval = _filtMEDIAN( vals, nvals );
  else 
    retval = ((1.0*nn)/n0) - 1.0;

  return(retval); 
}



double _filtMU3(double *vals, long nvals ) {
  double retval;
  long ii;
  double mean;
  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }
  mean=vals[0];
  for (ii=nvals;--ii;)
    mean += vals[ii];
  mean /= nvals; 
  retval = 0;
  for (ii=nvals;ii--;) {
    double diff = (vals[ii]-mean);
    retval += pow(diff, 3.0);
  }
  retval /= nvals;

  retval /= _filtSUM( vals, nvals );
  return(retval);
}


double _filtMU4(double *vals, long nvals ) {
  double retval;
  long ii;
  double mean;
  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }
  mean=vals[0];
  for (ii=nvals;--ii;)
    mean += vals[ii];
  mean /= nvals; 
  retval = 0;
  for (ii=nvals;ii--;) {
    double diff = (vals[ii]-mean);
    retval += pow(diff, 4.0);
  }
  retval /= nvals;
  retval /= _filtSUM( vals, nvals );
  return(retval);
}

double _filtJITTER(double *vals, long nvals ) {
  double retval;
  long idx;

  static int settime=0;

  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }

  if (0 == settime) {
    srand48( time( NULL ));
    settime=1;
  }

  retval = drand48() * nvals;

  idx = retval;
  if ( idx >= nvals ) idx = (nvals-1);

  retval=vals[idx];

  return(retval);

}


double _filtRMS(double *vals, long nvals ) { 
  double retval;
  long ii;

  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }

  retval = vals[0]*vals[0];
  for(ii=nvals;--ii;)
    retval += (vals[ii]*vals[ii]);

  retval /= nvals;
  retval = sqrt(retval);

  return(retval); 
}



double _filt3SIGMA_MEAN(double *vals, long nvals ) { 
  double retval;
  long ii;

  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }

  double *clip_vals = ( double*)calloc( nvals, sizeof(double));
  for (ii=nvals; ii--; ) {
    clip_vals[ii] = vals[ii];
  }
  
  long nclip = nvals;
  double mean;
  double sigm;

  long atiter = 0;

  while ( atiter < SIGMA_ITER ) {
    
    mean = _filtMEAN( clip_vals, nclip );
    sigm = _filtSIG( clip_vals ,nclip );
    sigm *= 3; /* 3 sigma */
    
    nclip = 0;
    for (ii=nvals; ii--; ) {
      if ( fabs( vals[ii] - mean ) < sigm ) {
	clip_vals[nclip] = vals[ii];
	nclip++;
      }
    }

    atiter++;
  } /* end while */

  retval = _filtMEAN( clip_vals, nclip );

  free( clip_vals );
  return(retval); 
}




double _filt3SIGMA_MEDIAN(double *vals, long nvals ) { 
  double retval;
  long ii;

  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }

  double *clip_vals = ( double*)calloc( nvals, sizeof(double));
  for (ii=nvals; ii--; ) {
    clip_vals[ii] = vals[ii];
  }
  
  long nclip = nvals;
  double mean;
  double sigm;

  long atiter = 0;

  while ( atiter < SIGMA_ITER ) {
    
    mean = _filtMEDIAN( clip_vals, nclip );
    sigm = _filtSIG( clip_vals ,nclip );
    sigm *= 3; /* 3 sigma */
    
    nclip = 0;
    for (ii=nvals; ii--; ) {
      if ( fabs( vals[ii] - mean ) < sigm ) {
	clip_vals[nclip] = vals[ii];
	nclip++;
      }
    }

    atiter++;
  } /* end while */

  retval = _filtMEDIAN( clip_vals, nclip );

  free( clip_vals );
  return(retval); 
}







static int _filtSLOPE_sort( const void* a1, const void *a2 )
{
  double *aa = (double*)a1;
  double *bb = (double*)a2;
  
  if ( *aa < *bb )
    return(-1);
  if ( *aa > *bb )
    return(1);

  return(0);
}


static double* __get_slope( double *vals, long nvals )
{
  qsort( vals, nvals, sizeof(double), _filtSLOPE_sort);
  double *slope = (double*)calloc( nvals, sizeof(double));
  long ii;
  
  slope[0] = (vals[1] - vals[0]);
  slope[nvals-1] = (vals[nvals-1] - vals[nvals-2]);
  
  for (ii=(nvals-2);ii>=1;ii--) {
    slope[ii] = (vals[ii+1]-vals[ii]);
  }
  
  return(slope);
}



double _filtMIN_SLOPE(double *vals, long nvals) {
  double retval;
  if ( 0 == nvals ) {
    ds_MAKE_DNAN(retval);
    return(retval);
  }
  
  if ( 1 == nvals ) {
    return(vals[0]);
  }

  double *slope = __get_slope( vals, nvals );
  double minslope = fabs(slope[nvals-1]);
  retval = vals[nvals-1];
  
  long ii;
  for (ii=nvals-1;ii--;) { 
    if ( fabs(slope[ii]) <= minslope ) {
      minslope = fabs(slope[ii]);
      retval = vals[ii];
    }
  }

  
  
  free(slope);

  return(retval);

}
