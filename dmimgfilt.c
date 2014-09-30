/*                                                                
**  Copyright (C) 2004-2008,2011  Smithsonian Astrophysical Observatory 
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


#include <dslib.h>
#include <dsnan.h>
#include <histlib.h>
#include <math.h>
#include <float.h>
#include <cxcregion.h>
#include <dmimgfilt.h>

typedef enum { filtUNKN, filtMIN, filtMAX, filtCOUNT,
	       filtMEAN, filtMEDIAN, filtMODE, filtSIG,
               filtEXTREME, filtLOCEQ, filtKUWAHARA,
               filtUNSHARP, filtRANGE, filtVARIANCE,
	       filtNORM_MODE, filtQUANTILE_25,
	       filtQUANTILE_33, filtQUANTILE_67,
	       filtQUANTILE_75, filtQUANTILE_XX,           
           filtMOST_COMMON,
               filtSUM, filtRCLIP, 
               filtPEAK, filtVALLEY,
               filtRIDGE, filtPLAIN,
               filtOLYMPIC, filtPMEAN, filtMID, filtMU3, filtMU4,
               filtJITTER, filtRMS, filtMIN_SLOPE,
	       filt3SIGMA_MEDIAN, filt3SIGMA_MEAN
} Method;

int dmimgfilter(void);

short *get_image_mask( dmBlock *inBlock, void *data, dmDataType dt, 
		       long *lAxes, regRegion *dss, long null, short has_null, 
		       dmDescriptor *xAxis, dmDescriptor *yAxis );


double get_image_value( void *data, dmDataType dt, 
			long xx, long yy, long *lAxes, 
			short *mask );

dmDataType get_image_data( dmBlock *inBlock, void **data, long **lAxes,
			   regRegion **dss, long *nullval, short *nullset );

short  get_image_wcs( dmBlock *imgBlock, dmDescriptor **xAxis, 
		      dmDescriptor **yAxis );

long evaluate_kernel( char *kernel, double **kx, double **ky );

Method get_method( char *func ) ;


double fill_data( Method mth, double *vals, long nvals ) ;




short *get_image_mask( dmBlock *inBlock, void *data, dmDataType dt, 
		       long *lAxes, regRegion *dss, long null, short has_null, 
		       dmDescriptor *xAxis, dmDescriptor *yAxis )
{
  long npix = lAxes[0] * lAxes[1];
  short *mask;
  long xx, yy;
  mask = (short*)calloc( npix, sizeof(short));

  
  for ( xx=lAxes[0]; xx--; ) {
    for ( yy=lAxes[1]; yy--; ) {
      double dat;
      long idx;
      idx = xx + ( yy * lAxes[0] );
      
      dat = get_image_value( data, dt, xx, yy, lAxes, NULL );
      
      /* Now ... if it is an integer data type, it could possibly have a
	 null value. Check for that */
      if ( ( has_null && ( dat == null ) ) ||
	   ds_dNAN( dat ) ) {
	continue;
      }
      
      
      
      /* If the image has a data sub space (aka a region filter applied)
	 then need to convert coords to physical and check */
      if ( dss && xAxis ) {
	double pos[2];
	double loc[2];
	pos[0]=xx+1;
	pos[1]=yy+1;
	
	if (yAxis) {  /* If no y axis, then xAxis has 2 components */
	  dmCoordCalc_d( xAxis, pos, loc );
	  dmCoordCalc_d( yAxis, pos+1, loc+1 );
	} else {
	  dmCoordCalc_d( xAxis, pos, loc );
	}
	if ( !regInsideRegion( dss, loc[0], loc[1] ) )
	  continue;
      }
      
      mask[idx] = 1;
    }
  }


  return(mask );
}



double get_image_value( void *data, dmDataType dt, 
			long xx, long yy, long *lAxes, 
			short *mask )
{

  long npix = xx + (yy * lAxes[0] );
  double retval;

  /* Okay, first get all the data from the different data types.  
     Cast everything to doubles */

  switch ( dt ) {
    
  case dmBYTE: {
    unsigned char *img = (unsigned char*)data;
    retval = img[npix];
    break;
  }
    
  case dmSHORT: {
    short *img = (short*)data;
    retval = img[npix];
    break;
  }
    
  case dmUSHORT: {
    unsigned short *img = (unsigned short*)data;
    retval = img[npix];
    break;
  }
    
  case dmLONG: {
    long *img = (long*)data;
    retval = img[npix];
    break;
  }
    
  case dmULONG: {
    unsigned long *img = (unsigned long*)data;
    retval = img[npix];
    break;
  }
    
  case dmFLOAT: {
    float *img = (float*)data;
    retval = img[npix];
    break;
  }
  case dmDOUBLE: {
    double *img = (double*)data;
    retval = img[npix];
    break;
  }
  default:
    ds_MAKE_DNAN( retval );

  }


  if ( mask ) {
    if ( !mask[npix] ) {
      ds_MAKE_DNAN( retval );
    }
  }


  return(retval);

}



/* Load the data into memory,  check for DSS, null values */
dmDataType get_image_data( dmBlock *inBlock, void **data, long **lAxes,
			   regRegion **dss, long *nullval, short *nullset )
{

  dmDescriptor *imgDesc;
  dmDataType dt;
  dmDescriptor *grp;
  dmDescriptor *imgdss;

  long naxes;
  long npix;
  char ems[1000];

  *nullval = INDEFL;
  *dss = NULL;
  *nullset = 0;
  
  imgDesc = dmImageGetDataDescriptor( inBlock );

  /* Sanity check, only 2D images */
  naxes = dmGetArrayDimensions( imgDesc, lAxes );
  if ( naxes != 2 ) {
    return( dmUNKNOWNTYPE );
  }
  npix = (*lAxes)[0] * (*lAxes)[1];
  dt = dmGetDataType( imgDesc );


  /* Okay, first lets get the image descriptor */
  grp = dmArrayGetAxisGroup( imgDesc, 1 );
  dmGetName( grp, ems, 1000);
  imgdss = dmSubspaceColOpen( inBlock, ems );
  if ( imgdss )
    *dss = dmSubspaceColGetRegion( imgdss);
  
  
  switch ( dt ) 
    {
    case dmBYTE:
      *data = ( void *)calloc( npix, sizeof(char ));
      dmGetArray_ub( imgDesc, (unsigned char*) *data, npix );
      if ( dmDescriptorGetNull_l( imgDesc, nullval) == 0 ) {
	*nullset=0;
      } else
	*nullset=1;
      break;
      
    case dmSHORT:
      *data = ( void *)calloc( npix, sizeof(short ));
      dmGetArray_s( imgDesc, (short*) *data, npix );
      if ( dmDescriptorGetNull_l( imgDesc, nullval) == 0 ) {
	*nullset=0;
      } else
	*nullset=1;
      break;
      
    case dmUSHORT:
      *data = ( void *)calloc( npix, sizeof(short ));
      dmGetArray_us( imgDesc, (unsigned short*) *data, npix );
      if ( dmDescriptorGetNull_l( imgDesc, nullval) == 0 ) {
	*nullset=0;
      } else
	*nullset=1;
      break;
      
    case dmLONG:
      *data = ( void *)calloc( npix, sizeof(long ));
      dmGetArray_l( imgDesc, (long*) *data, npix );
      if ( dmDescriptorGetNull_l( imgDesc, nullval) == 0 ) {
	*nullset=0;
      } else
	*nullset=1;
      break;
      
    case dmULONG:
      *data = ( void *)calloc( npix, sizeof(long ));
      dmGetArray_ul( imgDesc, (unsigned long*) *data, npix );
      if ( dmDescriptorGetNull_l( imgDesc, nullval) == 0 ) {
	*nullset=0;
      } else
	*nullset=1;
      break;
      
    case dmFLOAT:
      *data = ( void *)calloc( npix, sizeof(float ));
      dmGetArray_f( imgDesc, (float*) *data, npix );
      *nullset = 0;
      break;
      
    case dmDOUBLE:
      *data = ( void *)calloc( npix, sizeof(double ));
      dmGetArray_d( imgDesc, (double*) *data, npix );
      *nullset = 0;
      break;
      
    default:
      return( dmUNKNOWNTYPE );
    }

  return(dt);

}



long evaluate_kernel( char *kernel, double **kx, double **ky )
{
  regRegion *reg;
  
  double dx[2] = { 0, 0};
  double dy[2] = { 0, 0};
  long xx, yy;

  double fx[2] = { -FLT_MAX, FLT_MAX };
  double fy[2] = { -FLT_MAX, FLT_MAX };
  
  long retval=0;

  if ( NULL == (reg=dmRegParse( kernel ))) {
    err_msg("ERROR: problems parsing region string '%s'\n", kernel );
    return(-1);
  }

  regExtent( reg, fx, fy, dx, dy );

  xx = dx[1] - dx[0];
  yy = dy[1] - dy[0];


  dx[0] = ((short)dx[0]);
  dx[1] = ((short)dx[1]);
  dy[0] = ((short)dy[0]);
  dy[1] = ((short)dy[1]);

  if ( 0 != regRegionToList(reg, dx[0], dx[1], dy[0], dy[1],
			    1, kx, ky, &retval )) {
    return(-1);
  }

  return(retval);

}




/* Get the WCS descriptor */
short  get_image_wcs( dmBlock *imgBlock, dmDescriptor **xAxis, 
		      dmDescriptor **yAxis )
{
  

  dmDescriptor *imgData;
  long n_axis_groups;

  imgData = dmImageGetDataDescriptor( imgBlock );
  n_axis_groups = dmArrayGetNoAxisGroups( imgData );
  

  /* This is the usual trick ... can have 1 axis group w/ 
     dimensionality 2 (eg a vector column) or can have
     2 axis groups w/ dimensionaity 1 (eg 2 disjoint columns)*/

  if ( n_axis_groups == 1 ) {
    dmDescriptor *pos = dmArrayGetAxisGroup( imgData, 1 );
    dmDescriptor *xcol;
    long n_components;
    
    n_components = dmGetElementDim( pos );
    if ( n_components != 2 ) {
      err_msg("ERROR: could not find 2D image\n");
      return(-1);
    }
    
    xcol = dmGetCpt( pos, 1 );
    
    *xAxis = pos;
    *yAxis = NULL;
    
  } else if ( n_axis_groups == 2 ) {
    dmDescriptor *xcol;
    dmDescriptor *ycol;
  
    xcol = dmArrayGetAxisGroup( imgData, 1 );
    ycol = dmArrayGetAxisGroup( imgData, 2 );

    *xAxis = xcol;
    *yAxis = ycol;
    
  } else {
    err_msg("Invalid number of axis groups\n");
    *xAxis = NULL;
    *yAxis = NULL;
    return(-1);
  }

  return(0);

}

Method get_method( char *func ) 
{
  Method retval;
  
  if ( !strcmp( func, "min" ) ) {
    retval = filtMIN; 
  } else if ( !strcmp( func, "max" ) ) {
    retval = filtMAX; 
  } else if ( !strcmp( func, "mean" ) ) {
    retval = filtMEAN; 
  } else if ( !strcmp( func, "sum" ) ) {
    retval = filtSUM;
  } else if ( !strcmp( func, "median" ) ) {
    retval = filtMEDIAN; 
  } else if ( !strcmp( func, "mode" ) ) {
    retval = filtMODE; 
  } else if ( !strcmp( func, "count" ) ) {
    retval = filtCOUNT;
  } else if ( !strcmp( func, "sigma" ) ) {
    retval = filtSIG; 
  } else if ( !strcmp( func, "extreme" ) ) {
    retval = filtEXTREME; 
  } else if ( !strcmp( func, "locheq" ) ) {
    retval = filtLOCEQ; 
  } else if ( !strcmp( func, "kuwahara" ) ) {
    retval = filtKUWAHARA; 
  } else if ( !strcmp( func, "unsharp" ) ) {
    retval = filtUNSHARP; 
  } else if ( !strcmp( func, "range" ) ) {
    retval = filtRANGE; 
  } else if ( !strcmp( func, "variance" ) ) {
    retval = filtVARIANCE; 
  } else if ( !strcmp( func, "nmode" ) ) {
    retval = filtNORM_MODE; 
  } else if ( !strcmp( func, "q25" ) ) {
    retval = filtQUANTILE_25; 
  } else if ( !strcmp( func, "q33" ) ) {
    retval = filtQUANTILE_33; 
  } else if ( !strcmp( func, "q67" ) ) {
    retval = filtQUANTILE_67; 
  } else if ( !strcmp( func, "q75" ) ) {
    retval = filtQUANTILE_75; 
  } else if ( !strcmp( func, "qxx" ) ) {
    retval = filtQUANTILE_XX; 
  } else if ( !strcmp( func, "mcv" ) ) {
    retval = filtMOST_COMMON; 
  } else if ( !strcmp( func, "rclip" ) ) {
    retval = filtRCLIP; 
  } else if ( !strcmp( func, "peak" ) ) {
    retval = filtPEAK; 
  } else if ( !strcmp( func, "ridge" ) ) {
    retval = filtRIDGE; 
  } else if ( !strcmp( func, "valley" ) ) {
    retval = filtVALLEY; 
  } else if ( !strcmp( func, "plain" ) ) {
    retval = filtPLAIN; 
  } else if ( !strcmp( func, "olympic" ) ) {
    retval = filtOLYMPIC; 
  } else if ( !strcmp( func, "pmean" ) ) {
    retval = filtPMEAN; 
  } else if ( !strcmp( func, "mid" ) ) {
    retval = filtMID; 
  } else if ( !strcmp( func, "mu3" ) ) {
    retval = filtMU3; 
  } else if ( !strcmp( func, "mu4" ) ) {
    retval = filtMU4; 
  } else if ( !strcmp( func, "jitter" ) ) {
    retval = filtJITTER; 
  } else if ( !strcmp( func, "rms" ) ) {
    retval = filtRMS; 
  } else if ( !strcmp( func, "nslope" ) ) {
    retval = filtMIN_SLOPE; 
  } else if ( !strcmp( func, "3sigmean" ) ) {
    retval = filt3SIGMA_MEAN; 
  } else if ( !strcmp( func, "3sigmedian" ) ) {
    retval = filt3SIGMA_MEDIAN; 
  } else {
    retval = filtUNKN;
  }
  
  return(retval);

}



int dmimgfilter(void)
{

  char infile[DS_SZ_PATHNAME];
  char outfile[DS_SZ_PATHNAME];
  char func[10];
  char mask[DS_SZ_FNAME];
  short clobber;
  short verbose;

  void **data = NULL;
  long *lAxes = NULL;
  regRegion *dss = NULL;
  long null;
  short has_null;

  Stack instack;
  char *stkfile;
  long nfiles;
  long ii;
  long kk;
  Header_Type **hdr;
  char lookup[DS_SZ_PATHNAME];

  dmDataType *dt;
  dmBlock *inBlock;
  dmDescriptor *xdesc, *ydesc;
  
  dmBlock *outBlock = NULL;
  dmDescriptor *outDesc = NULL;

  long nkpix;
  double *kx, *ky;


  double *vals;
  long nvals;

  long xx, yy;

  double *outdata = NULL;
  Method mth;
  long nkvals;
  short **nullmask;

  unsigned long outpix;

  double (*nonlinear)( double *vals, long nvals ) = NULL;
  long niter;


  /* Get the parameters */
  clgetstr( "infile", infile, DS_SZ_PATHNAME );
  clgetstr( "outfile", outfile, DS_SZ_PATHNAME );
  clgetstr( "function", func, 10 );
  clgetstr( "mask", mask, DS_SZ_FNAME );
  niter = clgeti( "numiter" );
  clgetstr( "lookupTab", lookup, DS_SZ_PATHNAME );
  clobber = clgetb( "clobber" );
  verbose = clgeti( "verbose" );


  instack = stk_build( infile );
  if (( NULL == instack ) ||
      ( stk_count(instack) == 0 ) ||
      (( stk_count(instack) == 1 ) && ( strlen(stk_read_num(instack,1)) == 0 ) ) ) {
    err_msg("ERROR: Problems expanding stack of infiles '%s'\n", infile );
  }

  nfiles = stk_count(instack);
  dt = (dmDataType*)calloc(nfiles,sizeof(dmDataType));
  nullmask = ( short**)calloc(nfiles,sizeof(short*));
  data = ( void**)calloc(nfiles,sizeof(void*));
  hdr = (Header_Type**)calloc( nfiles, sizeof(Header_Type*));
  ii = 0;


  if ( ds_clobber( outfile, clobber, NULL ) != 0 ) {
    return(-1);
  }


  stk_rewind(instack);
  while ( NULL != (stkfile=stk_read_next(instack)) ) {
    long *locaxes;
    
    if ( NULL == ( inBlock = dmImageOpen( stkfile) ) ) {
      err_msg("ERROR: Cannot open image '%s'\n", infile );
      return(-1);
    }
    
    if ( dmUNKNOWNTYPE == ( dt[ii] = get_image_data( inBlock, &(data[ii]), &locaxes, 
						     &dss, &null, &has_null ))) {
      err_msg("ERROR: Cannot get image data or unknown image data-type for "
	      "file '%s'\n", infile);
      return(-1);
    }
    
    if ( ii == 0 ) {
      lAxes = locaxes;


      if ( NULL == ( outBlock = dmImageCreate( outfile, dmFLOAT,lAxes,2 ))){
	err_msg("ERROR: Cannot create output image '%s'\n", outfile );
	return(-1);
      }
      outDesc = dmImageGetDataDescriptor( outBlock );
      dmBlockCopy( inBlock, outBlock, "HEADER");
      ds_copy_full_header( inBlock, outBlock, "dmimgfilt", 0 );
      dmBlockCopyWCS( inBlock, outBlock);
      
    } else {
      /* check */
      if ( ( locaxes[0] != lAxes[0] ) ||
	   ( locaxes[1] != lAxes[1] )    ) {
	err_msg("ERROR: All images must be the same size\n");
	return(-1);
      }
    }
    
    if ( 0 != get_image_wcs( inBlock, &xdesc, &ydesc ) ) {
      err_msg("ERROR: Cannot load WCS for file '%s'\n", stkfile );
      return(-1);
    }
    
    if ( NULL == ( nullmask[ii] = get_image_mask( inBlock, data[ii], dt[ii], 
						  lAxes, dss, null, has_null, 
						  xdesc, ydesc ))){
      
    }
    hdr[ii] = getHdr( inBlock, hdrDM_FILE);
    ii++;
    dmImageClose( inBlock );
  } /* end while loop over stack */


  if ( -1 == ( nkpix = evaluate_kernel( mask, &kx, &ky ) ) ) {
    err_msg("ERROR: cannot parse mask function '%s'\n", mask );
    return(-1);
  }

  if ( 0 == nkpix ) {
    err_msg("ERROR: mask has no pixels in it\n");
    return(-1);
  }


  if ( filtUNKN == ( mth = get_method( func ) ) ) {
    err_msg("ERROR: function '%s' is not supported\n", func );
    return(-1);
  }


  /* Use function points so that we don't have to go through the switch
     statement for every pixel, also makes adding new f() easier.*/
  switch( mth ) {
  case filtCOUNT:       nonlinear = _filtCOUNT;       break;
  case filtMIN:         nonlinear = _filtMIN;         break;
  case filtMAX:         nonlinear = _filtMAX;         break;
  case filtSUM:         nonlinear = _filtSUM;         break;
  case filtEXTREME:     nonlinear = _filtEXTREME;     break;
  case filtRANGE:       nonlinear = _filtRANGE;       break;
  case filtRCLIP:       nonlinear = _filtRCLIP;       break;
  case filtLOCEQ:       nonlinear = _filtLOCEQ;       break;
  case filtMEAN:        nonlinear = _filtMEAN;        break;
  case filtMEDIAN:      nonlinear = _filtMEDIAN;      break;
  case filtMODE:        nonlinear = _filtMODE;        break;
  case filtNORM_MODE:   nonlinear = _filtNORM_MODE;   break;
  case filtSIG:         nonlinear = _filtSIG;         break;
  case filtKUWAHARA:    nonlinear = _filtKUWAHARA;    break;
  case filtUNSHARP:     nonlinear = _filtUNSHARP;     break;
  case filtVARIANCE:    nonlinear = _filtVARIANCE;    break;
  case filtQUANTILE_25: nonlinear = _filtQUANTILE_25; break;
  case filtQUANTILE_33: nonlinear = _filtQUANTILE_33; break;
  case filtQUANTILE_67: nonlinear = _filtQUANTILE_67; break;
  case filtQUANTILE_75: nonlinear = _filtQUANTILE_75; break;
  case filtQUANTILE_XX: nonlinear = _filtQUANTILE_XX; break;
  case filtMOST_COMMON: nonlinear = _filtMOST_COMMON; break;
  case filtPEAK:        nonlinear = _filtPEAK;        break;
  case filtRIDGE:       nonlinear = _filtRIDGE;       break;
  case filtVALLEY:      nonlinear = _filtVALLEY;      break;
  case filtPLAIN:       nonlinear = _filtPLAIN;       break;
  case filtOLYMPIC:     nonlinear = _filtOLYMPIC;     break;
  case filtPMEAN:       nonlinear = _filtPMEAN;       break;
  case filtMID:         nonlinear = _filtMID;         break;
  case filtMU3:         nonlinear = _filtMU3;         break;
  case filtMU4:         nonlinear = _filtMU4;         break;
  case filtJITTER:      nonlinear = _filtJITTER;      break;
  case filtRMS:         nonlinear = _filtRMS;         break;
  case filtMIN_SLOPE:   nonlinear = _filtMIN_SLOPE;   break;
  case filt3SIGMA_MEAN: nonlinear = _filt3SIGMA_MEAN; break;
  case filt3SIGMA_MEDIAN: nonlinear = _filt3SIGMA_MEDIAN; break;
  case filtUNKN:  /* fall thru intended */
  default:
    nonlinear = NULL;
    break;
  }


  if ( filtPMEAN == mth  ) {
    for (ii=0;ii<nfiles;ii++) {
      if (( dt[ii] == dmFLOAT ) || (dt[ii] == dmDOUBLE)) {
	err_msg("WARNING: File #%ld is real valued but 'pmean'"
		" method only useful for integer values, data"
		" will be cast to integer\n", ii );
      }
    }
  }

  minpixel = 0;
  maxpixel = 255;
  pixelrange = 256;


  nkvals = nkpix;


  if ( mth == filtKUWAHARA ) {
    for (xx=4;xx--;) {
      kuw_vals[xx] = (double*)calloc(nkvals*nfiles,sizeof(double));
    }
  }


  if ( mth == filtMOST_COMMON ) {
    num_hist_bins =sqrt( nkvals )+0.5;
    histogram = (long*)calloc(num_hist_bins, sizeof(long));
  }



  vals = (double*)calloc(nkvals*nfiles, sizeof(double));


  while ( niter ) {
  /* Instead of computing  xx + yy*lAxes[0], just keep a running
     counter of the current pixel number */


    if ( verbose>1) {
      fprintf(stderr, "Processing iteration %ld\n", niter );
    }

  outdata = (double*)calloc( lAxes[0]*lAxes[1], sizeof(double));
  outpix=( lAxes[0] * lAxes[1]) -1;
  

  /* Okay, first loop over yy and xx which are the image axes */
  for (yy=lAxes[1];yy--; ) {
    
    if ( verbose>2) {
      fprintf( stderr, "Processing %4.1f%% complete\r",
	       ((lAxes[1]-yy)/(0.01*lAxes[1])) );
    }


    for (xx=lAxes[0];xx--;) {
      long pix;      

      nvals = 0;
      pix = 0;
      kuw_num[0]=0; kuw_num[1] =0;kuw_num[2]=0;kuw_num[3]=0;

      /* Now loop over the number of files */
      for (ii=nfiles;ii--; ) {
	centerval = get_image_value( data[ii], dt[ii], xx, yy, 
				     lAxes, nullmask[ii] );
	if ( ds_dNAN( centerval ) ) {
	  continue;
	}

	pix=nkvals; /* Same game as w/ the output pixel */
	pix--;

	/* These are looping over the pixels in the kernel */

	for (kk=nkvals;kk--;) {
	  double dater;
	  long ay;
	  long ax;

	  ax=xx+kx[kk];
	  ay=yy+ky[kk];
	  
	  if ( (ax<0) || (ax >= lAxes[0])) {
	    pix--;
	    continue;
	  }
	  if ( (ay<0) || (ay >= lAxes[1])) {
	    pix--;
	    continue;
	  }
	  /* If the kernel pixel is non-zero then add more data to
	     the array */
	  dater = get_image_value( data[ii], dt[ii], ax, ay, lAxes, 
				   nullmask[ii] );
	  if ( ds_dNAN(dater) ) {
	    pix--;
	    continue;
	  }
	  
	  if ( mth == filtKUWAHARA ) {
	    long qx, qy;
	    qx = kx[kk];
	    qy = ky[kk];
	    
	    /* Figure out which quadrant we're in */
	    if ( ( qx >= 0 ) && ( qy >= 0 ) )
	      kuw_vals[0][ kuw_num[0]++] = dater;
	    if ( ( qx <= 0 ) && ( qy >= 0 ) )
	      kuw_vals[1][ kuw_num[1]++] = dater;
	    if ( ( qx <= 0 ) && ( qy <= 0 ) )
	      kuw_vals[2][ kuw_num[2]++] = dater;
	    if ( ( qx >= 0 ) && ( qy <= 0 ) )
	      kuw_vals[3][ kuw_num[3]++] = dater;
	    
	  } 
	  vals[nvals] = dater;
	  nvals++;
	  pix--;
	  
	} /* end for kk */

	
      } /* end loop over files, ii */

      
      if ( nvals == 0 ) {
	double nanval;
	ds_MAKE_DNAN(nanval);
	outdata[outpix] = nanval;
      } else {
	outdata[outpix] = nonlinear( vals, nvals ); /* 'nonlinear' is function pointer */
      }
      outpix--;
      
    } /* end for xx */

  } /* end for yy */

  if ( verbose>2) {
    fprintf(stderr,"\n");
  }
  /* Okay, if more than one iteration then we ignore the other files
     and just re-run outdata thru the algorithm again.

     One may think that instead the number of files would increase by
     one w/ the new data being a virtual dataset that should be added.
     This is one interp. that could also be done.
  */
  nfiles = 1;
  dt[0] = dmDOUBLE;
  free( data[0] );
  data[0] = outdata;
  nullmask[0] = NULL;
  
  niter --;
  } /* end while niter */


  nfiles = stk_count(instack); /* Re-count since while (niter) loop */
  if (( strlen(lookup) > 0 ) && 
      ( ds_strcmp_cis(lookup, "none" ) != 0 ) &&
      ( nfiles > 1 ) ) {
    Header_Type *mergehdr;
    mergehdr = mergeHdr( lookup, hdr, nfiles );
    putHdr( outBlock, hdrDM_FILE, mergehdr, ALL_STS, "dmimgfilt");
  }
  put_param_hist_info( outBlock, "dmimgfilt", NULL, 0 );


  dmSetArray_d( outDesc, outdata, (lAxes[0]*lAxes[1]));
  dmImageClose(outBlock );
  return(0);

}








