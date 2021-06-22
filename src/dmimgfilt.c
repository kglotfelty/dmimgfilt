/*                                                                
**  Copyright (C) 2004-2008,2011,2021  Smithsonian Astrophysical Observatory 
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
#include <time.h>




/*
 *  Load images using dmimgio routines
 */
Image* load_infile(char *infile)
{
    // Load image

    Image *image;
    if (NULL == (image = calloc(1,sizeof(Image)))) {
        err_msg("ERROR: Cannot allocate memory for image\n");
        return(NULL);
    }

    if (NULL == (image->block = dmImageOpen(infile))) {
        err_msg("ERROR: Cannot load infile '%s'\n",infile);
        return(NULL);
    }

    // dmimgio
    regRegion *dss = NULL;
    long null_value;
    short has_null;
    image->dt = get_image_data(image->block, &(image->data),
                    &(image->lAxes), &dss, &null_value, &has_null);
    get_image_wcs(image->block, &(image->xdesc), &(image->ydesc));
    image->mask = get_image_mask(image->block, image->data,
                    image->dt, image->lAxes, dss, null_value,
                    has_null, image->xdesc, image->ydesc);

    if (dss != NULL){
        regFree(dss);
        dss=NULL;
    }

    image->hdr = getHdr(image->block, hdrDM_FILE);

    return(image);
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
  char errStr[ dmMAXSTR + 1]= "" ;

  if ( NULL == (reg=dmRegParse( kernel ))) {
    dmGetErrorMessage( errStr, dmMAXSTR );
    err_msg("ERROR: %s.\n", errStr);
    return(-1);
  }

  else{
    /*Parsed kernel will be in Physical coordinates but image is in Logical 
    * coordinates and not currently compatible. Remove this when a region 
    * library method is created to set coordinate system. 
    */
    if ( NULL != ( strstr( kernel, "mask(" ))){
      err_msg("ERROR: mask syntax, '%s', not compatible with dmimgfilt.  ", kernel);
      return(-1);
    }
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




int dmimgfilter(void)
{



  long ii;
  long kk;
  
  dmBlock *outBlock = NULL;
  dmDescriptor *outDesc = NULL;

  double *kx, *ky;


  double *vals;
  long nvals;

  long xx, yy;

  double *outdata = NULL;
  long nkvals;
  short **nullmask;

  unsigned long outpix;

  double (*nonlinear)( double *vals, long nvals ) = NULL;

  char unit[DS_SZ_KEYWORD];
  memset( &unit[0], 0, DS_SZ_KEYWORD) ;


  /* Get the parameters */
  char infile[DS_SZ_PATHNAME];
  char outfile[DS_SZ_PATHNAME];
  char func[10];
  char mask[DS_SZ_FNAME];
  long niter;
  char lookup[DS_SZ_PATHNAME];
  short clobber;
  short verbose;
  long seed;

  clgetstr( "infile", infile, DS_SZ_PATHNAME );
  clgetstr( "outfile", outfile, DS_SZ_PATHNAME );
  clgetstr( "function", func, 10 );
  clgetstr( "mask", mask, DS_SZ_FNAME );
  niter = clgeti( "numiter" );
  clgetstr( "lookupTab", lookup, DS_SZ_PATHNAME );
  seed = clgeti("randseed");  
  verbose = clgeti( "verbose" );
  clobber = clgetb( "clobber" );

  if ( 0 == seed) {
      srand48(time(NULL));
  } else {
      srand48(seed);
  }

  /* Setup input stack */
  Stack instack;
  char *stkfile;
  long nfiles;

  instack = stk_build( infile );
  if (( NULL == instack ) ||
      ( stk_count(instack) == 0 ) ||
      (( stk_count(instack) == 1 ) && ( strlen(stk_read_num(instack,1)) == 0 ) ) ) {
    err_msg("ERROR: Problems expanding stack of infiles '%s'\n", infile );
  }

  nfiles = stk_count(instack);

  /* Load input images */
  Image **images = (Image**)calloc(nfiles,sizeof(Image*));
  if ( NULL == images) {
      err_msg("ERROR: problem allocing memory");
      return(-1);
  }


  stk_rewind(instack);
  ii = 0;
  while ( NULL != (stkfile=stk_read_next(instack)) ) {

    Image *in_image;
    if ( NULL == (in_image = load_infile( stkfile ))) {
        return(-1);
    }
    images[ii]=in_image;
    
    if ( ii == 0 ) { // First image, create output image

      if ( ds_clobber( outfile, clobber, NULL ) != 0 ) {
        return(-1);
      }

      if ( NULL == ( outBlock = dmImageCreate( outfile, dmFLOAT,in_image->lAxes,2 ))){
        err_msg("ERROR: Cannot create output image '%s'\n", outfile );
        return(-1);
      }
      dmDescriptor *inDesc = NULL;
      outDesc = dmImageGetDataDescriptor( outBlock );
      inDesc = dmImageGetDataDescriptor( in_image->block);
      dmGetUnit( inDesc, unit, DS_SZ_KEYWORD );
      dmSetUnit( outDesc, unit );

      dmBlockCopy( in_image->block, outBlock, "HEADER");
      ds_copy_full_header( in_image->block, outBlock, "dmimgfilt", 0 );
      dmBlockCopyWCS( in_image->block, outBlock);
      
    } else {
      /* check , make sure next image is compatible with 1st*/
      if ( ( images[0]->lAxes[0] != in_image->lAxes[0] ) ||
           ( images[0]->lAxes[1] != in_image->lAxes[1] )    ) {
        err_msg("ERROR: All images must be the same size\n");
        return(-1);
      }
    }
    
    ii++;
    dmImageClose( in_image->block );
  } /* end while loop over stack */


  /* Merge headers if lookupTab is not-none */
  stk_rewind(instack);
  nfiles = stk_count(instack); /* Re-count since while (niter) loop */
  if (( strlen(lookup) > 0 ) && 
      ( ds_strcmp_cis(lookup, "none" ) != 0 ) &&
      ( nfiles > 1 ) ) {
    Header_Type **all_headers;
    if (NULL == (all_headers=calloc(nfiles,sizeof(Header_Type*)))) {
        err_msg("ERROR allocing memory");
        return(1);
    }
    for (ii=0;ii<nfiles;ii++) {
        all_headers[ii] = images[ii]->hdr;
    }

    Header_Type *mergehdr;
    mergehdr = mergeHdr( lookup, all_headers, nfiles );
    putHdr( outBlock, hdrDM_FILE, mergehdr, ALL_STS, "dmimgfilt");
    free(all_headers);
  }

  /* Add history to output */
  put_param_hist_info( outBlock, "dmimgfilt", NULL, 0 );


  /* Setup the mask and the filter function */


  if ( -1 == ( nkvals = evaluate_kernel( mask, &kx, &ky ) ) ) {
    err_msg("ERROR: cannot parse mask function '%s'\n", mask );
    return(-1);
  }

  if ( 0 == nkvals ) {
    err_msg("ERROR: mask has no pixels in it\n");
    return(-1);
  }

  if ( NULL == ( nonlinear = get_method( func ) ) ) {
    err_msg("ERROR: function '%s' is not supported\n", func );
    return(-1);
  }


  /* Special cases for various filters */

  minpixel = 0;       // From dmfilters 
  maxpixel = 255;     // From dmfilters
  pixelrange = 256;   // From dmfilters

  if ( _filtPMEAN == nonlinear  ) {
    for (ii=0;ii<nfiles;ii++) {
      if (( images[ii]->dt == dmFLOAT ) || (images[ii]->dt == dmDOUBLE)) {
        err_msg("WARNING: File #%ld is real valued but 'pmean'"
                " method only useful for integer values, data"
                " will be cast to integer\n", ii );
      }
    }
  }

  if ( _filtKUWAHARA == nonlinear ) {
    for (xx=4;xx--;) {
      kuw_vals[xx] = (double*)calloc(nkvals*nfiles,sizeof(double));
    }
  }

  if ( _filtMOST_COMMON == nonlinear ) {
    num_hist_bins =sqrt( nkvals )+0.5;
    histogram = (long*)calloc(num_hist_bins, sizeof(long));
  }


  // Setup buffer for stats value 
  vals = (double*)calloc(nkvals*nfiles, sizeof(double));


  while ( niter ) {
  /* Instead of computing  xx + yy*lAxes[0], just keep a running
     counter of the current pixel number */


    if ( verbose>1) {
      fprintf(stderr, "Processing iteration %ld\n", niter );
    }

  outdata = (double*)calloc( images[0]->lAxes[0]*images[0]->lAxes[1], sizeof(double));
  outpix=( images[0]->lAxes[0] * images[0]->lAxes[1]) -1;
  

  /* Okay, first loop over yy and xx which are the image axes */
  for (yy=images[0]->lAxes[1];yy--; ) {
    
    if ( verbose>2) {
      fprintf( stderr, "Processing %4.1f%% complete\r",
               ((images[0]->lAxes[1]-yy)/(0.01*images[0]->lAxes[1])) );
    }


    for (xx=images[0]->lAxes[0];xx--;) {
      nvals = 0;
      kuw_num[0]=0; kuw_num[1] =0;kuw_num[2]=0;kuw_num[3]=0;

      /* Now loop over the number of files */
      for (ii=nfiles;ii--; ) {
        centerval = get_image_value( images[ii]->data, images[ii]->dt, 
                                     xx, yy, 
                                     images[ii]->lAxes, images[ii]->mask );

        if ( ds_dNAN( centerval ) ) {
          continue;
        }

        //pix=nkvals; /* Same game as w/ the output pixel */
        //pix--;

        /* These are looping over the pixels in the kernel */

        for (kk=nkvals;kk--;) {
          double dater;
          long ay;
          long ax;

          ax=xx+kx[kk];  // kx,ky are offsets from center of mask
          ay=yy+ky[kk];
          
          if ( (ax<0) || (ax >= images[ii]->lAxes[0])) {
            continue;
          }

          if ( (ay<0) || (ay >= images[ii]->lAxes[1])) {
            continue;
          }
          /* If the kernel pixel is finite then add more data to
             the array */
          dater = get_image_value( images[ii]->data, images[ii]->dt, 
                                     ax, ay, 
                                     images[ii]->lAxes, images[ii]->mask );
          if ( ds_dNAN(dater) ) {
            continue;
          }

          
          if (_filtKUWAHARA == nonlinear) {
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
  images[0]->dt = dmDOUBLE;
  free( images[0]->data);
  images[0]->data = outdata;
  images[0]->mask = NULL;
  
  niter --;
  } /* end while niter */



  dmSetArray_d( outDesc, outdata, (images[0]->lAxes[0]*images[0]->lAxes[1]));
  dmImageClose(outBlock );
  return(0);

}








