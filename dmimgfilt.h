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

extern double _filtCOUNT(double *vals, long nvals ); 
extern double _filtMIN(double *vals, long nvals ); 
extern double _filtMAX(double *vals, long nvals );
extern double _filtSUM(double *vals, long nvals );
extern double _filtEXTREME(double *vals, long nvals );
extern double _filtRANGE(double *vals, long nvals );
extern double _filtRCLIP(double *vals, long nvals );
extern double _filtLOCEQ(double *vals, long nvals );
extern double _filtMEAN(double *vals, long nvals );
extern double _filtMEDIAN(double *vals, long nvals );
extern double _filtMODE(double *vals, long nvals );
extern double _filtNORM_MODE(double *vals, long nvals );
extern double _filtSIG(double *vals, long nvals );
extern double _filtKUWAHARA(double *vals, long nvals );
extern double _filtUNSHARP(double *vals, long nvals );
extern double _filtVARIANCE(double *vals, long nvals );
extern double _filtQUANTILE_25(double *vals, long nvals );
extern double _filtQUANTILE_33(double *vals, long nvals );
extern double _filtQUANTILE_67(double *vals, long nvals );
extern double _filtQUANTILE_75(double *vals, long nvals );
extern double _filtQUANTILE_XX(double *vals, long nvals );
extern double _filtMOST_COMMON(double *vals, long nvals );
extern double _filtPEAK(double *vals, long nvals );
extern double _filtVALLEY(double *vals, long nvals );
extern double _filtRIDGE(double *vals, long nvals );
extern double _filtPLAIN(double *vals, long nvals );
extern double _filtOLYMPIC(double *vals, long nvals );
extern double _filtPMEAN(double *vals, long nvals );
extern double _filtMID(double *vals, long nvals );
extern double _filtMU3(double *vals, long nvals );
extern double _filtMU4(double *vals, long nvals );
extern double _filtJITTER(double *vals, long nvals );
extern double _filtRMS(double *vals, long nvals );
extern double _filtMIN_SLOPE(double *vals, long nvals );
extern double _filt3SIGMA_MEDIAN(double *vals, long nvals );
extern double _filt3SIGMA_MEAN(double *vals, long nvals );


extern double minpixel;
extern double maxpixel;
extern double pixelrange;

extern double *kuw_vals[4];
extern long kuw_num[4];

extern long *histogram;
extern long num_hist_bins;

extern double centerval;

