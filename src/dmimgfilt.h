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

#include "dmfilters.h"
#include "dmimgio.h"
#include "hdrlib2.h"

/* Hold info for an input image */
typedef struct {
    void *data;        // pixel values
    dmDataType dt;     // pixel datatype
    long *lAxes;       // axis lenghts
    short *mask;        // mask of valid pixels
    dmDescriptor *xdesc;  // X (or sky) coordinate descriptor
    dmDescriptor *ydesc;  // Y coordinate descriptor
    dmBlock *block; // The block image came from
    Header_Type *hdr;   // Header keywords
} Image;

extern Image* load_infile(char *infile);
extern int dmimgfilter(void);
extern long evaluate_kernel( char *kernel, double **kx, double **ky );
