#! /bin/sh

# 30 January 2002

# This is the official template for pipetool regression test scripts.
# In addition to supporting the "SHORTTEST" option, this script also
# allows the user to run individual subtests from the command line.
# The script will accept a series of test identifiers, generally of
# the form "test1" "test2" ... which are to be run.

# Portions of the script which must be customized are marked with "!!",
# below.

# The complete list of tests must be placed in "alltests", !!3, below.
# The test[s] for the SHORTTEST must be placed in "shortlist", !!4 below.


# !!1
# dmimgfilt.t
# test script for dmimgfilt


# !!2
# syntax:
# dmimgfilt.t [<testid> ... ]
 



######################################################################
# subroutine
# error_exit <message>
# Fatal error exit

error_exit()
{
  echo "$1" | tee -a $LOGFILE
  echo "${toolname} : FAIL" | tee -a $LOGFILE
  exit 1
}

######################################################################
# Initialization

# !!3
toolname="dmimgfilt"

# set up list of tests
# !!4
alltests="test_min test_max test_mean test_median test_mode test_sigma test_extreme test_locheq test_kuwahara test_unsharp  min_calcgti test_mask_empty test_mask_path"

# "short" test to run
# !!5
shortlist="$alltests"


# compute date string for log file
DT=`date +'%d%b%Y_%T'`

# check for directory environment variables
if test "x${TESTIN}" = "x" -o "x${TESTOUT}" = "x" -o "x${TESTSAV}" = "x" \
   -o "x${TESTLOG}" = "x" ; then
  error_exit "one or more of TESTIN/TESTOUT/TESTSAV/TESTLOG not defined" 
fi


# convenience definitions
OUTDIR=$TESTOUT/$toolname
SAVDIR=$TESTSAV/$toolname
INDIR=$TESTIN/$toolname
LOGDIR=$TESTLOG/$toolname

# set up log file name
LOGFILE=$LOGDIR/${toolname}_log.$DT

#get rid of old logs
rm -f $LOGDIR/${toolname}_log.*

# Any tests specified on command line?
if test $# -gt 0; then
  # yes, do those tests
  testlist=$*
else
  # No, see if we are to do "short" test
  if test "x$SHORTTEST" = "x" ; then
    # No, do everything
    testlist=$alltests
  else
    # yes, do short test
    testlist=$shortlist
  fi
fi


# Make sure we have a log directory
if test -d $LOGDIR ; then
 :
else
  mkdir -p $LOGDIR 
  if test $? -ne 0 ; then
    error_exit ""
  fi
fi


# Make sure we have an output directory
if test -d $OUTDIR ; then
 :
else
  mkdir -p $OUTDIR >> $LOGFILE 2>&1
  if test $? -ne 0 ; then
    error_exit "can't create output directory $OUTDIR"
  fi
fi

# announce ourselves
echo ""
echo "${toolname} regression" | tee $LOGFILE
echo ""

# All parameters except verbose should be set anyway, but clear them
# to be safe.
bose=`pget $toolname verbose`
punlearn $toolname
pset $toolname verbose=$bose

lkTab1=${ASCDS_CALIB}/dmmerge_header_lookup.txt
lkTab2=${INDIR}/calcgti_dmmerge_header_lookup.txt

script_succeeded=0

######################################################################
# Begin per-test loop

for testid in $testlist
do
    
  # delete old outputs
  rm -f $OUTDIR/${testid}*

  # Set up file names
  if [ ${testid} = "test_mask_empty" ] || [ ${testid} = "test_mask_path" ] ; then
    outfile=$OUTDIR/${testid}.fits
  else
    outfile=$OUTDIR/${testid}.fits"[IMAGE]"
  fi
  savfile=$SAVDIR/${testid}.fits

  echo "running $testid" >> $LOGFILE

  ####################################################################
  # run the tool
  case ${testid} in

    # !!6
    test_min ) test1_string="dmimgfilt infile=$INDIR/img.fits outfile="${outfile}" func=min mask='box(0,0,3,3)' clob+ lookupTab=${lkTab1}"
            ;;

    #
    ## same as test_min except using calcGTI/lkTab2
    #
    min_calcgti ) test1_string="dmimgfilt infile=$INDIR/img.fits outfile="${outfile}" func=min mask='box(0,0,3,3)' clob+ lookupTab=${lkTab2}"
            ;;

    test_max ) test1_string="dmimgfilt infile=$INDIR/img.fits outfile="${outfile}" func=max mask='circle(0,0,5)' clob+ lookupTab=${lkTab1}"
            ;;
    test_mean ) test1_string="dmimgfilt infile=$INDIR/img.fits outfile="${outfile}" func=mean mask='box(0,0,7,7)-circle(0,0,3)' clob+ lookupTab=${lkTab1}"
            ;;
    test_median ) test1_string="dmimgfilt infile=$INDIR/img.fits outfile="${outfile}" func=median mask='annulus(0,0,3,7)' clob+ lookupTab=${lkTab1}"
            ;;
    test_mode ) test1_string="dmimgfilt infile=$INDIR/img.fits outfile="${outfile}" func=mode mask='box(0,0,5,5)' clob+ lookupTab=${lkTab1}"
            ;;
    test_sigma ) test1_string="dmimgfilt infile=$INDIR/img.fits outfile="${outfile}" func=sigma mask='circle(0,0,7)' clob+ lookupTab=${lkTab1}"
            ;;

    test_extreme ) test1_string="dmimgfilt infile=$INDIR/img.fits outfile="${outfile}" func=extreme mask='box(0,0,5,1)+box(0,0,1,5)' clob+ lookupTab=${lkTab1}"
            ;;

    test_locheq ) test1_string="dmimgfilt infile=$INDIR/img.fits outfile="${outfile}" func=locheq mask='circle(0,0,7)' clob+ lookupTab=${lkTab1}"
            ;;

    test_kuwahara ) test1_string="dmimgfilt infile=$INDIR/img.fits outfile="${outfile}" func=kuwahara mask='box(0,0,5,5)' clob+ lookupTab=${lkTab1}"
            ;;

    test_unsharp ) test1_string="dmimgfilt infile=$INDIR/img.fits outfile="${outfile}" func=unsharp mask='box(0,0,5,5,45)' clob+ lookupTab=${lkTab1}"
            ;;

    ##SAT-39: Test to catch error for new dmRegParse error handling from SL-17
    #Test with empty mask string. Should fail
    test_mask_empty ) test1_string="dmimgfilt infile=$INDIR/img.fits outfile="${outfile}" func=min mask='mask()' clob+ lookupTab=${lkTab1}"
            ;;
    #Test with mask(path/to/mask/file). Should fail for now, until region library method created to set coordinate system outside of region library
    test_mask_path ) test1_string="dmimgfilt infile=$INDIR/img.fits outfile="${outfile}" func=min mask='mask($INDIR/img_mask.fits)' clob+ lookupTab=${lkTab1}"
            ;;


  esac
  if [ ${testid} = "test_mask_empty" ] || [ ${testid} = "test_mask_path" ]; then
      echo $test1_string | tee -a  $LOGFILE
      eval $test1_string
      mystatus=$?
      echo " -- status of call was $mystatus"  >>  $LOGFILE
      echo " -- status of call was $mystatus"  > ${outfile}.status
  else
    echo $test1_string | tee -a  $LOGFILE 
    eval $test1_string  | tee -a  $LOGFILE  2>&1
  fi
 


  ####################################################################
  # check the outputs

  # Init per-test error flag
  mismatch=1

  # if different tests need different kinds of comparisons, use a 
  #  case ${testid} in...  here

  ####################################################################
  # compare
  # !!16
  if [ ${testid} = "test_mask_empty" ] || [ ${testid} = "test_mask_path" ]; then 
     # this test checks for expected errors so diff the exit status of run
     diff $outfile.status $savfile.status  > /dev/null 2>>$LOGFILE
     if  test $? -ne 0 ; then
       echo "ERROR: MISMATCH in $outfile.status" >> $LOGFILE
       mismatch=0
     fi
  else
   dmdiff "${outfile}" $savfile tol=$SAVDIR/tolerance verb=0 > \
         /dev/null 2>>$LOGFILE
   if  test $? -ne 0 ; then
     echo "ERROR: HEADER MISMATCH in ${outfile}" >> $LOGFILE
     mismatch=0
   fi
  fi

  ####################################################################
  # Did we get an error?
  if test $mismatch -eq 0 ; then
    # Yes
    echo "${testid} NOT-OK"
    script_succeeded=1
  else
    # No
    echo "${testid} OK"
  fi

done
# end per-test loop
######################################################################


######################################################################
# report results

# blank line
echo ""

if test $script_succeeded -eq 0; then
    echo "${toolname} : PASS" | tee -a $LOGFILE
else
    echo "${toolname} : FAIL" | tee -a $LOGFILE
fi

echo "log file in ${LOGFILE}"


exit $script_succeeded
