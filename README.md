

This is a fork of the CIAO 4.6 code.

It adds 3 new statistics:

  - ridge: same as 'peak' but uses < rather than <= so values if they
    are the max value are included.
    
  - plain: same as 'valley' but uses > rather than >= so all values 
    if they are the min value are included.

  - qxx : Use the environment variable QXX to extract an arbitary
    quantile.  If the env. variable is not set then the 
    output will all be NaN.

