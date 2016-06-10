##############################################################################


MK_TOP = /export/ciao_from_source/ciao-4.8/src
KJG = /export/ciao

include $(MK_TOP)/Makefile.master
include $(MK_TOP)/include/Makefile.scidev

EXEC              = dmimgfilt
LIB_FILES         =
PAR_FILES         = dmimgfilt.par
INC_FILES         =
XML_FILES         = dmimgfilt.xml

SRCS	= dmimgfilt.c t_dmimgfilt.c filters.c

OBJS	= $(SRCS:.c=.o)

MAKETEST_SCRIPT = dmimgfilt.t


include $(MK_TOP)/Makefile.all

#-----------------------------------------------------------------------
# 			MAKEFILE DEPENDENCIES	
#-----------------------------------------------------------------------

$(EXEC): $(OBJS)
	$(LINK)
	@echo

kjg: $(EXEC)
	/bin/cp -f $(EXEC) $(KJG)/binexe/
	/bin/cp -f $(KJG)/bin/dmlist $(KJG)/bin/$(EXEC)
	/bin/cp -f $(PAR_FILES) $(KJG)/param/$(PAR_FILES)

announce1:
	@echo "   /----------------------------------------------------------\ "
	@echo "   |             Building dmimgfilt DM host tool           | "
	@echo "   \----------------------------------------------------------/ "

