tool = dmimgfilt

binexedir = $(prefix)/binexe
paramdir = $(prefix)/param
ahelpdir = $(prefix)/share/doc/xml

binexe_PROGRAMS = $(tool)

dmimgfilt_SOURCES = dmimgfilt.c t_dmimgfilt.c
dmimgfilt_CPPFLAGS = $(CIAO_CFLAGS)
dmimgfilt_LDADD = $(CIAO_LIBS)
dmimgfilt_LINK = $(CXX) -o $@ -Wl,-rpath,$(prefix)/lib -Wl,-rpath,$(prefix)/ots/lib 

# problems with libstdc++ and cxcparam on Linux :(
if LINUX
dmimgfilt_LDADD += -L$(prefix)/ots/lib -lstdc++
endif


dist_param_DATA = $(tool).par
dist_ahelp_DATA = $(tool).xml

install-data-hook:
	chmod a-w $(paramdir)/$(dist_param_DATA)
	chmod a-w $(ahelpdir)/$(dist_ahelp_DATA)
