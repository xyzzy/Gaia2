ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST = LICENSE ChangeLog.md

MAINTAINERCLEANFILES = 
DISTCLEANFILES =
CLEANFILES = 

bin_PROGRAMS = gaia2 merge simulTame4

AM_CXXFLAGS = -std=c++98

gaia2_LDADD = /usr/lib/libgd.a /usr/lib/libpng.a /usr/lib/libjpeg.a /usr/lib/libz.a -lpthread
gaia2_SOURCES = gaia2.cc

merge_LDADD = /usr/lib/libgd.a /usr/lib/libpng.a /usr/lib/libjpeg.a /usr/lib/libz.a -lpthread
merge_SOURCES = merge.cc

simulTame4_LDADD = /usr/lib/libgd.a /usr/lib/libpng.a /usr/lib/libjpeg.a /usr/lib/libz.a -lpthread
simulTame4_SOURCES = simulTame4.cc
