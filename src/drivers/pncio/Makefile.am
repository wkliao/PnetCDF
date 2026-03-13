#
# Copyright (C) 2025, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# @configure_input@

SUFFIXES = .a .o .c .m4 .h

AM_CPPFLAGS  = -I${top_srcdir}/src/include
AM_CPPFLAGS += -I${top_builddir}/src/include
AM_CPPFLAGS += -I${top_srcdir}/src/drivers/include
AM_CPPFLAGS += -I${top_builddir}/src/drivers/include

if PNETCDF_DEBUG
   AM_CPPFLAGS += -DPNETCDF_DEBUG
endif

noinst_LTLIBRARIES = libpncio.la

H_SRCS = pncio.h

C_SRCS = pncio_open.c \
         pncio_close.c \
         pncio_fstype.c \
         pncio_aggregate.c \
         pncio_read_coll.c \
         pncio_read_indep.c \
         pncio_write_coll.c \
         pncio_write_indep.c \
		 pncio_ufs_rd_indep.c \
		 pncio_ufs_rd_coll.c \
		 pncio_ufs_wr_indep.c \
		 pncio_ufs_wr_coll.c \
         pncio_lustre_open.c \
         pncio_lustre_wr_coll.c \
         pncio_lustre_wr_indep.c \
         pncio_lock.c \
         pncio_set_size.c \
         pncio_sync.c \
         pncio_delete.c \
         pncio_set_view.c \
         pncio_hints.c \
         pncio_utils.c


libpncio_la_SOURCES = $(C_SRCS) $(H_SRCS)

CLEANFILES = core core.* *.gcda *.gcno *.gcov gmon.out

tests-local: all

