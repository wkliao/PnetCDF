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

noinst_LTLIBRARIES = libromio.la

H_SRCS = adio.h

C_SRCS = ad_read.c \
         ad_write.c \
         ad_open.c \
         ad_close.c \
         ad_fstype.c \
         ad_aggregate.c \
         ad_read_str.c \
         ad_read_coll.c \
         ad_read_str_naive.c \
         ad_write_coll.c \
         ad_write_str.c \
         ad_write_str_naive.c \
         ad_subarray.c \
         ad_darray.c \
         ad_utils.c \
         ad_lustre_open.c \
         ad_lustre_wrcoll.c \
         ad_lustre_wrstr.c \
         ad_iscontig.c \
         ad_flatten.c \
         ad_lock.c \
         ad_set_size.c \
         ad_sync.c \
         ad_delete.c \
         ad_set_view.c \
         ad_hints.c


libromio_la_SOURCES = $(C_SRCS) $(H_SRCS)

CLEANFILES = core core.* *.gcda *.gcno *.gcov gmon.out

tests-local: all

