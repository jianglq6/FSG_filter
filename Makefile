################################################################################
#  Makefile for FSG
#
################################################################################

# known issues:
#  .h dependency is not included, use make cleanall

#-------------------------------------------------------------------------------
# compiler
#-------------------------------------------------------------------------------

CC     :=  gcc
CXX    :=  g++
NETCDF := /export/apps/gnu-4.8.5/disable-netcdf-4.4.1

#-- 
CFLAGS := -I$(NETCDF)/include -I./lib/ -I./forward/ -I./media/  $(CFLAGS)

#- debug
CFLAGS   := -g -std=c99 $(CFLAGS)
#CPPFLAGS := -g -std=c++11 $(CPPFLAGS)
#- O3
CFLAGS   := -O3 -std=c99 $(CFLAGS)
CPPFLAGS := -O2 -std=c++11 $(CPPFLAGS)

#- static
#LDFLAGS := $(NETCDF)/lib/libnetcdf.a -lm -static $(LDFLAGS)
LDFLAGS := -lm  $(LDFLAGS) $(NETCDF)/lib/libnetcdf.a
#- dynamic
#LDFLAGS := -L$(NETCDF)/lib -lnetcdf -lm $(LDFLAGS)

#-------------------------------------------------------------------------------
# target
#-------------------------------------------------------------------------------

# special vars:
# 	$@ The file name of the target
# 	$< The names of the first prerequisite
#   $^ The names of all the prerequisites 


main_curv_fsg_2d: \
		cJSON.o sacLib.o fdlib_mem.o fdlib_math.o  \
		fd_t.o par_t.o \
		media_utility.o \
		media_layer2model.o \
		media_geometry2d.o\
		media_read_file.o \
		gd_info.o gd_t.o md_t.o wav_t.o \
		blk_t.o \
		bdry_abs.o \
		src_t.o \
		filter_fsg.o \
		sv_eq1st_curv_fsg_el_aniso.o \
		main_curv_fsg_2d.o
	$(CXX) -o $@ $^ $(LDFLAGS)


media_geometry2d.o: media/media_geometry2d.cpp 
	${CXX} -c -o $@ $(CPPFLAGS) $<
media_utility.o: media/media_utility.cpp 
	${CXX} -c -o $@ $(CPPFLAGS) $<
media_layer2model.o: media/media_layer2model.cpp
	${CXX} -c -o $@ $(CPPFLAGS) $<
media_read_file.o: media/media_read_file.cpp
	${CXX} -c -o $@ $(CPPFLAGS) $<
cJSON.o: lib/cJSON.c
	${CC} -c -o $@ $(CFLAGS) $<
sacLib.o: lib/sacLib.c
	${CC} -c -o $@ $(CFLAGS) $<
fdlib_mem.o: lib/fdlib_mem.c
	${CC} -c -o $@ $(CFLAGS) $<
fdlib_math.o: lib/fdlib_math.c
	${CC} -c -o $@ $(CFLAGS) $<
fd_t.o: forward/fd_t.c
	${CC} -c -o $@ $(CFLAGS) $<
par_t.o: forward/par_t.c
	${CC} -c -o $@ $(CFLAGS) $<
gd_info.o: forward/gd_info.c
	${CC} -c -o $@ $(CFLAGS) $<
gd_t.o: forward/gd_t.c
	${CC} -c -o $@ $(CFLAGS) $<
md_t.o: forward/md_t.c
	${CC} -c -o $@ $(CFLAGS) $<
wav_t.o: forward/wav_t.c
	${CC} -c -o $@ $(CFLAGS) $<
blk_t.o: forward/blk_t.c
	${CC} -c -o $@ $(CFLAGS) $<
sv_eq1st_curv_fsg_el_aniso.o:   forward/sv_eq1st_curv_fsg_el_aniso.c
	${CC} -c -o $@ $(CFLAGS) $<
bdry_abs.o:   forward/bdry_abs.c
	${CC} -c -o $@ $(CFLAGS) $<
src_t.o:   forward/src_t.c
	${CC} -c -o $@ $(CFLAGS) $<
filter_fsg.o:   forward/filter_fsg.c
	${CC} -c -o $@ $(CFLAGS) $<

main_curv_fsg_2d.o: forward/main_curv_fsg_2d.c
	${CC} -c -o $@ $(CFLAGS) $<

cleanexe:
	rm -f main_curv_fsg_2d

cleanobj:
	rm -f *.o
cleanall: cleanexe cleanobj
	echo "clean all"
distclean: cleanexe cleanobj
	echo "clean all"
