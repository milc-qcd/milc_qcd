#  Makefile to look at the propagating form factors
#
#   make dump_3pt.x  --- program to dump all the propagating heavy --> light three
#                       point functions to standard output.
#
#   make dump_2pt.x  --- writes to standard output the various two point functions
#                        required in the form factor project.
#

MAKEFILE = Makefile

#COMPILER = gcc  -pg
COMPILER = cc
#COMPILER = gcc
FLAGS =  -lm
#CFLAGS = 
#CFLAGS = -g -DDEBUG
#CFLAGS = -Wall -pedantic
#CFLAGS = -g
#CFLAGS = -O -DSHORT32     # T3E
CFLAGS = -O                # Others  

LISTA =    \
        byterev_array.c\
	bsd_sum.c\
        dump_3pt.c\
        main_dump_3pt.c\
        names_of_three_operators.c\
        read_3pt.c\
        titles_hl_propagating.c

LISTB= \
       byterev_array.c\
       bsd_sum.c\
       dump_2pt.c\
       main_dump_2pt.c\
       names_of_two_operators.c\
       read_2pt.c\
       titles_two_pt.c

LISTC =    \
        bsd_sum.c\
        byterev_array.c\
        calc_jack_samples.c\
        get_input_param.c\
        jack_routines.c\
        main_propagating_matrix2.c\
        names_of_three_operators.c\
        names_of_two_operators.c\
        select_3pt.c\
        select_2pt.c\
        write_matrix_element.c 


LISTD= \
       byterev_array.c\
       bsd_sum.c\
       dump_2pt.c\
       get_select_input_param.c\
       main_select_2pt.c\
       names_of_two_operators.c\
       select_2pt.c\
       titles_two_pt.c

LISTE= \
       byterev_array.c\
       bsd_sum.c\
       dump_3pt.c\
       get_select_input_param.c\
       main_select_3pt.c\
       names_of_three_operators.c\
       select_3pt.c\
       titles_hl_propagating.c

LISTF= \
       byterev_array.c\
       bsd_sum.c\
       dump_2pt.c\
       main_split_3pt.c\
       names_of_two_operators.c\
       read_3pt.c\
       write_prop_3pt_onemom.c \
       titles_hl_propagating.c

LISTG= \
       byterev_array.c\
       bsd_sum.c\
       dump_2pt.c\
       main_merge_3pt.c\
       names_of_two_operators.c\
       merge_3pt.c\
       write_3pt.c\
       titles_hl_propagating.c

LISTH= \
       byterev_array.c\
       bsd_sum.c\
       main_norm_2pt.c\
       names_of_two_operators.c\
       read_2pt.c\
       write_2pt.c\
       titles_two_pt.c

LISTI= \
       byterev_array.c\
       bsd_sum.c\
       dump_2pt.c\
       main_split_3pt2.c\
       names_of_two_operators.c\
       read_3pt.c\
       write_prop_3pt_onemom.c \
       titles_hl_propagating.c

LISTJ= \
       byterev_array.c\
       bsd_sum.c\
       get_multiselect_input_param.c\
       main_multiselect_3pt.c\
       multiselect_3pt.c\
       names_of_three_operators.c

LISTK= \
       byterev_array.c\
       bsd_sum.c\
       get_multiselect_input_param.c\
       main_multiselect_2pt.c\
       multiselect_2pt.c\
       names_of_two_operators.c

LISTL= \
       byterev_array.c\
       bsd_sum.c\
       dump_2pt.c\
       main_split_spect_2pt.c\
       names_of_two_operators.c\
       read_2pt.c\
       write_prop_2pt_onespect.c \
       titles_two_pt.c

LISTM= \
       byterev_array.c\
       bsd_sum.c\
       dump_3pt.c\
       main_split_spect_3pt2.c\
       names_of_three_operators.c\
       read_3pt.c\
       write_prop_3pt_onespect.c \
       titles_hl_propagating.c

LISTN= \
       byterev_array.c\
       bsd_sum.c\
       dump_3pt.c\
       main_split_spectmom_3pt.c\
       names_of_three_operators.c\
       read_3pt.c\
       write_prop_3pt_onemom.c \
       titles_hl_propagating.c

INCS = \
	read_hl_form.h\
	prop_form_utilities.h


OBJSA =  $(LISTA:.c=.o)
OBJSB =  $(LISTB:.c=.o)
OBJSC =  $(LISTC:.c=.o)
OBJSD =  $(LISTD:.c=.o)
OBJSE =  $(LISTE:.c=.o)
OBJSF =  $(LISTF:.c=.o)
OBJSG =  $(LISTG:.c=.o)
OBJSH =  $(LISTH:.c=.o)
OBJSI =  $(LISTI:.c=.o)
OBJSJ =  $(LISTJ:.c=.o)
OBJSK =  $(LISTK:.c=.o)
OBJSL =  $(LISTL:.c=.o)
OBJSM =  $(LISTM:.c=.o)
OBJSN =  $(LISTN:.c=.o)


.SUFFIXES:

.SUFFIXES:	.o .c

.c.o:
	$(COMPILER) -c  $< $(CFLAGS)


$(OBJSA): $(INCS) $(MAKEFILE)
$(OBJSB): $(INCS) $(MAKEFILE)
$(OBJSC): $(INCS) $(MAKEFILE)
$(OBJSD): $(INCS) $(MAKEFILE)
$(OBJSE): $(INCS) $(MAKEFILE)
$(OBJSF): $(INCS) $(MAKEFILE)
$(OBJSG): $(INCS) $(MAKEFILE)
$(OBJSH): $(INCS) $(MAKEFILE)
$(OBJSI): $(INCS) $(MAKEFILE)
$(OBJSJ): $(INCS) $(MAKEFILE)
$(OBJSK): $(INCS) $(MAKEFILE)
$(OBJSL): $(INCS) $(MAKEFILE)
$(OBJSM): $(INCS) $(MAKEFILE)
$(OBJSN): $(INCS) $(MAKEFILE)

dump_3pt.x : $(LISTA) $(OBJSA) $(INCS) $(MAKEFILE) 
	$(COMPILER) -o $@   $(OBJSA)   $(FLAGS)

dump_2pt.x : $(LISTB) $(OBJSB) $(INCS) $(MAKEFILE) 
	$(COMPILER) -o $@   $(OBJSB)   $(FLAGS)

mat_el.x : $(LISTC) $(OBJSC) $(INCS) $(MAKEFILE) 
	$(COMPILER) -o $@   $(OBJSC)   $(FLAGS)

select_2pt.x : $(LISTD) $(OBJSD) $(INCS) $(MAKEFILE) 
	$(COMPILER) -o $@   $(OBJSD)   $(FLAGS)

select_3pt.x : $(LISTE) $(OBJSE) $(INCS) $(MAKEFILE) 
	$(COMPILER) -o $@   $(OBJSE)   $(FLAGS)

split_3pt.x : $(LISTF) $(OBJSF) $(INCS) $(MAKEFILE) 
	$(COMPILER) -o $@   $(OBJSF)   $(FLAGS)

merge_3pt.x : $(LISTG) $(OBJSG) $(INCS) $(MAKEFILE) 
	$(COMPILER) -o $@   $(OBJSG)   $(FLAGS)

renorm_2pt.x : $(LISTH) $(OBJSH) $(INCS) $(MAKEFILE) 
	$(COMPILER) -o $@   $(OBJSH)   $(FLAGS)

split_3pt2.x: $(LISTI) $(OBJSI) $(INCS) $(MAKEFILE) 
	$(COMPILER) -o $@   $(OBJSI)   $(FLAGS)

multiselect_3pt.x : $(LISTJ) $(OBJSJ) $(INCS) $(MAKEFILE) 
	$(COMPILER) -o $@   $(OBJSJ)   $(FLAGS)

multiselect_2pt.x : $(LISTK) $(OBJSK) $(INCS) $(MAKEFILE) 
	$(COMPILER) -o $@   $(OBJSK)   $(FLAGS)

split_2pt_spect.x : $(LISTL) $(OBJSL) $(INCS) $(MAKEFILE) 
	$(COMPILER) -o $@   $(OBJSL)   $(FLAGS)

split_3pt_spect2.x : $(LISTM) $(OBJSM) $(INCS) $(MAKEFILE) 
	$(COMPILER) -o $@   $(OBJSM)   $(FLAGS)

split_3pt_spectmom.x : $(LISTN) $(OBJSN) $(INCS) $(MAKEFILE) 
	$(COMPILER) -o $@   $(OBJSN)   $(FLAGS)
