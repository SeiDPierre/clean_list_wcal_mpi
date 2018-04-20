include ../flags.mk

OBJ	=	$(O)/wcal_az_mpi.o \
		$(O)/subs_acal_A30_anis_cl.o \
		$(O)/legndr.o \
		$(O)/pathprojectn1.o \
		$(O)/excit_q1.o \
		$(O)/geteigys_ellrot.o

LIB	=	$(L)/mathlib.$(ARCH).a \
		$(L)/sphere.$(ARCH).a \
		$(L)/iolib.$(ARCH).a \
		$(L)/prem.$(ARCH).a \
		$(L)/swaprecord.$(ARCH).a -lm

BIN     =       $(B)/wcal_az_mpi.x

all     :       $(BIN)

.PHONY: clean
clean:
	rm -rf $(O) $(B)

# Main executable
$(BIN)	:	$(O) $(B) $(OBJ)
	$(FC) $(OBJ) $(LIB) -o $(@)

# Generic rules, first acal subs, then main code
$(O)/%.o:	$(DIR)/%.f
	$(FC) $(CFLAGS) $< -o $@

$(O)/%.o:	$(DIR)/%.c
	$(CC) -std=c99 $(CFLAGS) $< -o $@

$(O)/%.o:	%.c
	$(CC) -std=c99 $(CFLAGS) $< -o $@

$(O)	:
	mkdir -p $(O)
$(B)	:
	mkdir -p $(B)
