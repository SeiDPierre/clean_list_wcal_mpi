include ../flags.mk

OBJ	= 	$(O)/listing_hdr_az.o \
		$(O)/subs_acal_A30_anis_cl.o \
		$(O)/legndr.o \
		$(O)/excit_q1.o \
		$(O)/midpoint.o \
		$(O)/geteigys_ellrot.o \
		$(O)/util.o

LIB 	=	$(L)/mathlib.$(ARCH).a \
		$(L)/sphere.$(ARCH).a \
		$(L)/iolib.$(ARCH).a \
		$(L)/prem.$(ARCH).a \
		$(L)/swaprecord.$(ARCH).a -lm

BIN	=	$(B)/listing_hdr_az_mpi.x

all	:	$(BIN)

.PHONY: clean
clean:
	rm -rf $(O) $(B)

# Main executable
$(BIN)	:	$(O) $(B) $(OBJ)
	$(FC) $(OBJ) $(LIB) -o $(@)

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
