include ../flags.mk
OBJ	=	$(O)/wcal_clean_az.o

LIB	=	$(L)/mathlib.$(ARCH).a \
		$(L)/sphere.$(ARCH).a \
		$(L)/iolib.$(ARCH).a \
		$(L)/prem.$(ARCH).a -lm

BIN	=	$(B)/wcal_clean_az.x
# Targets and rules
all	:	$(BIN)

# Main executable
$(BIN)	:	$(O) $(B) $(OBJ)
	$(FC) $(OBJ) $(LIB) -o $(@)

.PHONY: clean
clean:
	rm -rf $(B) $(O)

# Rules
$(O)/%.o: %.c
	$(CC) -std=c99 $(CFLAGS) $< -o $@

# Generic where to build
$O :
	mkdir -p $O
$B :
	mkdir -p $B
