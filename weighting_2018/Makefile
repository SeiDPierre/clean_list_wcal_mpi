SRC 	=	clean_weight listing_hdr wcal wcal_mpi listing_hdr_mpi

.PHONY : clean $(SRC)

all	:	$(SRC)

$(SRC)	:
	make -C $@

clean	:	$(patsubst %,clean.%,$(SRC))
clean.% :
	make -C $* clean
