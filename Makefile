PREFIX=/usr/local
LIBDIR=$(PREFIX)/lib
BINDIR=$(PREFIX)/bin

CFLAGS=-O3 -DTFLIB=\"$(LIBDIR)\"

SOURCES=txfasta.c errors.c translate.c

all: txfasta

txfasta: $(SOURCES)
	$(CC) -o $@ $(CFLAGS) $(SOURCES)

test: all
	./txfasta < test/test.fa > test-out.fa
	cmp test-out.fa test/testp.fa

install: txfasta
	mkdir -p $(LIBDIR)
	mkdir -p $(BINDIR)
	cp txfasta $(BINDIR)

clean:
	rm -f core *.o txfasta *~ test-out.fa
