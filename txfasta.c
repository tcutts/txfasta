/*
	TPatterns

	Copyright (C) T J R Cutts 1998-2002

	timc@chiark.greenend.org.uk

	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tplib.h"

#define LINELENGTH 60

void txfasta(void)
{
	char *s, *t, *p;
	char *frm[6];
	char buf[1024];
	char title[1024];
	int sz = 8192;
	int i, n, len;
	int lengths[6];

	s = (char *)malloc(sz);
	t = (char *)malloc(sz);
	frm[0] = (char *)malloc(sz);
	frm[1] = (char *)malloc(sz);
	frm[2] = (char *)malloc(sz);
	frm[3] = (char *)malloc(sz);
	frm[4] = (char *)malloc(sz);
	frm[5] = (char *)malloc(sz);

	if (!frm[0] || !frm[1] || !frm[2] ||
		!frm[3] || !frm[4] || !frm[5] ||
		!s || !t)
	{
		tp_error(TPERR_MEM);
		exit(1);
	}

	*s = '\0';
	len = 0;

	while (!feof(stdin))
	{
		p = fgets(buf, sizeof(buf), stdin);
		if (p)
		{
			if (*p == '>')
			{
				/* We have found a comment line */
				if (*s != '\0')
				{
					s[len] = '\0';
					/* We now have a complete sequence, so
					   process it */
					translate(s, frm, lengths);
					for (n = 0; n < 6; n++)
					{
						printf("%s (Frame %d)\n",
							   title,
							   n);
						i = LINELENGTH;
						p = frm[n];
						while (i < lengths[n])
						{
							fwrite(p, sizeof(char), LINELENGTH, stdout);
							fputc('\n', stdout);
							i += LINELENGTH;
							p += LINELENGTH;
						}
						printf("%s\n", p);
					}
					len = 0;
					*s = '\0';
				}
				buf[strlen(buf) - 1] = '\0';
				strcpy(title, buf);
			}
			else
			{
				/* This is a sequence line, so add it to s */
				n = strlen(p) - 1;
				if ((len + n + 3) > sz)
				{
					sz <<= 1;
					s = (char *)realloc(s, sz);
					t = (char *)realloc(t, sz);
					frm[0] = (char *)realloc(frm[0], sz);
					frm[1] = (char *)realloc(frm[1], sz);
					frm[2] = (char *)realloc(frm[2], sz);
					frm[3] = (char *)realloc(frm[3], sz);
					frm[4] = (char *)realloc(frm[4], sz);
					frm[5] = (char *)realloc(frm[5], sz);
					if (!frm[0] || !frm[1] || !frm[2] ||
						!frm[3] || !frm[4] || !frm[5] ||
						!s || !t)
					{
						tp_error(TPERR_MEM);
						exit(1);
					}
				}
				memcpy(&s[len], buf, n);
				len += n;
			}
		}
	}

	if (*s != '\0')
	{
		s[len] = '\0';
		/* We now have a complete sequence, so
	   process it */
		translate(s, frm, lengths);
		for (n = 0; n < 6; n++)
		{
			printf("%s (Frame %d)\n",
				   title,
				   n);
			i = LINELENGTH;
			p = frm[n];
			while (i < lengths[n])
			{
				fwrite(p, sizeof(char), LINELENGTH, stdout);
				fputc('\n', stdout);
				i += LINELENGTH;
				p += LINELENGTH;
			}
			printf("%s\n", p);
		}
		len = 0;
		*s = '\0';
	}
	free(frm[0]);
	free(frm[1]);
	free(frm[2]);
	free(frm[3]);
	free(frm[4]);
	free(frm[5]);
	free(t);
	free(s);
}

/* usage() prints a usage message and returns a non-zero
   exit status */

void usage(char *a, char *b)
{
	fprintf(stderr,
			"%s\n\nUsage:\n\t%s [-m translation] < infile > outfile\n",
			b, a);
	exit(1);
}

/* options() processes the command line arguments,
   returning the index of the first argument that
   is not an option */

int options(int argc, char *argv[])
{
	char *opts = "m:";
	extern char *optarg;
	extern int optind;
	int c;

	while ((c = getopt(argc, argv, opts)) != -1)
	{
		switch (c)
		{
		case 'm':
			compilemx(argv[optind - 1]);
			break;
		default:
			usage(argv[0], "Unknown option");
		}
	}

	return optind;
}

int main(int argc, char *argv[])
{
	options(argc, argv);
	initbasebits();
	make_revmatrix();
	txfasta();
	return 0;
}
