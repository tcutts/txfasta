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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tplib.h"

/* The following matrix encodes the Universal genetic code.  During
   execution any matrix supplied via the -m option will overwrite it */

txmatrix matrix = {{{"KNKN", "TTTT", "RSRS", "IIMI"},
                    {"QHQH", "PPPP", "RRRR", "LLLL"},
                    {"EDED", "AAAA", "GGGG", "VVVV"},
                    {"*Y*Y", "SSSS", "*CWC", "LFLF"}},
                   {{"XXXX", "XXXX", "XXXX", "XXXX"},
                    {"XXXX", "XXXX", "XXXX", "XXXX"},
                    {"XXXX", "XXXX", "XXXX", "XXXX"},
                    {"XXXX", "XXXX", "XXXX", "XXXX"}},
                   {{"XXXX", "XXXX", "XXXX", "XXXX"},
                    {"XXXX", "XXXX", "XXXX", "XXXX"},
                    {"XXXX", "XXXX", "XXXX", "XXXX"},
                    {"XXXX", "XXXX", "XXXX", "XXXX"}},
                   {{"XXXX", "XXXX", "XXXX", "XXXX"},
                    {"XXXX", "XXXX", "XXXX", "XXXX"},
                    {"XXXX", "XXXX", "XXXX", "XXXX"},
                    {"XXXX", "XXXX", "XXXX", "XXXX"}}};

txmatrix revmatrix;

/* make_revmatrix() creates the reverse complemented matrix by
   swapping bits 4 and 5 with bits 0 and 1, inverting the bits, and
   then keeping only bits 0-5 */

void make_revmatrix(void)
{
  int n, r;

  for (n = 0; n < 64; n++)
  {
    r = (~((n & 0x0C) |
           ((n & 0x03) << 4) |
           ((n >> 4) & 0x03))) &
        0x3F;

    ((char *)revmatrix)[r] = ((char *)matrix)[n];
  }

  for (n = 64; n < 256; n++)
    ((char *)revmatrix)[n] = 'X';
}

/* compilemx() reads a file of codon -> amino acid translations and
   places them in matrix, above, for later use by the translate()
   function */

int compilemx(char *filename)
{
  char buf[80];
  char path[PATH_MAX];

  FILE *f;
  int res;
  char ta, tb, tc, tx, *r;

  if ((f = fopen(filename, "r")) == NULL)
  {
    /* Don't bother with the environment if an absolute path was
 given */
    if (filename[0] != '/')
    {
      if ((r = getenv("TP_TXMDIR")) != NULL)
      {
        sprintf(path, "%s/%s", r, filename);
        f = fopen(path, "r");
      }
      if (f == NULL)
      {
        sprintf(path, "%s/%s", TFLIB, filename);
        f = fopen(path, "r");
      }
    }
  }

  if (f == NULL)
  {
    perror("Could not open translation matrix file:");
    exit(1);
  }

  for (int a = 0; a < 4; a++)
  {
    for (int b = 0; b < 4; b++)
    {
      for (int c = 0; c < 4; c++)
      {
        do
        {
          r = fgets(buf, sizeof(buf), f);
          if (r == NULL)
          {
            /* Premature end of file encountered */
            tp_error(TPERR_MXEOF);
            exit(1);
          }
          res = sscanf(buf,
                       "%c%c%c %c",
                       &ta, &tb, &tc, &tx);
        } while (res != 4);

        matrix[0][a][b][c] = tx;
      }
    }
  }

  fclose(f);
  return 0;
}

/* initbasebits() initialises the lookup table for the translate()
   function */

int basebits[256];

#ifdef NVSN
int comphash[256];
#endif

void initbasebits(void)
{

  /* First of all the hash for N->P translation */

  memset(basebits, 255, sizeof(int) * 256);
  basebits['A'] = 0;
  basebits['C'] = 1;
  basebits['G'] = 2;
  basebits['T'] = 3;
  basebits['U'] = 3;
  basebits['a'] = 0;
  basebits['c'] = 1;
  basebits['g'] = 2;
  basebits['t'] = 3;
  basebits['u'] = 3;

#ifdef NVSN
  /* And now the hash for reverse complementation */
  memset(comphash, 'x', sizeof(int) * 256);
  comphash['A'] = 'T';
  comphash['B'] = 'V';
  comphash['C'] = 'G';
  comphash['D'] = 'H';
  comphash['G'] = 'C';
  comphash['H'] = 'D';
  comphash['K'] = 'M';
  comphash['M'] = 'K';
  comphash['N'] = 'N';
  comphash['R'] = 'Y';
  comphash['S'] = 'S';
  comphash['T'] = 'A';
  comphash['U'] = 'A';
  comphash['V'] = 'B';
  comphash['X'] = 'X';
  comphash['Y'] = 'Y';
  comphash['.'] = '.';
  comphash['~'] = '~';
#endif
}

/* translate() a nucleic acid sequence in all six reading frames
   simultaneously.  output_aa_sequence should be an array of six (char *) pointers,
   all of which must be large enough to hold the translated sequence.
   No sanity checking is performed here. output_aa_seq_length is an array of six
   integers, which will take the lengths of the translated sequences.

   This code is ugly as hell, but it is *very* fast, for the following reasons:

   Threefold loop unrolling, and very fast bit shift and comparison
   operators results in only four conditionals executed per loop
   iteration, but generate 6 translated amino acids.  This keeps the CPU
   pipeline as full as possible.  The lookup table is small enough (1K)
   to remain in the L1 cache of pretty much any CPU, too.

   Thanks to Peter Benie for help optimising the code */

void translate(char *input_dna_sequence, char **output_aa_sequence, int *output_aa_seq_length)
{
  int input_sequence_length, codon_remainder, index;
  char *p, *input_sequence_end, *rev_frame_3, *rev_frame_4, *rev_frame_5;
  char *r0 = output_aa_sequence[0];
  char *r1 = output_aa_sequence[1];
  char *r2 = output_aa_sequence[2];

  input_sequence_length = strlen(input_dna_sequence);
  input_sequence_end = &input_dna_sequence[input_sequence_length] - 5;

  /* How much of a part codon do we have at the end */
  codon_remainder = input_sequence_length % 3;

  /* How many residues long will the frame 0 sequence be? */
  input_sequence_length /= 3;

  output_aa_seq_length[0] = output_aa_seq_length[1] = output_aa_seq_length[2] = output_aa_seq_length[3] = output_aa_seq_length[4] = output_aa_seq_length[5] = input_sequence_length;

  /* Position the reverse reading frame pointers where the
     protein sequences will end */

  rev_frame_3 = &output_aa_sequence[3][input_sequence_length];
  rev_frame_4 = &output_aa_sequence[4][input_sequence_length];
  rev_frame_5 = &output_aa_sequence[5][input_sequence_length];

  /* If there aren't two bases at the end, the translations
     will not all be the same length; one or two of them
     will be one residue shorter */

  if (codon_remainder != 2)
  {
    rev_frame_5--;
    output_aa_seq_length[5]--;
    output_aa_seq_length[2]--;
    if (codon_remainder == 0)
    {
      rev_frame_4--;
      output_aa_seq_length[4]--;
      output_aa_seq_length[1]--;
    }
  }

  *rev_frame_3 = *rev_frame_4 = *rev_frame_5 = '\0';
  rev_frame_3--;
  rev_frame_4--;
  rev_frame_5--;

  index = (basebits[input_dna_sequence[0]] << 2) | basebits[input_dna_sequence[1]];

  /* OK, that's all the preamble, now lets get to the
   * main loop! */

  for (p = input_dna_sequence; p <= input_sequence_end; p += 3)
  {

    /* Frame 0/3 */
    index = ((index << 2) & 0x3C) | basebits[p[2]];

    *r0++ = ((char *)matrix)[index];
    *rev_frame_3-- = ((char *)revmatrix)[index];

    /* Frame 1/4 */

    index = ((index << 2) & 0x3C) | basebits[p[3]];

    *r1++ = ((char *)matrix)[index];
    *rev_frame_4-- = ((char *)revmatrix)[index];

    /* Frame 2/5 */

    index = ((index << 2) & 0x3C) | basebits[p[4]];

    *r2++ = ((char *)matrix)[index];
    *rev_frame_5-- = ((char *)revmatrix)[index];
  }

  /* We may, at this point, have one or two codons still untranslated,
     so we need to deal with those */

  if (codon_remainder < 2)
  {
    index = ((index << 2) & 0x3C) | basebits[p[2]];

    *r0++ = ((char *)matrix)[index];
    *rev_frame_3-- = ((char *)revmatrix)[index];

    if (codon_remainder == 1)
    {
      index = ((index << 2) & 0x3C) | basebits[p[3]];

      *r1++ = ((char *)matrix)[index];
      *rev_frame_4-- = ((char *)revmatrix)[index];
    }
  }

  /* And that's it... */

  *r0 = *r1 = *r2 = '\0';
}

#ifdef NVSN
/* rev_comp() produces the reverse complement of a nucleic acid
   sequence.  It is only used when comparing a nucleotide pattern
   against a nucleotide db */

void rev_comp(char *input_dna_sequence, char *out)
{
  register char *p = input_dna_sequence;
  register char *q;

  q = &out[strlen(p)];
  *q-- = '\0';

  while (q >= out)
    *q-- = comphash[*p++];
}
#endif
