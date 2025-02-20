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

#include "tplib.h"

#include <stdio.h>

const char *messages[] = {
  "Success", /* TF_SUCCESS */
  "Memory allocation failure", /* TF_MEM */
  "Premature end of matrix file", /* TF_MXEOF */
  "Could not open file for output"
};

void tp_error(tperror code)
{
  fputs(messages[code], stderr);
}
