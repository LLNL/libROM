#!/usr/local/bin/perl
##
## Copyright (c) 2013-2016, Lawrence Livermore National Security, LLC.
## Produced at the Lawrence Livermore National Laboratory
## Written by William Arrighi wjarrighi@llnl.gov
## CODE-686965
## All rights reserved.
##
## This file is part of libROM.
## For details, see https://computation.llnl.gov/librom
## Please also read README_BSD_NOTICE.
##
## Redistribution and use in source and binary forms, with or without
## modifications, are permitted provided that the following conditions are met:
##
##    o Redistributions of source code must retain the above copyright notice,
##      this list of conditions and the disclaimer below.
##    o Redistribution in binary form must reproduce the above copyright
##      notice, this list of conditions and the disclaimer (as noted below) in
##      the documentation and/or other materials provided with the
##      distribution.
##    o Neither the name of the LLNS/LLNL nor the names of its contributors may
##      be used to endorse or promote products derived from this software
##      without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
## AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
## IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
## ARE DISCLAIMED.  IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
## LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
## INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
## (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
## LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
## ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
## THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OR SUCH DAMAGE.
##

## Usage: cmp.pl <file1> <file2>
##

$ANAME = shift(@ARGV);
$BNAME = shift(@ARGV);

open(AFILE, "$ANAME") || die "Cannot open input file $ANAME...";
open(BFILE, "$BNAME") || die "Cannot open input file $BNAME...";

while (!eof(AFILE) && !eof(BFILE)) {
   $ALINE = <AFILE>;
   $BLINE = <BFILE>;
   $_ = $ALINE;

   if (!/^(\/\/|c|C|#|##| \*)[ ]*(Release:[\t ]*\$Name|Revision:[\t ]*\$LastChangedRevision|Modified:[\t ]*\$LastChangedDate):[^\$]*\$/o) {
      if ($ALINE ne $BLINE) {
         exit 1;
      }
   }
}

exit 0 if (eof(AFILE) && eof(BFILE));
exit 1;
