#!/usr/bin/perl


BEGIN { unshift(@INC,"/home/hqin/lib/perl/");   }

use strict; use warnings; use Util;

my $in_file     = "kaeberlein.tmp.tab";
my $in_id_file  = "Orf2Name_SGD.tab";

my $out_file = "kaeberlein.tab";
my $debug = 1;

#temp variables
my (@lines, $line, @all_lines ) = ();

# set up the lookup table
my %name_2_orf = ();
open (IN, "<$in_id_file");
while ( $line = <IN> ) {
 chomp $line;
 my ( $orf, $name, @rests ) = split( /\t/, $line ) ;
 $name_2_orf{ $name } = $orf;
}
close (IN);
showHashTable(\%name_2_orf);

open (OUT, ">$out_file");
open (IN, "<$in_file");
while( my $line= <IN> ) {
   chomp $line;
   my @tokens = split( /[\s\t]+/, $line);
   print @tokens;
   my $id = $tokens[0];
   print OUT "$name_2_orf{$id}\t$line\n";
}
close (OUT);
close (IN);

#
# DONE


