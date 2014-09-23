#!/usr/bin/perl
# merge kaeb and gasch data

BEGIN { unshift(@INC,"/home/hqin/lib/perl/");   }

use strict; use warnings; use Util;

my $k_file  = "kaeberlein.tab";
my $g_file  = "gasch00.tab";

my $out_file = "k.and.gasch.tab";
my $debug = 1;

#temp variables
my (@lines, $line, @all_lines ) = ();

# read in k_file
my %k_orfs = ();
open (IN, "<$k_file"); 
while ( $line = <IN> ) {
 chomp $line;
 my ( $orf, $name, @rests ) = split( /\t/, $line ) ;
 $k_orfs{ $orf } = $name;
}
close (IN);

# take gash records
my %g_lines = ();
open (IN, "<$g_file");
while( my $line= <IN> ) {
   chomp $line;
   my ( $orf, @tokens) = split( /[\s\t]+/, $line);
   #print OUT "$name_2_orf{$id}\t$line\n";
   if (exists $k_orfs{ $orf } ) {
      $g_lines { $orf  } = $line;
   }
}
close (IN);

#merge the file
open (IN, "<$k_file"); 
open (OUT, ">$out_file");
while ( $line = <IN> ) {
 chomp $line;
 my ( $orf, $name, @rests ) = split( /\t/, $line ) ;
 print OUT "$line\t$g_lines{$orf}\n";
}
close (IN);




#
# DONE


