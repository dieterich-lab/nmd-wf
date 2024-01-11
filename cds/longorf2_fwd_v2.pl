#!/usr/bin/perl -w
# longorf.pl v0208020920
# (c) Dan Kortschak 2002

use vars qw($USAGE);

use strict;
use Getopt::Long;
use Bio::SeqIO;

$USAGE = "longorf [--help] [--notstrict] [--verbose] [--graph] [--width printwidth] [--format seqformat] --input seqfile\n";

my ($sequencefile,$sequenceformat,$notstrict,$startcodon,$graph,$verb,$printwidth,$help) =
   (undef,        'fasta',        undef,      undef, undef,1,50,         undef);

&GetOptions('input|i=s'              => \$sequencefile,
            'format|f=s'             => \$sequenceformat,
            'notstrict|n'            => \$notstrict,
            'startcodon|s=s'            => \$startcodon,
            'width|w=s'              => \$printwidth,
            'graph|g'                => \$graph,
            'verbose|v'              => \$verb,
            'help|h'                 => \$help,
            );

if ($help) {
   exec('perldoc', $0);
   die;
}

if (!defined $sequencefile) {
    die($USAGE . "\nPlease specify an input filename.\n");
}

sub longestORF {
   my $best=0;
   my ($bests,$beste,$beststrand)=(-1,-1,0);
   my $bestorf="";

   my $relaxed=$_[1];

   my $dna=Bio::Seq->new(-seq => $_[0]);
   my %strand=('+'=>$dna->seq);

   foreach my $direction (keys %strand) {
      my @starts=();
      my @ends=();
      if ($relaxed) {
         for (my $frame=0;$frame<3;$frame++) {
            unless ($strand{$direction}=~m/^.{$frame}(taa|tga|tag)/i) {
               push @starts,$frame+1;
            }
         }
      }
      while ($strand{$direction}=~m/(atg)/gi) {
         push @starts,pos($strand{$direction})-2;
      }
      #add starts as found in riboseq
      
      while ($strand{$direction}=~m/(taa|tga|tag)/gi) {
         push @ends,pos($strand{$direction})-2;
      }
      push @ends,($dna->length-2,$dna->length-1,$dna->length);

      for my $s (@starts) {
         for my $e (@ends) {
            if ($e%3==$s%3 and $e>$s) {
               if ($e-$s>$best) {
                  $best=$e-$s;
                  ($bests,$beste,$beststrand)=($s,$e,$direction);

                  $bestorf=Bio::Seq->new(-seq=>$strand{$direction})->subseq($s,$e);
		  #$bestorf.="TAA";
               }
               last
            } else {
               next
            }
         }
      }
   }
   return ($best,$bests,$beste,$beststrand,$bestorf);
}


my $seqio = new Bio::SeqIO('-format' => $sequenceformat,
                           '-file'   => $sequencefile );

my ($length,$start,$end,$direction,$sequence);
my $count=0;
my @lengths;
my $totallength=0;
my %atgHash;

#print $startcodon,"\n";
if($startcodon){open(IN,$startcodon); while(my $line =<IN>){chomp $line; my @tmp = split(/\t+/,$line);  $atgHash{$tmp[6]}{$tmp[7]}=1; } close(IN);}
while (my $dna = $seqio->next_seq) {
  $count++;
  if($atgHash{$dna->display_id()}){foreach my $key (keys %{$atgHash{$dna->display_id()}}){
    next if ($key<=0);
   ($length,$start,$end,$direction,$sequence)=longestORF($dna->subseq($key,$dna->length()),$notstrict);
   if ($verb && length($sequence)>=18) {
      print ">",$dna->display_id," ",$dna->desc,": "; 
       my $s = new Bio::Seq(-seq => $sequence);
      print "$length, ",($start+$key-1),", ",($end+$key-1)," ($direction / ATG: $key),",$dna->length,"\n",$s->seq(),"\n\n"; 
  }
				 }
   
   }
   #$totallength+=$length;
   #$lengths[$length/3]++;
}

#print "Average ORF length: ", $totallength/$count,"\n\n";

#print "Length distribution is:\n";

# oha .. add start codons


__END__

=head1 NAME

longorf - perl script to find the longest ORF of a sequence

=head1 SYNOPSIS

% longorf [-h] [-n] [-v] [-g] [-w printwidth] [-f seqformat] -i seqfile

=head1 DESCRIPTION

This script will examine a set of nucleotide sequences and determine
the longest ORF in each sequence. ORFs may start at the canonical ATG
or at the beginning of the sequence if the notstrict option is chosen.
The script will output a list of the longest ORF lengths, starts, ends
and strands with the ORF and amino acid sequence if the verbose option
is chosen. A histogram of the longest ORFs in the input set may be
printed by choosing the graph option.

=head1 FEEDBACK

This script is not supported by anyone, but requests can be made to the
author.

=head1 AUTHOR - Dan Kortschak <kortschak@rsbs.anu.spanner.edu.au>

=cut
