#!usr/perl/bin
use strict;

my $infile = $ARGV[0]; # Input file Either Ibuprofen or Tamoxifen

print "Please type 0 for default cut-offs and 1 for standard cutoff inputs\n" ;
my $Eopt = <stdin>; # 0 or 1: O is for Default options and 1 is for additional options. 

if ($Eopt == 0)
{
	system ("perl XenoD1.pl $infile Hybrid Upsampling");
	system ("perl XenoD2.pl $infile Hybrid Upsampling");
	system ("perl XenoD3.pl $infile Hybrid");
}
if ($Eopt == 1)
{
	print "Please enter prediction probability threshold at first level: Range 0.4 to 0.8 \n" ;
	chomp (my $FirstP = <stdin>);
	print "Please enter prediction probability threshold at Second level: Range 0.3 to 0.7 \n" ;
	chomp (my $SecondP = <stdin>);
	print "Please enter Tanimoto index valuue for similarity search at Third level: Range 0.3 to 1 \n" ;
	chomp (my $TMI = <stdin>);
	print "Please enter Identiy cut-off to pick the best hit at Third level: Range 40 to 100 \n" ;
	chomp (my $Ide = <stdin>);
	print "Please enter Query coverage cut-off to pick the best hit at Third level: Range 80 to 100 \n" ;
	chomp (my $Qcov = <stdin>);
	print "Please enter E-value cut-off to pick the best hit at Third level: Min: 15 and Max: 0.0 \n" ;
	chomp (my $Eval = <stdin>);
	

	print "$TMI\t$Ide\t$Qcov\t$Eval";
	system ("perl XenoD1.pl $infile Hybrid Upsampling $FirstP");
	system ("perl XenoD2.pl $infile Hybrid Upsampling $SecondP");
	system ("perl XenoD3.pl $infile Hybrid $TMI $Ide $Qcov $Eval");
}
