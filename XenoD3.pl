#! /usr/bin/perl
use strict;
use Data::Dumper; 
#---------------------------------------------------------------------------------------------------------------------------------------------------------#
# This programme is for calculating similarity coefficient and on the basis of similar molecule hit, extract EC, Microbes and Their respective enzymes
# Made by: Ashok K. Sharma, Date 27nov-5Dec 2015
#---------------------------------------------------------------------------------------------------------------------------------------------------------#
my $Query = $ARGV[0]; # Query sdf file #
my $FP = $ARGV[1]; # Fingerprint option for the simialarity calculation#
my $TI1 = $ARGV[2]; # Tanimoto Index Value
my $IDE1 = $ARGV[3]; # Identity cut-off
my $QC1 = $ARGV[4]; # Qcoverage
my $EV1 = $ARGV[5]; # Evalue

my $TI; my $IDE; my $QC; my $EV;

if ($TI1 == "blank") { $TI = 0;} else { $TI = $TI1; }
if ($IDE1 == "blank") { $IDE = 40 ;} else { $IDE = $IDE1; }
if ($QC1 == "blank") { $QC = 80 ;} else { $QC = $QC1; }
if ($EV1 == "blank") { $EV = 15 ;} else { $EV = $EV1; }

#print "$TI\n$IDE\n$QC\n$EV\n"; <stdin>;

#------------------To store two digit EC information obtained after second model ---------------------#
my $Eclass = `ls *response2 | head -1 | cut -f2`;
open (IF, "$Eclass") || die "cant open predicted class at second level:\n";
chomp (my $line = <IF>);
my @arr = split ("\t", $line); my $ec;
if ($arr[1] =~ /EC/){ $ec = $arr[1]; }
else{ $ec = $arr[2]; }
$ec=~ tr/A-Z/a-z/;

#--------------- To store blast output criteria of each hit in hash ----------------------------------#
my $blastfile = "Databases_Sample/Total_HMP_Genomes.tophit.besthit";
open (BF, "$blastfile") || die "can't open complete BLAST best hit best file:\n";

my $hash = {}; 
while (chomp (my $hit = <BF>))
{
# ACDG_00003_23S_rRNA_(uracil-5-)-methyltransferase_RumA_[Acidaminococcus_sp._D21]	sp|Q97LN4|Y523_CLOAB	46.85	429	225	2	34	460	33	460	1e-151	 441	460
# ACDG_00006_acetyltransferase,_GNAT_family_[Acidaminococcus_sp._D21]	sp|Q18B70|ATSE_PEPD6	52.29	109	52	0	7	115	51	159	6e-33	 114	115
	
	my @arr = split ("\t", $hit); my @arr1 = split(/_\[/, $arr[0]);
	my $anot = $arr1[0]; my $ide = $arr[2]; my $Qlen = $arr[12]; my $Q_Alen = $arr[7]-$arr[6]; my $Qcov = ($Q_Alen/$Qlen)*100; my $rounded = sprintf "%.2f", $Qcov;
	my $evalue = $arr[10];
	$hash->{$anot}="$ide\t$rounded\t$evalue";
}

#-------------------Open a file that contain molecules list-------------------#
my $file1 = "Databases_Sample/$ec.database/order_$ec";
open(IF1,"$file1")|| die "Not open:\n";

#1.---------------------------------------------------------------------- Hybrid Approach ---------------------------------------------------------------------# 
if ($FP =~ /Hybrid/)
{
	while(chomp($line=<IF1>))
	{
		system ("babel Databases_Sample/$ec.database/$line.mol.smi $Query -ofpt -xfMACCS > $line.MACCSFP.out && sed -i '1d' $line.MACCSFP.out");
		open (IF2, "$line.MACCSFP.out");
		open (OF2, ">>$line.MACCSFP.out1");
		chomp (my $in = <IF2>); print OF2 "$line\t$in\n";

		system ("babel Databases_Sample/$ec.database/$line.mol.smi $Query -ofpt -xfFP2 > $line.FP2.out && sed -i '1d' $line.FP2.out");
		open (IF3, "$line.FP2.out");
		open (OF3, ">>$line.FP2.out2");
		chomp (my $in = <IF3>); print OF3 "$line\t$in\n";

		system ("babel Databases_Sample/$ec.database/$line.mol.smi $Query -ofpt -xfFP4 > $line.FP4.out && sed -i '1d' $line.FP4.out");
		open (IF4, "$line.FP4.out");
		open (OF4, ">>$line.FP4.out3");
		chomp (my $in = <IF4>); print OF4 "$line\t$in\n";	    
	}
	system ("cat *out1 > MACCSFP.out.Final && sed -i 's/ = /\t/g' MACCSFP.out.Final && sort -k3 -r MACCSFP.out.Final > MACCSFP.out.Final.Sort.tmp");
	system ("cat *out2 > FP2.out.Final && sed -i 's/ = /\t/g' FP2.out.Final && sort -k3 -r FP2.out.Final > FP2.out.Final.Sort.tmp");
	system ("cat *out3 > FP4.out.Final && sed -i 's/ = /\t/g' FP4.out.Final && sort -k3 -r FP4.out.Final > FP4.out.Final.Sort.tmp");

	my $maccsfp = "MACCSFP.out.Final.Sort.tmp";
	open (IF5, "$maccsfp") || die "cant open file predicted class after similarity search at third level using MACCSFP:\n";
	chomp (my $c1 = <IF5>);
	my @arr1 = split ("\t", $c1); my $C1 = $arr1[0];
	
	my $fp2 = "FP2.out.Final.Sort.tmp";
	open (IF6, "$fp2") || die "cant open file predicted class after similarity search at third level using FP2:\n";
	chomp (my $c2 = <IF6>);
	my @arr2 = split ("\t", $c2); my $C2 = $arr2[0];
	
	my $fp4 = "FP4.out.Final.Sort.tmp";
	open (IF7, "$fp4") || die "cant open file predicted class after similarity search at third level using FP4:\n";
	chomp (my $c3 = <IF7>);
	my @arr3 = split ("\t", $c3); my $C3 = $arr3[0];

	open (OF7, ">All.Hybrid.Final") || die "cant open file ";
	if ($C1 =~ /$C2/ && $C2 =~ /$C3/ && $C3 =~ /$C1/) { print OF7 "$c1\n$c2\n$c3\n"; }
	elsif ($C1 =~ /$C2/ && $C2 !=~ /$C3/) { print OF7 "$c1\n$c2\n"; }
	elsif ($C2 =~ /$C3/ && $C2 !=~ /$C1/) { print OF7 "$c2\n$c3\n"; }
	elsif ($C3 =~ /$C1/ && $C1 !=~ /$C2/) { print OF7 "$c1\n$c3\n"; }
	else 
	{
		my @arr4 = split ("\t", $c1); my $value1 = $arr4[2];  my @arr5 = split ("\t", $c2); my $value2 = $arr5[2]; my @arr6 = split ("\t", $c3); my $value3 = $arr6[2];
		my $max = 0; 
		if ($value1 >= $value2)
		{
			 $max == $value1;
			 if ($max >= $value3) { print OF7 "$c1\tMACCSFP\n"; }
			 else { print OF7 "$c3\tFP4\n";}
		}
		else 
		{
			$max == $value2;
			if ($max >= $value3) { print OF7 "$c2\tFP2\n"; }
			else { print OF7 "$c3\tFP4\n"; }
		}     
	}
system ("sed -i 's/ = /\t/g' All.Hybrid.Final && sort -k3 -r All.Hybrid.Final > All.Hybrid.Final.Sort");
system ("rm *out *out1 *out2 *out3 *Final *Sort.tmp");
close IF1; close IF2; close OF2; close IF3; close OF3; close IF4; close OF4; close IF5; close IF6, close IF7; close OF7;
}

#2.----------------------------------------------------------------------------MACCSFP---------------------------------------------------------------------------------#
if ($FP =~ /MACCSFP/)
{
	while(chomp($line=<IF1>))
	{
		system ("babel Databases_Sample/$ec.database/$line.mol.smi $Query -ofpt -xfMACCS > $line.MACCSFP.out && sed -i '1d' $line.MACCSFP.out");
		open (IF2, "$line.MACCSFP.out");
		open (OF2, ">>$line.MACCSFP.out1");
		chomp (my $in = <IF2>); print OF2 "$line\t$in\n";
	}
	system ("cat *out1 > MACCSFP.out.Final && sed -i 's/ = /\t/g' MACCSFP.out.Final && sort -k3 -r MACCSFP.out.Final > MACCSFP.out.Final.Sort && rm *out *out1 *Final");
close IF1; close IF2; close OF2;
}
#3.-----------------------------------------------------------------------------FP2-------------------------------------------------------------------------------------#
if ($FP =~ /FP2/)
{
	while(chomp($line=<IF1>))
	{
		system ("babel Databases_Sample/$ec.database/$line.mol.smi $Query -ofpt -xfFP2 > $line.FP2.out && sed -i '1d' $line.FP2.out");
		open (IF2, "$line.FP2.out");
		open (OF2, ">>$line.FP2.out1");
		chomp (my $in = <IF2>); print OF2 "$line\t$in\n";
	}
	system ("cat *out1 > FP2.out.Final && sed -i 's/ = /\t/g' FP2.out.Final && sort -k3 -r FP2.out.Final > FP2.out.Final.Sort && rm *out *out1 *Final");
close IF1; close IF2; close OF2;
}
#4.-----------------------------------------------------------------------------FP4-------------------------------------------------------------------------------------#
if ($FP =~ /FP4/)
{
	while(chomp($line=<IF1>))
	{
		system ("babel Databases_Sample/$ec.database/$line.mol.smi $Query -ofpt -xfFP4 > $line.FP4.out && sed -i '1d' $line.FP4.out");
		open (IF2, "$line.FP4.out");
		open (OF2, ">>$line.FP4.out1");
		chomp (my $in = <IF2>); print OF2 "$line\t$in\n";
	}
	system ("cat *out1 > FP4.out.Final && sed -i 's/ = /\t/g' FP4.out.Final && sort -k3 -r FP4.out.Final > FP4.out.Final.Sort && rm *out *out1 *Final");
close IF1; close IF2; close OF2;
}


#---Part 2: To use Final Sort output from previous methods: Hybrid, MACCSFP, FP2 and FP4. File contains Cnumber and their tanimoto coefficient values in third column---#

my $filename = `ls *Final.Sort | head -1 | cut -f2`;
open (IF3, "$filename") || die "cant open file predicted class after similarity search at third level:\n";

chomp (my $cnumber = <IF3>);
my @arr1 = split ("\t", $cnumber); my $C = $arr1[0]; my $TIcutoff = $arr1[2]; #print "$TI\t$TIcutoff"; <stdin>;

if ($TIcutoff < $TI)
{
	 open (RNF, ">Error"); print RNF "No Significant hit found. Please Use Less Tanimoto Index Cut-off."; 
	system ("rm $ec.proteins.Final.Result.Genus $ec.proteins.Final.Result $ec.proteins.Final.Result1 Error1"); exit;
}
else
{
	system ("rm Error");
	system ("grep -w $C \"Databases_Sample/$ec.Detail.Final\" > $ec.temporary");
	open (IF4, "$ec.temporary") || die "cant open temporary file contain C id:\n";
	while (chomp (my $enumber = <IF4>))
	{
		my @arr2 = split("\t", $enumber); my $E = $arr2[0]; 
		#my @Esplit = split(/\./, $E); 
		system ("grep -F -w \'$E\' \"Databases_Sample/492_GenomesAnot_EC.out.New\" >> $ec.proteins");
		#system ("grep -w \"$Esplit[0].$Esplit[1].$Esplit[2].$Esplit[3]\" \"Databases_Sample/492_GenomesAnot_EC.out.New\" >> $ec.proteins");
		#system ("grep -w $E \"Databases_Sample/492_GenomesAnot_EC.out.New\" >> $ec.proteins");
	}
		my $hash1 = {}; # New hash to store the microbes and their proteins info in EC. 
		open (IF5, "$ec.proteins") || die "Error: opening $ec.proteins:\n";
		open (OF5, ">$ec.proteins.Final_1") || die "Error: opening :\n";
		
		while (chomp (my $line5 = <IF5>))
		{
		# 1.14.14.1	HMPREF1207_00390_hypothetical_protein_[Paenibacillus_sp._HGH0039];HMPREF9412_3285_bifunctional_P-450/NADPH-P450_reductase_[Paenibacillus_sp._HGF5];HMPREF1211_00829_hypothetical_protein_[Streptomyces_sp._HGB0020];HMPREF0989_04160_hypothetical_protein_[Ralstonia_sp._5_2_56FAA];HMPREF1004_00113_FAD_binding_domain_protein_[Ralstonia_sp._5_7_47FAA];HMPREF9413_3023_FAD_binding_domain_protein_[Paenibacillus_sp._HGF7];HMPREF1014_01714_hypothetical_protein_[Bacillus_sp._7_6_55CFAA_CT2];
		# 2.3.1.5	HMPREF1211_07660_hypothetical_protein_[Streptomyces_sp._HGB0020];

			my @arr5 = split("\t", $line5); my @arr6 = split(/;/, $arr5[1]); my $ec_new = $arr5[0];
			foreach my $detail (@arr6)
			{
				my @arr7 = split(/_\[/, $detail); my $microbe = $arr7[1]; $microbe =~ s/\]//g; my $enzyme = $arr7[0];
				push(@{$hash1->{$ec_new}->{$microbe}->{$enzyme}});
			}
		}

		foreach my $val ( keys %{$hash1})
		{
				foreach my $val1(keys %{$hash1->{$val}})
				{	
					foreach my $val2(keys %{$hash1->{$val}->{$val1}})
					{
						#--------- For printing of Ide, Qcov and Evalues for each hit from the very first BLAST Output hash---------#
						if (exists $hash->{$val2})
						{
						    print OF5 "$val\t$val1\t$val2\t".$hash->{$val2}."\n"; 
						    #print OF5 "\t\t$val2\t".$hash->{$val2}."\n";
						}
					}
				}
		}
	close IF4; close IF5; close OF5;

		open (IF6, "$ec.proteins.Final_1") || die "Error: opening $ec.proteins:\n";
		open (OF6, ">$ec.proteins.Final.Result1") || die "Error: opening :\n";
		print OF6 "EC\tMicrobes\tProteins\tIdentity\tQ-coverage\tE-value\n";
		while (chomp (my $line6 = <IF6>))
		{
# 2.3.1.-	Prevotella_salivae_DSM_15606	HMPREF9420_2097_bacterial_transferase_hexapeptide_repeat_protein	45.31	92.68	8e-43
# 2.3.1.-	Prevotella_salivae_DSM_15606	HMPREF9420_0478_UDP-3-O-[3-hydroxymyristoyl]_glucosamine_N-acyltransferase	63.77	99.87	8e-41
			    #print $line6; <stdin>;
			    my @arr = split ("\t", $line6); #print $arr[4]; <stdin>;
			    my $x = $arr[3]; my $y = $arr[4]; my $evalue = $arr[5]; my @eval = split("-", $evalue); my $z = $eval[1];

			    if ($EV == 0)
			    {
					if (($x >= $IDE) && ($y >= $QC) && ($evalue == 0.0)) 
					{
						  print OF6 "$line6\n";
					}
			    }
			    else
			    {
					if (($x >= $IDE) && ($y >= $QC) && ($z >= $EV)) 
					{
						  print OF6 "$line6\n";
					}
					if (($x >= $IDE) && ($y >= $QC) && ($evalue == 0.0)) 
					{
						  print OF6 "$line6\n";
					}
			    }
		}

	system ("sed 's/_/ /g' $ec.proteins.Final.Result1 > $ec.proteins.Final.Result1.temporary"); 
	system ("sed 's/\t/;/g' $ec.proteins.Final.Result1.temporary > $ec.proteins.Final.Result");

	system ("rm *temporary $ec.proteins $ec.proteins.Final_1"); 

	#------------------------------------------ To Extract Genus Information --------------------------------------------------------------------------------------
	my $hashG = {};
	open(GI, "Databases_Sample/Species_Taxonnomy_492"); # always coply this file to the directory where we will run this program
	while (chomp(my $lineG = <GI>)) 
	{
		my @arrG = split('\t', $lineG);
		$hashG->{$arrG[0]}="$arrG[1]";
	}
	open(RF, "$ec.proteins.Final.Result1");
	open(ROF, ">>genus.out");

	while (chomp(my $lineR = <RF>))
	{
	# 2.3.1.-	Prevotella_salivae_DSM_15606	HMPREF9420_2097_bacterial_transferase_hexapeptide_repeat_protein	Identity:45.31	QueryCoverage:92.68	Evalue:8e-43
	# 2.3.1.-	Prevotella_salivae_DSM_15606	HMPREF9420_0478_UDP-3-O-[3-hydroxymyristoyl]_glucosamine_N-acyltransferase	Identity:63.77	QueryCoverage:99.42	Evalue:1e-162	
		
		my @arrR = split ("\t", $lineR);
		my $sp = $arrR[1]; 
		if (exists ($hashG->{$sp}))
		{
			print ROF $hashG->{$sp}."\n";
		}
	}
	system("sort genus.out | uniq -c > final.genus && sort -nr final.genus > final.genus1 && rm final.genus genus.out");
	system ("sed -i 's/^  *//g' final.genus1 && cut -d' ' -f2 final.genus1 > $ec.proteins.Final.Result.Genus && rm final.genus1");
	
	my $Rlines = `wc -l $ec.proteins.Final.Result1 | cut -d" " -f1`; #print $Rlines; <stdin>;
	if ($Rlines == 1)
	{
		open (EF, ">Error1"); print EF "No Significant hit found. Try Using Less Cut-off for Identity, Q-coverage and E-value."; 
		system ("rm $ec.proteins.Final.Result.Genus $ec.proteins.Final.Result $ec.proteins.Final.Result1 Error"); exit;
	}
	else { system ("rm Error1"); }
	

	#`tr '\\n' ',' <$ec.proteins.Final.Result.Genus | sed "s/,/, /g" >$ec.proteins.Final.Result.Genus1; mv $ec.proteins.Final.Result.Genus1 $ec.proteins.Final.Result.Genus`;
system ("rm *Sort *log *Result *Pred* && mv $ec.proteins.Final.Result1 $Query.proteins.Final.Result && mv $ec.proteins.Final.Result.Genus $Query.Genus.Final.Result");
}
#----------------------------------------------------------------------------------End ------------------------------------------------------------#
