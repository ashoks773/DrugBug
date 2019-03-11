#!/usr/bin/perl

use strict;

#----------------------------------------XXXXXXXXXXXXXX------------------------------------------------
# This is the pipeline for the prediction of EC Class at second level for a query molecule#
# Date: 23,24-Nov-2015
# Written by: Ashok K Sharma
#----------------------------------------XXXXXXXXXXXXXX------------------------------------------------

my $inputfile = $ARGV[0]; # input molecule file (.mol) or Cid
my $model = $ARGV[1]; # input will be either a fingerprint/descripter or Hybrid
my $sample = $ARGV[2]; # input will be either normal or upsampling
my $PP = $ARGV[3]; # Additional option Of prediction probabilities
my $Thr; my $threshold; 
if ($PP == "blank") { $Thr = 0.5;}
else { $Thr = $PP; }
#print $Thr; <stdin>;

my $classname = `ls *response | head -1 | cut -f2`;
open (IF, "$classname") || die "cant open predicted class at second level:\n";
chomp (my $line = <IF>);
my @arr = split ("\t", $line);
my $ec;

if ($arr[1] =~ /EC/){ $ec = $arr[1]; }
else{ $ec = $arr[2]; }
$ec=~ tr/A-Z/a-z/;


#*-----------------------------------------------------------------------------------------------EC1-----------------------------------------------------------------------------------------------------------------------------#

if ($ec =~ /ec1/)
{
	                                                    #---------------------------- Hybrid--------------------------------------#
	  if ($model =~ /Hybrid/)
	  {
		    system ("cut -d\",\" -f2-10209 Hybrid_fingerprints.out > Hybrid_fingerprints.Out");
		    system ("cut -d\",\" -f17,18,87,92,122,128,159,193,200,212,221,239,243,254,258,266,288,315,356,359,394,412,464,477,495,508,512,514,525,529,549,562,598,633,658,729,744,801,802,832,837,853,884,929,937,949,974,997,1002,1005,1007 Hybrid_fingerprints.Out > 1");
		    system ("cut -d\",\" -f1032,1035,1042,1045,1048,1057,1058,1060 Hybrid_fingerprints.Out > 2");
		    system ("cut -d\",\" -f1133,1195,1196,1253,1259,1269,1331,1345,1367,1389,1409,1410,1412,1431,1436,1465,1801,1825,1826,1828,1834,1855,1937,1989,2023,2085,2108 Hybrid_fingerprints.Out > 3");
		    system ("cut -d\",\" -f2140,2153,2161,2170,2176,2177,2180,2184,2211,2217,2225,2233,2237,2250,2252,2253,2263,2264,2267,2269,2273,2282,2285  Hybrid_fingerprints.Out > 4");
		    system ("cut -d\",\" -f2307,2309,2314,2315,2475,2480,2549,2617,2627,2659,2674,2699,2726,2730,2734,2745,2756,2822,2908,2953,2956,2958,2980,3112 Hybrid_fingerprints.Out > 5");
		    system ("cut -d\",\" -f3190,3192,3221,3222,3223,3249,3262,3317,3378,3449,3455,3456,3470 Hybrid_fingerprints.Out > 6");
		    system ("cut -d\",\" -f3489,3497,3670,4158,4887,5047,6028,6045,6601,6641,6749,6863,6869,7135,7149,7180,7210,7225,7270,7287,7301,7302,7363,7438,7496,7603,7674,7718,7776,8176,8292,8298,8302,8307 Hybrid_fingerprints.Out > 7");
		    system ("cut -d\",\" -f8343,8433,8443,8677,8679,8744,8911,8989 Hybrid_fingerprints.Out > 8");
		    system ("cut -d\",\" -f9124,9126,9137,9168,9169,9170,9209,9264,9325,9395,9396,9402,9403,9416,9417,9423 Hybrid_fingerprints.Out > 9");
		    system ("cut -d\",\" -f9429,9430,9431,9508,9520,9530,9585,9608,9667,9764,9831,9842,9910,9998,10076 Hybrid_fingerprints.Out > 10");
		    system ("paste -d\",\" 1 2 3 4 5 6 7 8 9 10 > Hybrid_fingerprints.Out.Final");
		    system ("rm 1 2 3 4 5 6 7 8 9 10 Hybrid_fingerprints.Out");

		    if ($sample =~ /Upsampling/) { system ("sh ec1RscriptHybridUp.sh"); }
		    if ($sample =~ /Normal/) { system ("sh ec1RscriptHybrid.sh"); }
		    system ("sed -i '1d' Hybrid_Pred2");

		    my $Rfpred = "Hybrid_Pred2"; # prediction file with probabilities obtained from the random forest commnad.
		    $threshold = $Thr;        # $ARGV[1]; # cutoff thresold for the prediction 
		    open (PR, $Rfpred) || die "Cant open the file:\n";
		    open (PRR, ">$Rfpred.response2") || die "Cant open the file:\n";

		    while (chomp (my $line = <PR>))
		    {
	            #	 "EC1.1"	"EC1.10"	"EC1.11"	"EC1.13"	"EC1.14"	"EC1.15"	"EC1.16"	"EC1.17"	"EC1.18"	"EC1.2"	"EC1.20"	"EC1.3"	"EC1.4"	"EC1.5"	"EC1.6"	"EC1.7"	"EC1.8"	"EC1.97"
	            #"1"  0.14	0	0	0.064	0.172	0	0.008	0.012	0	0.122	0	0.084	0.036	0.328	0	0.026	0.008	0
			    
			    my @arr = split ("\t", $line);
			    
			    my %hash; $hash{$arr[1]}="EC1.1"; $hash{$arr[2]}="EC1.10"; $hash{$arr[3]}="EC1.11"; $hash{$arr[4]}="EC1.13"; $hash{$arr[5]}="EC1.14"; $hash{$arr[6]}="EC1.15"; $hash{$arr[7]}="EC1.16";
			    $hash{$arr[8]}="EC1.17"; $hash{$arr[9]}="EC1.18"; $hash{$arr[10]}="EC1.2"; $hash{$arr[11]}="EC1.20"; $hash{$arr[12]}="EC1.3"; $hash{$arr[13]}="EC1.4"; $hash{$arr[14]}="EC1.5"; $hash{$arr[15]}="EC1.6";
			    $hash{$arr[16]}="EC1.7"; $hash{$arr[17]}="EC1.8"; $hash{$arr[18]}="EC1.97";
			    
			    my $x = 0; foreach my $y (@arr){ $x = $y if $y > $x; }
			    
			    if ($arr[1] > $threshold) { print PRR $arr[0]."\tEC1.1\n"; }
			    elsif ($arr[2] > $threshold) { print PRR $arr[0]."\tEC1.10\n"; }
			    elsif ($arr[3] > $threshold) { print PRR $arr[0]."\tEC1.11\n"; }
			    elsif ($arr[4] > $threshold) { print PRR $arr[0]."\tEC1.13\n"; }
			    elsif ($arr[5] > $threshold) { print PRR $arr[0]."\tEC1.14\n"; }
			    elsif ($arr[6] > $threshold) { print PRR $arr[0]."\tEC1.15\n"; }
			    elsif ($arr[7] > $threshold) { print PRR $arr[0]."\tEC1.16\n"; }
			    elsif ($arr[8] > $threshold) { print PRR $arr[0]."\tEC1.17\n"; }
			    elsif ($arr[9] > $threshold) { print PRR $arr[0]."\tEC1.18\n"; }
			    elsif ($arr[10] > $threshold) { print PRR $arr[0]."\tEC1.2\n"; }
			    elsif ($arr[11] > $threshold) { print PRR $arr[0]."\tEC1.20\n"; }
			    elsif ($arr[12] > $threshold) { print PRR $arr[0]."\tEC1.3\n"; }
			    elsif ($arr[13] > $threshold) { print PRR $arr[0]."\tEC1.4\n"; }
			    elsif ($arr[14] > $threshold) { print PRR $arr[0]."\tEC1.5\n"; }
			    elsif ($arr[15] > $threshold) { print PRR $arr[0]."\tEC1.6\n"; }
			    elsif ($arr[16] > $threshold) { print PRR $arr[0]."\tEC1.7\n"; }
			    elsif ($arr[17] > $threshold) { print PRR $arr[0]."\tEC1.8\n"; }
			    elsif ($arr[18] > $threshold) { print PRR $arr[0]."\tEC1.97\n";}
			    else { print PRR $arr[0]."\tNP (Probability below thresold)\t$hash{$x}\n"; }# We have to decide wheteher this class should included or not

		    }
	  close PR; close PRR;
	  }
	                                                                  #---------------------------------------- KRFP------------------------------------------------#
	  if ($model =~ /KRFP/)
	  {
		    system ("cut -d\",\" -f3483-8342 Hybrid_fingerprints.out > KRFP.Out.Final");
		    #system ("rm Fingerprinter.out");

		    if ($sample =~ /Upsampling/) { system ("sh ec1RscriptKRFPUp.sh"); }
		    if ($sample =~ /Normal/) { system ("sh ec1RscriptKRFP.sh"); }
		    system ("sed -i '1d' KRFP_Pred2");

		    my $Rfpred = "KRFP_Pred2"; # prediction file with probabilities obtained from the random forest commnad.
		    $threshold = $Thr;        # $ARGV[1]; # cutoff thresold for the prediction 
		    open (PR, $Rfpred) || die "Cant open the file:\n";
		    open (PRR, ">$Rfpred.response2") || die "Cant open the file:\n";

		    while (chomp (my $line = <PR>))
		    {
		    #	"EC1.1"	"EC1.10"	"EC1.11"	"EC1.13"	"EC1.14"	"EC1.15"	"EC1.16"	"EC1.17"	"EC1.18"	"EC1.2"	"EC1.20"	"EC1.3"	"EC1.4"	"EC1.5"	"EC1.6"	"EC1.7"	"EC1.8"	"EC1.97"
		#"1"	0.14	0	0	0.064	0.172	0	0.008	0.012	0	0.122	0	0.084	0.036	0.328	0	0.026	0.008	0
			    
			    my @arr = split ("\t", $line);
			    
			    my %hash; $hash{$arr[1]}="EC1.1"; $hash{$arr[2]}="EC1.10"; $hash{$arr[3]}="EC1.11"; $hash{$arr[4]}="EC1.13"; $hash{$arr[5]}="EC1.14"; $hash{$arr[6]}="EC1.15"; $hash{$arr[7]}="EC1.16";
			    $hash{$arr[8]}="EC1.17"; $hash{$arr[9]}="EC1.18"; $hash{$arr[10]}="EC1.2"; $hash{$arr[11]}="EC1.20"; $hash{$arr[12]}="EC1.3"; $hash{$arr[13]}="EC1.4"; $hash{$arr[14]}="EC1.5"; $hash{$arr[15]}="EC1.6";
			    $hash{$arr[16]}="EC1.7"; $hash{$arr[17]}="EC1.8"; $hash{$arr[18]}="EC1.97";
			    
			    my $x = 0; foreach my $y (@arr){ $x = $y if $y > $x; }
			    
			    if ($arr[1] > $threshold) { print PRR $arr[0]."\tEC1.1\n"; }
			    elsif ($arr[2] > $threshold) { print PRR $arr[0]."\tEC1.10\n"; }
			    elsif ($arr[3] > $threshold) { print PRR $arr[0]."\tEC1.11\n"; }
			    elsif ($arr[4] > $threshold) { print PRR $arr[0]."\tEC1.13\n"; }
			    elsif ($arr[5] > $threshold) { print PRR $arr[0]."\tEC1.14\n"; }
			    elsif ($arr[6] > $threshold) { print PRR $arr[0]."\tEC1.15\n"; }
			    elsif ($arr[7] > $threshold) { print PRR $arr[0]."\tEC1.16\n"; }
			    elsif ($arr[8] > $threshold) { print PRR $arr[0]."\tEC1.17\n"; }
			    elsif ($arr[9] > $threshold) { print PRR $arr[0]."\tEC1.18\n"; }
			    elsif ($arr[10] > $threshold) { print PRR $arr[0]."\tEC1.2\n"; }
			    elsif ($arr[11] > $threshold) { print PRR $arr[0]."\tEC1.20\n"; }
			    elsif ($arr[12] > $threshold) { print PRR $arr[0]."\tEC1.3\n"; }
			    elsif ($arr[13] > $threshold) { print PRR $arr[0]."\tEC1.4\n"; }
			    elsif ($arr[14] > $threshold) { print PRR $arr[0]."\tEC1.5\n"; }
			    elsif ($arr[15] > $threshold) { print PRR $arr[0]."\tEC1.6\n"; }
			    elsif ($arr[16] > $threshold) { print PRR $arr[0]."\tEC1.7\n"; }
			    elsif ($arr[17] > $threshold) { print PRR $arr[0]."\tEC1.8\n"; }
			    elsif ($arr[18] > $threshold) { print PRR $arr[0]."\tEC1.97\n";}
			    else { print PRR $arr[0]."\tNP (Probability below thresold)\t$hash{$x}\n"; }# We have to decide wheteher this class should included or not
		    }
	  close PR; close PRR;
	}
	                                                                         #---------------------------------MACCSFP------------------------------------------------#
	  if ($model =~ /MACCSFP/)
	  {
		    system ("cut -d\",\" -f2129-2294 Hybrid_fingerprints.out  > MACCSFP.Out.Final");

		    if ($sample =~ /Upsampling/) { system ("sh ec1RscriptMACCSFPUp.sh"); }
		    if ($sample =~ /Normal/) { system ("sh ec1RscriptMACCSFP.sh"); }
		    system ("sed -i '1d' MACCSFP_Pred2");

		    my $Rfpred = "MACCSFP_Pred2"; # prediction file with probabilities obtained from the random forest commnad.
		    $threshold = $Thr;        # $ARGV[1]; # cutoff thresold for the prediction 
		    open (PR, $Rfpred) || die "Cant open the file:\n";
		    open (PRR, ">$Rfpred.response2") || die "Cant open the file:\n";

		    while (chomp (my $line = <PR>))
		    {
		    #	"EC1.1"	"EC1.10"	"EC1.11"	"EC1.13"	"EC1.14"	"EC1.15"	"EC1.16"	"EC1.17"	"EC1.18"	"EC1.2"	"EC1.20"	"EC1.3"	"EC1.4"	"EC1.5"	"EC1.6"	"EC1.7"	"EC1.8"	"EC1.97"
		#"1"	0.14	0	0	0.064	0.172	0	0.008	0.012	0	0.122	0	0.084	0.036	0.328	0	0.026	0.008	0
			    
			    my @arr = split ("\t", $line);
			    
			    my %hash; $hash{$arr[1]}="EC1.1"; $hash{$arr[2]}="EC1.10"; $hash{$arr[3]}="EC1.11"; $hash{$arr[4]}="EC1.13"; $hash{$arr[5]}="EC1.14"; $hash{$arr[6]}="EC1.15"; $hash{$arr[7]}="EC1.16";
			    $hash{$arr[8]}="EC1.17"; $hash{$arr[9]}="EC1.18"; $hash{$arr[10]}="EC1.2"; $hash{$arr[11]}="EC1.20"; $hash{$arr[12]}="EC1.3"; $hash{$arr[13]}="EC1.4"; $hash{$arr[14]}="EC1.5"; $hash{$arr[15]}="EC1.6";
			    $hash{$arr[16]}="EC1.7"; $hash{$arr[17]}="EC1.8"; $hash{$arr[18]}="EC1.97";
			    
			    my $x = 0; foreach my $y (@arr){ $x = $y if $y > $x; }
			    
			    if ($arr[1] > $threshold) { print PRR $arr[0]."\tEC1.1\n"; }
			    elsif ($arr[2] > $threshold) { print PRR $arr[0]."\tEC1.10\n"; }
			    elsif ($arr[3] > $threshold) { print PRR $arr[0]."\tEC1.11\n"; }
			    elsif ($arr[4] > $threshold) { print PRR $arr[0]."\tEC1.13\n"; }
			    elsif ($arr[5] > $threshold) { print PRR $arr[0]."\tEC1.14\n"; }
			    elsif ($arr[6] > $threshold) { print PRR $arr[0]."\tEC1.15\n"; }
			    elsif ($arr[7] > $threshold) { print PRR $arr[0]."\tEC1.16\n"; }
			    elsif ($arr[8] > $threshold) { print PRR $arr[0]."\tEC1.17\n"; }
			    elsif ($arr[9] > $threshold) { print PRR $arr[0]."\tEC1.18\n"; }
			    elsif ($arr[10] > $threshold) { print PRR $arr[0]."\tEC1.2\n"; }
			    elsif ($arr[11] > $threshold) { print PRR $arr[0]."\tEC1.20\n"; }
			    elsif ($arr[12] > $threshold) { print PRR $arr[0]."\tEC1.3\n"; }
			    elsif ($arr[13] > $threshold) { print PRR $arr[0]."\tEC1.4\n"; }
			    elsif ($arr[14] > $threshold) { print PRR $arr[0]."\tEC1.5\n"; }
			    elsif ($arr[15] > $threshold) { print PRR $arr[0]."\tEC1.6\n"; }
			    elsif ($arr[16] > $threshold) { print PRR $arr[0]."\tEC1.7\n"; }
			    elsif ($arr[17] > $threshold) { print PRR $arr[0]."\tEC1.8\n"; }
			    elsif ($arr[18] > $threshold) { print PRR $arr[0]."\tEC1.97\n";}
			    else { print PRR $arr[0]."\tNP (Probability below thresold)\t$hash{$x}\n"; }# We have to decide wheteher this class should included or not
		    }
	  close PR; close PRR;
	  }
system ("rm *Out*"); 
}
#*-----------------------------------------------------------------------------------------------EC2-----------------------------------------------------------------------------------------------------------------------------#

if ($ec =~ /ec2/)
{
	                                                     #--------------------------------- Hybrid-----------------------------------#
	  if ($model =~ /Hybrid/)
	  {
		    system ("cut -d\",\" -f2-10209 Hybrid_fingerprints.out > Hybrid_fingerprints.Out");		    
		    system ("cut -d\",\" -f17,18,87,92,122,128,159,193,200,212,221,239,243,254,258,266,288,315,356,359,394,412,464,477,495,508,512,514,525,529,549,562,598,633,658,729,744,801,802,832,837,853,884,929,937,949,974,997,1002,1005,1007 Hybrid_fingerprints.Out > 1");
		    system ("cut -d\",\" -f1032,1035,1042,1045,1048,1057,1058,1060 Hybrid_fingerprints.Out > 2");
		    system ("cut -d\",\" -f1133,1195,1196,1253,1259,1269,1331,1345,1367,1389,1409,1410,1412,1431,1436,1465,1801,1825,1826,1828,1834,1855,1937,1989,2023,2085,2108 Hybrid_fingerprints.Out > 3");
		    system ("cut -d\",\" -f2140,2153,2161,2170,2176,2177,2180,2184,2211,2217,2225,2233,2237,2250,2252,2253,2263,2264,2267,2269,2273,2282,2285  Hybrid_fingerprints.Out > 4");
		    system ("cut -d\",\" -f2307,2309,2314,2315,2475,2480,2549,2617,2627,2659,2674,2699,2726,2730,2734,2745,2756,2822,2908,2953,2956,2958,2980,3112 Hybrid_fingerprints.Out > 5");
		    system ("cut -d\",\" -f3190,3192,3221,3222,3223,3249,3262,3317,3378,3449,3455,3456,3470 Hybrid_fingerprints.Out > 6");
		    system ("cut -d\",\" -f3489,3497,3670,4158,4887,5047,6028,6045,6601,6641,6749,6863,6869,7135,7149,7180,7210,7225,7270,7287,7301,7302,7363,7438,7496,7603,7674,7718,7776,8176,8292,8298,8302,8307 Hybrid_fingerprints.Out > 7");
		    system ("cut -d\",\" -f8343,8433,8443,8677,8679,8744,8911,8989 Hybrid_fingerprints.Out > 8");
		    system ("cut -d\",\" -f9124,9126,9137,9168,9169,9170,9209,9264,9325,9395,9396,9402,9403,9416,9417,9423 Hybrid_fingerprints.Out > 9");
		    system ("cut -d\",\" -f9429,9430,9431,9508,9520,9530,9585,9608,9667,9764,9831,9842,9910,9998,10076 Hybrid_fingerprints.Out > 10");

		    system ("paste -d\",\" 1 2 3 4 5 6 7 8 9 10 > Hybrid_fingerprints.Out.Final");
		    system ("rm 1 2 3 4 5 6 7 8 9 10 Hybrid_fingerprints.Out");

		    if ($sample =~ /Upsampling/) { system ("sh ec2RscriptHybridUp.sh");}
		    if ($sample =~ /Normal/) { system ("sh ec2RscriptHybrid.sh"); }
		    system ("sed -i '1d' Hybrid_Pred2");

		    my $Rfpred = "Hybrid_Pred2"; # prediction file with probabilities obtained from the random forest commnad.
		    $threshold = $Thr;        # $ARGV[1]; # cutoff thresold for the prediction 
		    open (PR, $Rfpred) || die "Cant open the file:\n";
		    open (PRR, ">$Rfpred.response2") || die "Cant open the file:\n";

		    while (chomp (my $line = <PR>))
		    {
	  #	       "EC2.1"	"EC2.10"	"EC2.2"	"EC2.3"	"EC2.4"	"EC2.5"	"EC2.6"	"EC2.7"	"EC2.8"	"EC2.9"
                  #"1"	0.334	0.02	0	0.166	0.208	0.136	0.014	0.094	0.022	0.006
			    
			    my @arr = split ("\t", $line);

			    my %hash; $hash{$arr[1]}="EC2.1"; $hash{$arr[2]}="EC2.10"; $hash{$arr[3]}="EC2.2"; $hash{$arr[4]}="EC2.3"; $hash{$arr[5]}="EC2.4"; $hash{$arr[6]}="EC2.5"; $hash{$arr[7]}="EC2.6";
			    $hash{$arr[8]}="EC2.7"; $hash{$arr[9]}="EC2.8"; $hash{$arr[10]}="EC2.9"; 
			    
			    my $x = 0; foreach my $y (@arr){ $x = $y if $y > $x; }
			    
			    if ($arr[1] > $threshold) { print PRR $arr[0]."\tEC2.1\n"; }
			    elsif ($arr[2] > $threshold) { print PRR $arr[0]."\tEC2.10\n"; }
			    elsif ($arr[3] > $threshold) { print PRR $arr[0]."\tEC2.2\n"; }
			    elsif ($arr[4] > $threshold) { print PRR $arr[0]."\tEC2.3\n"; }
			    elsif ($arr[5] > $threshold) { print PRR $arr[0]."\tEC2.4\n"; }
			    elsif ($arr[6] > $threshold) { print PRR $arr[0]."\tEC2.5\n"; }
			    elsif ($arr[7] > $threshold) { print PRR $arr[0]."\tEC2.6\n"; }
			    elsif ($arr[8] > $threshold) { print PRR $arr[0]."\tEC2.7\n"; }
			    elsif ($arr[9] > $threshold) { print PRR $arr[0]."\tEC2.8\n"; }
			    elsif ($arr[10] > $threshold){ print PRR $arr[0]."\tEC2.9\n"; }
			    else { print PRR $arr[0]."\tNP (Probability below thresold)\t$hash{$x}\n"; } # We have to decide wheteher this class should included or not
		    }
	  close PR; close PRR;
	  }
	                                                                   #--------------------------MACCSFP--------------------------------------------#
	  if ($model =~ /MACCSFP/)
	  {
		    system ("cut -d\",\" -f2129-2294 Hybrid_fingerprints.out  > MACCSFP.Out.Final");
		    #system ("rm Fingerprinter.out");

		    if ($sample =~ /Upsampling/) { system ("sh ec2RscriptMACCSFPUp.sh"); }
		    if ($sample =~ /Normal/) { system ("sh ec2RscriptMACCSFP.sh"); }
		    system ("sed -i '1d' MACCSFP_Pred2");

		    my $Rfpred = "MACCSFP_Pred2"; # prediction file with probabilities obtained from the random forest commnad.
		    $threshold = $Thr;        # $ARGV[1]; # cutoff thresold for the prediction 
		    open (PR, $Rfpred) || die "Cant open the file:\n";
		    open (PRR, ">$Rfpred.response2") || die "Cant open the file:\n";

		    while (chomp (my $line = <PR>))
		    {
		   #"EC2.1"	"EC2.10"	"EC2.2"	"EC2.3"	"EC2.4"	"EC2.5"	"EC2.6"	"EC2.7"	"EC2.8"	"EC2.9"
                  #"1"	0.334	0.02	0	0.166	0.208	0.136	0.014	0.094	0.022	0.006
			    
			    my @arr = split ("\t", $line);
			    
			    my %hash; $hash{$arr[1]}="EC2.1"; $hash{$arr[2]}="EC2.10"; $hash{$arr[3]}="EC2.2"; $hash{$arr[4]}="EC2.3"; $hash{$arr[5]}="EC2.4"; $hash{$arr[6]}="EC2.5"; $hash{$arr[7]}="EC2.6";
			    $hash{$arr[8]}="EC2.7"; $hash{$arr[9]}="EC2.8"; $hash{$arr[10]}="EC2.9"; 
			    
			    my $x = 0; foreach my $y (@arr){ $x = $y if $y > $x; }
			    
			    if ($arr[1] > $threshold) { print PRR $arr[0]."\tEC2.1\n"; }
			    elsif ($arr[2] > $threshold) { print PRR $arr[0]."\tEC2.10\n"; }
			    elsif ($arr[3] > $threshold) { print PRR $arr[0]."\tEC2.2\n"; }
			    elsif ($arr[4] > $threshold) { print PRR $arr[0]."\tEC2.3\n"; }
			    elsif ($arr[5] > $threshold) { print PRR $arr[0]."\tEC2.4\n"; }
			    elsif ($arr[6] > $threshold) { print PRR $arr[0]."\tEC2.5\n"; }
			    elsif ($arr[7] > $threshold) { print PRR $arr[0]."\tEC2.6\n"; }
			    elsif ($arr[8] > $threshold) { print PRR $arr[0]."\tEC2.7\n"; }
			    elsif ($arr[9] > $threshold) { print PRR $arr[0]."\tEC2.8\n"; }
			    elsif ($arr[10] > $threshold){ print PRR $arr[0]."\tEC2.9\n"; }
			    else { print PRR $arr[0]."\tNP (Probability below thresold)\t$hash{$x}\n"; } # We have to decide wheteher this class should included or not
		    }
	  close PR; close PRR;
	  }
	                                                      #-----------------------------------------subFPC-------------------------------------#
	  if ($model =~ /subFPC/)
	  {
		    system ("cut -d\",\" -f9123-9429 Hybrid_fingerprints.out  > subFPC.Out.Final");

		    if ($sample =~ /Upsampling/) { system ("sh ec2RscriptsubFPCUp.sh"); }
		    if ($sample =~ /Normal/) { system ("sh ec2RscriptsubFPC.sh"); }
		    system ("sed -i '1d' subFPC_Pred2");

		    my $Rfpred = "subFPC_Pred2"; # prediction file with probabilities obtained from the random forest commnad.
		    $threshold = $Thr;        # $ARGV[1]; # cutoff thresold for the prediction 
		    open (PR, $Rfpred) || die "Cant open the file:\n";
		    open (PRR, ">$Rfpred.response2") || die "Cant open the file:\n";

		    while (chomp (my $line = <PR>))
		    {
		   #"EC2.1"	"EC2.10"	"EC2.2"	"EC2.3"	"EC2.4"	"EC2.5"	"EC2.6"	"EC2.7"	"EC2.8"	"EC2.9"
                  #"1"	0.334	0.02	0	0.166	0.208	0.136	0.014	0.094	0.022	0.006
			    
			    my @arr = split ("\t", $line);
			    
			    my %hash; $hash{$arr[1]}="EC2.1"; $hash{$arr[2]}="EC2.10"; $hash{$arr[3]}="EC2.2"; $hash{$arr[4]}="EC2.3"; $hash{$arr[5]}="EC2.4"; $hash{$arr[6]}="EC2.5"; $hash{$arr[7]}="EC2.6";
			    $hash{$arr[8]}="EC2.7"; $hash{$arr[9]}="EC2.8"; $hash{$arr[10]}="EC2.9"; 

			    my $x = 0; foreach my $y (@arr){ $x = $y if $y > $x; }

			    if ($arr[1] > $threshold) { print PRR $arr[0]."\tEC2.1\n"; }
			    elsif ($arr[2] > $threshold) { print PRR $arr[0]."\tEC2.10\n"; }
			    elsif ($arr[3] > $threshold) { print PRR $arr[0]."\tEC2.2\n"; }
			    elsif ($arr[4] > $threshold) { print PRR $arr[0]."\tEC2.3\n"; }
			    elsif ($arr[5] > $threshold) { print PRR $arr[0]."\tEC2.4\n"; }
			    elsif ($arr[6] > $threshold) { print PRR $arr[0]."\tEC2.5\n"; }
			    elsif ($arr[7] > $threshold) { print PRR $arr[0]."\tEC2.6\n"; }
			    elsif ($arr[8] > $threshold) { print PRR $arr[0]."\tEC2.7\n"; }
			    elsif ($arr[9] > $threshold) { print PRR $arr[0]."\tEC2.8\n"; }
			    elsif ($arr[10] > $threshold){ print PRR $arr[0]."\tEC2.9\n"; }
			    else { print PRR $arr[0]."\tNP (Probability below thresold)\t$hash{$x}\n"; } # We have to decide wheteher this class should included or not
		    }
	  close PR; close PRR;
	  }
system ("rm *Out*"); 
}
#*-----------------------------------------------------------------------------------------------EC3-----------------------------------------------------------------------------------------------------------------------------#

if ($ec =~ /ec3/)
{
	                                                     #--------------------------------- Hybrid-----------------------------------#
	  if ($model =~ /Hybrid/)
	  {
		    system ("cut -d\",\" -f2-10209 Hybrid_fingerprints.out > Hybrid_fingerprints.Out");
		    system ("cut -d\",\" -f17,18,87,92,122,128,159,193,200,212,221,239,243,254,258,266,288,315,356,359,394,412,464,477,495,508,512,514,525,529,549,562,598,633,658,729,744,801,802,832,837,853,884,929,937,949,974,997,1002,1005,1007 Hybrid_fingerprints.Out > 1");
		    system ("cut -d\",\" -f1032,1035,1042,1045,1048,1057,1058,1060 Hybrid_fingerprints.Out > 2");
		    system ("cut -d\",\" -f1133,1195,1196,1253,1259,1269,1331,1345,1367,1389,1409,1410,1412,1431,1436,1465,1801,1825,1826,1828,1834,1855,1937,1989,2023,2085,2108 Hybrid_fingerprints.Out > 3");
		    system ("cut -d\",\" -f2140,2153,2161,2170,2176,2177,2180,2184,2211,2217,2225,2233,2237,2250,2252,2253,2263,2264,2267,2269,2273,2282,2285  Hybrid_fingerprints.Out > 4");
		    system ("cut -d\",\" -f2307,2309,2314,2315,2475,2480,2549,2617,2627,2659,2674,2699,2726,2730,2734,2745,2756,2822,2908,2953,2956,2958,2980,3112 Hybrid_fingerprints.Out > 5");
		    system ("cut -d\",\" -f3190,3192,3221,3222,3223,3249,3262,3317,3378,3449,3455,3456,3470 Hybrid_fingerprints.Out > 6");
		    system ("cut -d\",\" -f3489,3497,3670,4158,4887,5047,6028,6045,6601,6641,6749,6863,6869,7135,7149,7180,7210,7225,7270,7287,7301,7302,7363,7438,7496,7603,7674,7718,7776,8176,8292,8298,8302,8307 Hybrid_fingerprints.Out > 7");
		    system ("cut -d\",\" -f8343,8433,8443,8677,8679,8744,8911,8989 Hybrid_fingerprints.Out > 8");
		    system ("cut -d\",\" -f9124,9126,9137,9168,9169,9170,9209,9264,9325,9395,9396,9402,9403,9416,9417,9423 Hybrid_fingerprints.Out > 9");
		    system ("cut -d\",\" -f9429,9430,9431,9508,9520,9530,9585,9608,9667,9764,9831,9842,9910,9998,10076 Hybrid_fingerprints.Out > 10");

		    system ("paste -d\",\" 1 2 3 4 5 6 7 8 9 10 > Hybrid_fingerprints.Out.Final");
		    system ("rm 1 2 3 4 5 6 7 8 9 10 Hybrid_fingerprints.Out");

		    if ($sample =~ /Upsampling/) { system ("sh ec3RscriptHybridUp.sh"); }
		    if ($sample =~ /Normal/) { system ("sh ec3RscriptHybrid.sh"); }
		    system ("sed -i '1d' Hybrid_Pred2");

		    my $Rfpred = "Hybrid_Pred2"; # prediction file with probabilities obtained from the random forest commnad.
		    $threshold = $Thr;        # $ARGV[1]; # cutoff thresold for the prediction 
		    open (PR, $Rfpred) || die "Cant open the file:\n";
		    open (PRR, ">$Rfpred.response2") || die "Cant open the file:\n";

		    while (chomp (my $line = <PR>))
		    {
	  #	       "EC3.1"	"EC3.11"	"EC3.2"	"EC3.3"	"EC3.4"	"EC3.5"	"EC3.6"	"EC3.7"	"EC3.8"
#            "1"	0.654	0	0.15	0.05	0.008	0.094	0.006	0.038	0
			    
			    my @arr = split ("\t", $line);
			    
			    my %hash; $hash{$arr[1]}="EC3.1"; $hash{$arr[2]}="EC3.11"; $hash{$arr[3]}="EC3.2"; $hash{$arr[4]}="EC3.3"; $hash{$arr[5]}="EC3.4"; $hash{$arr[6]}="EC3.5"; $hash{$arr[7]}="EC3.6";
			    $hash{$arr[8]}="EC3.7"; $hash{$arr[9]}="EC3.8";

			    my $x = 0; foreach my $y (@arr){ $x = $y if $y > $x; }
			    
			    if ($arr[1] > $threshold) { print PRR $arr[0]."\tEC3.1\n"; }
			    elsif ($arr[2] > $threshold) { print PRR $arr[0]."\tEC3.11\n"; }
			    elsif ($arr[3] > $threshold) { print PRR $arr[0]."\tEC3.2\n";  }
			    elsif ($arr[4] > $threshold) { print PRR $arr[0]."\tEC3.3\n";  }
			    elsif ($arr[5] > $threshold) { print PRR $arr[0]."\tEC3.4\n";  }
			    elsif ($arr[6] > $threshold) { print PRR $arr[0]."\tEC3.5\n";  }
			    elsif ($arr[7] > $threshold) { print PRR $arr[0]."\tEC3.6\n";  }
			    elsif ($arr[8] > $threshold) { print PRR $arr[0]."\tEC3.7\n";  }
			    elsif ($arr[9] > $threshold) { print PRR $arr[0]."\tEC3.8\n";  }
			    else { print PRR $arr[0]."\tNP (Probability below thresold)\t$hash{$x}\n"; }# We have to decide wheteher this class should included or not
		    }
	  close PR; close PRR;
	  }
	                                                                   #--------------------------Fingerprinter--------------------------------------------#
	  if ($model =~ /Fingerprinter/)
	  {
		    system ("cut -d\",\" -f2-1025 Hybrid_fingerprints.out  > Fingerprinter.Out.Final");
		    #system ("rm Fingerprinter.out");

		    if ($sample =~ /Upsampling/) { system ("sh ec3RscriptFingerprinterUp.sh"); }
		    if ($sample =~ /Normal/)     { system ("sh ec3RscriptFingerprinter.sh");   }
		    system ("sed -i '1d' Fingerprinter_Pred2");

		    my $Rfpred = "Fingerprinter_Pred2"; # prediction file with probabilities obtained from the random forest commnad.
		    $threshold = $Thr;        # $ARGV[1]; # cutoff thresold for the prediction 
		    open (PR, $Rfpred) || die "Cant open the file:\n";
		    open (PRR, ">$Rfpred.response2") || die "Cant open the file:\n";

		    while (chomp (my $line = <PR>))
		    {
	         	#      "EC3.1"	"EC3.11"	"EC3.2"	"EC3.3"	"EC3.4"	"EC3.5"	"EC3.6"	"EC3.7"	"EC3.8"
#                      "1"	0.654	0	0.15	0.05	0.008	0.094	0.006	0.038	0
			    
			    my @arr = split ("\t", $line);
			    
			    my %hash; $hash{$arr[1]}="EC3.1"; $hash{$arr[2]}="EC3.11"; $hash{$arr[3]}="EC3.2"; $hash{$arr[4]}="EC3.3"; $hash{$arr[5]}="EC3.4"; $hash{$arr[6]}="EC3.5"; $hash{$arr[7]}="EC3.6";
			    $hash{$arr[8]}="EC3.7"; $hash{$arr[9]}="EC3.8";
			    
			    my $x = 0; foreach my $y (@arr){ $x = $y if $y > $x; }
			    
			    if ($arr[1] > $threshold) { print PRR $arr[0]."\tEC3.1\n"; }
			    elsif ($arr[2] > $threshold) { print PRR $arr[0]."\tEC3.11\n"; }
			    elsif ($arr[3] > $threshold) { print PRR $arr[0]."\tEC3.2\n";  }
			    elsif ($arr[4] > $threshold) { print PRR $arr[0]."\tEC3.3\n";  }
			    elsif ($arr[5] > $threshold) { print PRR $arr[0]."\tEC3.4\n";  }
			    elsif ($arr[6] > $threshold) { print PRR $arr[0]."\tEC3.5\n";  }
			    elsif ($arr[7] > $threshold) { print PRR $arr[0]."\tEC3.6\n";  }
			    elsif ($arr[8] > $threshold) { print PRR $arr[0]."\tEC3.7\n";  }
			    elsif ($arr[9] > $threshold) { print PRR $arr[0]."\tEC3.8\n";  }
			    else { print PRR $arr[0]."\tNP (Probability below thresold)\t$hash{$x}\n"; }# We have to decide wheteher this class should included or not
		    }
	  close PR; close PRR;
	  }
	                                                      #-----------------------------------------GraphFP-------------------------------------#
	  if ($model =~ /GraphFP/)
	  {
		    system ("cut -d\",\" -f1105-2128 Hybrid_fingerprints.out  > GraphFP.Out.Final");

		    if ($sample =~ /Upsampling/) { system ("sh ec3RscriptGraphFPUp.sh"); }
		    if ($sample =~ /Normal/)     { system ("sh ec3RscriptGraphFP.sh");   }
		    system ("sed -i '1d' GraphFP_Pred2");

		    my $Rfpred = "GraphFP_Pred2"; # prediction file with probabilities obtained from the random forest commnad.
		    $threshold = $Thr;        # $ARGV[1]; # cutoff thresold for the prediction 
		    open (PR, $Rfpred) || die "Cant open the file:\n";
		    open (PRR, ">$Rfpred.response2") || die "Cant open the file:\n";

		    while (chomp (my $line = <PR>))
		    {
		      #        "EC3.1"	"EC3.11"	"EC3.2"	"EC3.3"	"EC3.4"	"EC3.5"	"EC3.6"	"EC3.7"	"EC3.8"
#                     "1"	0.654	0	0.15	0.05	0.008	0.094	0.006	0.038	0
	    
			    my @arr = split ("\t", $line);
			    
			    my %hash; $hash{$arr[1]}="EC3.1"; $hash{$arr[2]}="EC3.11"; $hash{$arr[3]}="EC3.2"; $hash{$arr[4]}="EC3.3"; $hash{$arr[5]}="EC3.4"; $hash{$arr[6]}="EC3.5"; $hash{$arr[7]}="EC3.6";
			    $hash{$arr[8]}="EC3.7"; $hash{$arr[9]}="EC3.8";
			    
			    my $x = 0; foreach my $y (@arr){ $x = $y if $y > $x; }

			    if ($arr[1] > $threshold) { print PRR $arr[0]."\tEC3.1\n"; }
			    elsif ($arr[2] > $threshold) { print PRR $arr[0]."\tEC3.11\n"; }
			    elsif ($arr[3] > $threshold) { print PRR $arr[0]."\tEC3.2\n";  }
			    elsif ($arr[4] > $threshold) { print PRR $arr[0]."\tEC3.3\n";  }
			    elsif ($arr[5] > $threshold) { print PRR $arr[0]."\tEC3.4\n";  }
			    elsif ($arr[6] > $threshold) { print PRR $arr[0]."\tEC3.5\n";  }
			    elsif ($arr[7] > $threshold) { print PRR $arr[0]."\tEC3.6\n";  }
			    elsif ($arr[8] > $threshold) { print PRR $arr[0]."\tEC3.7\n";  }
			    elsif ($arr[9] > $threshold) { print PRR $arr[0]."\tEC3.8\n";  }
			    else { print PRR $arr[0]."\tNP (Probability below thresold)\t$hash{$x}\n"; }# We have to decide wheteher this class should included or not
		    }
	  close PR; close PRR;
	  }
system ("rm *Out*"); 
}
#*-----------------------------------------------------------------------------------------------EC4-----------------------------------------------------------------------------------------------------------------------------#

if ($ec =~ /ec4/)
{
	                                                     #--------------------------------- Hybrid-----------------------------------#
	  if ($model =~ /Hybrid/)
	  {
		    system ("cut -d\",\" -f2-10209 Hybrid_fingerprints.out > Hybrid_fingerprints.Out");    
		    system ("cut -d\",\" -f17,18,87,92,122,128,159,193,200,212,221,239,243,254,258,266,288,315,356,359,394,412,464,477,495,508,512,514,525,529,549,562,598,633,658,729,744,801,802,832,837,853,884,929,937,949,974,997,1002,1005,1007 Hybrid_fingerprints.Out > 1");
		    system ("cut -d\",\" -f1032,1035,1042,1045,1048,1057,1058,1060 Hybrid_fingerprints.Out > 2");
		    system ("cut -d\",\" -f1133,1195,1196,1253,1259,1269,1331,1345,1367,1389,1409,1410,1412,1431,1436,1465,1801,1825,1826,1828,1834,1855,1937,1989,2023,2085,2108 Hybrid_fingerprints.Out > 3");
		    system ("cut -d\",\" -f2140,2153,2161,2170,2176,2177,2180,2184,2211,2217,2225,2233,2237,2250,2252,2253,2263,2264,2267,2269,2273,2282,2285  Hybrid_fingerprints.Out > 4");
		    system ("cut -d\",\" -f2307,2309,2314,2315,2475,2480,2549,2617,2627,2659,2674,2699,2726,2730,2734,2745,2756,2822,2908,2953,2956,2958,2980,3112 Hybrid_fingerprints.Out > 5");
		    system ("cut -d\",\" -f3190,3192,3221,3222,3223,3249,3262,3317,3378,3449,3455,3456,3470 Hybrid_fingerprints.Out > 6");
		    system ("cut -d\",\" -f3489,3497,3670,4158,4887,5047,6028,6045,6601,6641,6749,6863,6869,7135,7149,7180,7210,7225,7270,7287,7301,7302,7363,7438,7496,7603,7674,7718,7776,8176,8292,8298,8302,8307 Hybrid_fingerprints.Out > 7");
		    system ("cut -d\",\" -f8343,8433,8443,8677,8679,8744,8911,8989 Hybrid_fingerprints.Out > 8");
		    system ("cut -d\",\" -f9124,9126,9137,9168,9169,9170,9209,9264,9325,9395,9396,9402,9403,9416,9417,9423 Hybrid_fingerprints.Out > 9");
		    system ("cut -d\",\" -f9429,9430,9431,9508,9520,9530,9585,9608,9667,9764,9831,9842,9910,9998,10076 Hybrid_fingerprints.Out > 10");

		    system ("paste -d\",\" 1 2 3 4 5 6 7 8 9 10 > Hybrid_fingerprints.Out.Final");
		    system ("rm 1 2 3 4 5 6 7 8 9 10 Hybrid_fingerprints.Out");

		    if ($sample =~ /Upsampling/)  { system ("sh ec4RscriptHybridUp.sh"); }
		    if ($sample =~ /Normal/)      { system ("sh ec4RscriptHybrid.sh");   }
		    system ("sed -i '1d' Hybrid_Pred2");

		    my $Rfpred = "Hybrid_Pred2"; # prediction file with probabilities obtained from the random forest commnad.
		    $threshold = $Thr;        # $ARGV[1]; # cutoff thresold for the prediction 
		    open (PR, $Rfpred) || die "Cant open the file:\n";
		    open (PRR, ">$Rfpred.response2") || die "Cant open the file:\n";

		    while (chomp (my $line = <PR>))
		    {
	  #	      "EC4.1"	"EC4.2"	"EC4.3"	"EC4.4"	"EC4.5"	"EC4.6"	"EC4.7"	"EC4.99"
#                "1"	0.142	0.452	0.174	0.106	0.008	0.05	0.012	0.056
		    
			    my @arr = split ("\t", $line);

			    my %hash; $hash{$arr[1]}="EC4.1"; $hash{$arr[2]}="EC4.2"; $hash{$arr[3]}="EC4.3"; $hash{$arr[4]}="EC4.4"; $hash{$arr[5]}="EC4.5"; $hash{$arr[6]}="EC4.6"; $hash{$arr[7]}="EC4.7";
			    $hash{$arr[8]}="EC4.99";
			    
			    my $x = 0; foreach my $y (@arr){ $x = $y if $y > $x; }

			    if ($arr[1] > $threshold) { print PRR $arr[0]."\tEC4.1\n";}
			    elsif ($arr[2] > $threshold){ print PRR $arr[0]."\tEC4.2\n"; }
			    elsif ($arr[3] > $threshold){ print PRR $arr[0]."\tEC4.3\n"; }
			    elsif ($arr[4] > $threshold){ print PRR $arr[0]."\tEC4.4\n"; }
			    elsif ($arr[5] > $threshold){ print PRR $arr[0]."\tEC4.5\n"; }
			    elsif ($arr[6] > $threshold){ print PRR $arr[0]."\tEC4.6\n"; }
			    elsif ($arr[7] > $threshold){ print PRR $arr[0]."\tEC4.7\n"; }
			    elsif ($arr[8] > $threshold){ print PRR $arr[0]."\tEC4.99\n";}
			    else { print PRR $arr[0]."\tNP (Probability below thresold)\t$hash{$x}\n"; } # We have to decide wheteher this class should included or not
		    }
	  close PR; close PRR;
	  }
	                                                                   #--------------------------KRFP --------------------------------------------#
	  if ($model =~ /KRFP/)
	  {
		    system ("cut -d\",\" -f3483-8342 Hybrid_fingerprints.out > KRFP.Out.Final");

		    if ($sample =~ /Upsampling/) { system ("sh ec4RscriptKRFPUp.sh"); }
		    if ($sample =~ /Normal/)     { system ("sh ec4RscriptKRFP.sh");   }
		    system ("sed -i '1d' KRFP_Pred2");

		    my $Rfpred = "KRFP_Pred2"; # prediction file with probabilities obtained from the random forest commnad.
		    $threshold = $Thr;        # $ARGV[1]; # cutoff thresold for the prediction 
		    open (PR, $Rfpred) || die "Cant open the file:\n";
		    open (PRR, ">$Rfpred.response2") || die "Cant open the file:\n";

		    while (chomp (my $line = <PR>))
		    {
		#       "EC4.1"	"EC4.2"	"EC4.3"	"EC4.4"	"EC4.5"	"EC4.6"	"EC4.7"	"EC4.99"
#                   "1"	0.142	0.452	0.174	0.106	0.008	0.05	0.012	0.056
		    
			    my @arr = split ("\t", $line);
			    
			    my %hash; $hash{$arr[1]}="EC4.1"; $hash{$arr[2]}="EC4.2"; $hash{$arr[3]}="EC4.3"; $hash{$arr[4]}="EC4.4"; $hash{$arr[5]}="EC4.5"; $hash{$arr[6]}="EC4.6"; $hash{$arr[7]}="EC4.7";
			    $hash{$arr[8]}="EC4.99";
			    
			    my $x = 0; foreach my $y (@arr){ $x = $y if $y > $x; }
			    
			    if ($arr[1] > $threshold) { print PRR $arr[0]."\tEC4.1\n";}
			    elsif ($arr[2] > $threshold){ print PRR $arr[0]."\tEC4.2\n"; }
			    elsif ($arr[3] > $threshold){ print PRR $arr[0]."\tEC4.3\n"; }
			    elsif ($arr[4] > $threshold){ print PRR $arr[0]."\tEC4.4\n"; }
			    elsif ($arr[5] > $threshold){ print PRR $arr[0]."\tEC4.5\n"; }
			    elsif ($arr[6] > $threshold){ print PRR $arr[0]."\tEC4.6\n"; }
			    elsif ($arr[7] > $threshold){ print PRR $arr[0]."\tEC4.7\n"; }
			    elsif ($arr[8] > $threshold){ print PRR $arr[0]."\tEC4.99\n";}
			    else { print PRR $arr[0]."\tNP (Probability below thresold)\t$hash{$x}\n"; } # We have to decide wheteher this class should included or not
		    }
	  close PR; close PRR;
	  }
	                                                      #-----------------------------------------Pubchem-------------------------------------#
	  if ($model =~ /Pubchem/)
	  {
		    system ("cut -d\",\" -f2295-3175 Hybrid_fingerprints.out  > Pubchem.Out.Final");

		    if ($sample =~ /Upsampling/) { system ("sh ec4RscriptPubchemUp.sh"); }
		    if ($sample =~ /Normal/)     { system ("sh ec4RscriptPubchem.sh");   }
		    system ("sed -i '1d' Pubchem_Pred2");

		    my $Rfpred = "Pubchem_Pred2"; # prediction file with probabilities obtained from the random forest commnad.
		    $threshold = $Thr;        # $ARGV[1]; # cutoff thresold for the prediction 
		    open (PR, $Rfpred) || die "Cant open the file:\n";
		    open (PRR, ">$Rfpred.response2") || die "Cant open the file:\n";

		    while (chomp (my $line = <PR>))
		    {
		#    "EC4.1"	"EC4.2"	"EC4.3"	"EC4.4"	"EC4.5"	"EC4.6"	"EC4.7"	"EC4.99"
#                "1"	0.142	0.452	0.174	0.106	0.008	0.05	0.012	0.056
		    
			    my @arr = split ("\t", $line);
			    
			    my %hash; $hash{$arr[1]}="EC4.1"; $hash{$arr[2]}="EC4.2"; $hash{$arr[3]}="EC4.3"; $hash{$arr[4]}="EC4.4"; $hash{$arr[5]}="EC4.5"; $hash{$arr[6]}="EC4.6"; $hash{$arr[7]}="EC4.7";
			    $hash{$arr[8]}="EC4.99";
			    
			    my $x = 0; foreach my $y (@arr){ $x = $y if $y > $x; }

			    if ($arr[1] > $threshold) { print PRR $arr[0]."\tEC4.1\n";}
			    elsif ($arr[2] > $threshold){ print PRR $arr[0]."\tEC4.2\n"; }
			    elsif ($arr[3] > $threshold){ print PRR $arr[0]."\tEC4.3\n"; }
			    elsif ($arr[4] > $threshold){ print PRR $arr[0]."\tEC4.4\n"; }
			    elsif ($arr[5] > $threshold){ print PRR $arr[0]."\tEC4.5\n"; }
			    elsif ($arr[6] > $threshold){ print PRR $arr[0]."\tEC4.6\n"; }
			    elsif ($arr[7] > $threshold){ print PRR $arr[0]."\tEC4.7\n"; }
			    elsif ($arr[8] > $threshold){ print PRR $arr[0]."\tEC4.99\n";}
			    else { print PRR $arr[0]."\tNP (Probability below thresold)\t$hash{$x}\n"; } # We have to decide wheteher this class should included or not
		    }
	  close PR; close PRR;
	  }
system ("rm *Out*"); 
}

#*-----------------------------------------------------------------------------------------------EC5-----------------------------------------------------------------------------------------------------------------------------#

if ($ec =~ /ec5/)
{	                                                    
                                                                #--------------------------------- Hybrid-----------------------------------#
	  if ($model =~ /Hybrid/)
	  {
		    system ("cut -d\",\" -f2-10209 Hybrid_fingerprints.out > Hybrid_fingerprints.Out");	    
		    system ("cut -d\",\" -f17,18,87,92,122,128,159,193,200,212,221,239,243,254,258,266,288,315,356,359,394,412,464,477,495,508,512,514,525,529,549,562,598,633,658,729,744,801,802,832,837,853,884,929,937,949,974,997,1002,1005,1007 Hybrid_fingerprints.Out > 1");
		    system ("cut -d\",\" -f1032,1035,1042,1045,1048,1057,1058,1060 Hybrid_fingerprints.Out > 2");
		    system ("cut -d\",\" -f1133,1195,1196,1253,1259,1269,1331,1345,1367,1389,1409,1410,1412,1431,1436,1465,1801,1825,1826,1828,1834,1855,1937,1989,2023,2085,2108 Hybrid_fingerprints.Out > 3");
		    system ("cut -d\",\" -f2140,2153,2161,2170,2176,2177,2180,2184,2211,2217,2225,2233,2237,2250,2252,2253,2263,2264,2267,2269,2273,2282,2285  Hybrid_fingerprints.Out > 4");
		    system ("cut -d\",\" -f2307,2309,2314,2315,2475,2480,2549,2617,2627,2659,2674,2699,2726,2730,2734,2745,2756,2822,2908,2953,2956,2958,2980,3112 Hybrid_fingerprints.Out > 5");
		    system ("cut -d\",\" -f3190,3192,3221,3222,3223,3249,3262,3317,3378,3449,3455,3456,3470 Hybrid_fingerprints.Out > 6");
		    system ("cut -d\",\" -f3489,3497,3670,4158,4887,5047,6028,6045,6601,6641,6749,6863,6869,7135,7149,7180,7210,7225,7270,7287,7301,7302,7363,7438,7496,7603,7674,7718,7776,8176,8292,8298,8302,8307 Hybrid_fingerprints.Out > 7");
		    system ("cut -d\",\" -f8343,8433,8443,8677,8679,8744,8911,8989 Hybrid_fingerprints.Out > 8");
		    system ("cut -d\",\" -f9124,9126,9137,9168,9169,9170,9209,9264,9325,9395,9396,9402,9403,9416,9417,9423 Hybrid_fingerprints.Out > 9");
		    system ("cut -d\",\" -f9429,9430,9431,9508,9520,9530,9585,9608,9667,9764,9831,9842,9910,9998,10076 Hybrid_fingerprints.Out > 10");

		    system ("paste -d\",\" 1 2 3 4 5 6 7 8 9 10 > Hybrid_fingerprints.Out.Final");
		    system ("rm 1 2 3 4 5 6 7 8 9 10 Hybrid_fingerprints.Out");

		    if ($sample =~ /Upsampling/) { system ("sh ec5RscriptHybridUp.sh");}
		    if ($sample =~ /Normal/)     { system ("sh ec5RscriptHybrid.sh");  }
		    system ("sed -i '1d' Hybrid_Pred2");

		    my $Rfpred = "Hybrid_Pred2"; # prediction file with probabilities obtained from the random forest commnad.
		    $threshold = $Thr;        # $ARGV[1]; # cutoff thresold for the prediction 
		    open (PR, $Rfpred) || die "Cant open the file:\n";
		    open (PRR, ">$Rfpred.response2") || die "Cant open the file:\n";

		    while (chomp (my $line = <PR>))
		    {
	  #	        EC5.1"	"EC5.2"	"EC5.3"	"EC5.4"	"EC5.5"
  #          "1"	0.418	0.104	0.314	0.108	0.056
		    
			    my @arr = split ("\t", $line);
			     
			    my %hash; $hash{$arr[1]}="EC5.1"; $hash{$arr[2]}="EC5.2"; $hash{$arr[3]}="EC5.3"; $hash{$arr[4]}="EC5.4"; $hash{$arr[5]}="EC5.5";
			    
			    my $x = 0; foreach my $y (@arr){ $x = $y if $y > $x; }
			    
			    if ($arr[1] > $threshold) { print PRR $arr[0]."\tEC5.1\n";}
			    elsif ($arr[2] > $threshold) { print PRR $arr[0]."\tEC5.2\n";}
			    elsif ($arr[3] > $threshold) { print PRR $arr[0]."\tEC5.3\n";}
			    elsif ($arr[4] > $threshold) { print PRR $arr[0]."\tEC5.4\n";}
			    elsif ($arr[5] > $threshold) { print PRR $arr[0]."\tEC5.5\n";}
			    else { print PRR $arr[0]."\tNP (Probability below thresold)\t$hash{$x}\n"; } # We have to decide wheteher this class should included or not
		    }
	  close PR; close PRR;
	  }
                                                                   #--------------------------KRFP --------------------------------------------#
	  if ($model =~ /KRFP/)
	  {
		    system ("cut -d\",\" -f3483-8342 Hybrid_fingerprints.out > KRFP.Out.Final");
		    #system ("rm Fingerprinter.out");

		    if ($sample =~ /Upsampling/) { system ("sh ec5RscriptKRFPUp.sh");}
		    if ($sample =~ /Normal/)     { system ("sh ec5RscriptKRFP.sh");  }
		    system ("sed -i '1d' KRFP_Pred2");

		    my $Rfpred = "KRFP_Pred2"; # prediction file with probabilities obtained from the random forest commnad.
		    $threshold = $Thr;        # $ARGV[1]; # cutoff thresold for the prediction 
		    open (PR, $Rfpred) || die "Cant open the file:\n";
		    open (PRR, ">$Rfpred.response2") || die "Cant open the file:\n";

		    while (chomp (my $line = <PR>))
		    {
		#     EC5.1"	"EC5.2"	"EC5.3"	"EC5.4"	"EC5.5"
  #          "1"	0.418	0.104	0.314	0.108	0.056
		    
			    my @arr = split ("\t", $line);
			    
			    my %hash; $hash{$arr[1]}="EC5.1"; $hash{$arr[2]}="EC5.2"; $hash{$arr[3]}="EC5.3"; $hash{$arr[4]}="EC5.4"; $hash{$arr[5]}="EC5.5";
			    
			    my $x = 0; foreach my $y (@arr){ $x = $y if $y > $x; }
			    
			    if ($arr[1] > $threshold) { print PRR $arr[0]."\tEC5.1\n";}
			    elsif ($arr[2] > $threshold) { print PRR $arr[0]."\tEC5.2\n";}
			    elsif ($arr[3] > $threshold) { print PRR $arr[0]."\tEC5.3\n";}
			    elsif ($arr[4] > $threshold) { print PRR $arr[0]."\tEC5.4\n";}
			    elsif ($arr[5] > $threshold) { print PRR $arr[0]."\tEC5.5\n";}
			    else { print PRR $arr[0]."\tNP (Probability below thresold)\t$hash{$x}\n"; } # We have to decide wheteher this class should included or not
		    }
	  close PR; close PRR;
	  }
	                                                      #-----------------------------------------Fingerprinter-------------------------------------#
	  if ($model =~ /Fingerprinter/)
	  {

		    system ("cut -d\",\" -f2-1025 Hybrid_fingerprints.out  > Fingerprinter.Out.Final");

		    if ($sample =~ /Upsampling/) { system ("sh ec5RscriptFingerprinterUp.sh"); }
		    if ($sample =~ /Normal/)     { system ("sh ec5RscriptFingerprinter.sh");  }
		    system ("sed -i '1d' Fingerprinter_Pred2");

		    my $Rfpred = "Fingerprinter_Pred2"; # prediction file with probabilities obtained from the random forest commnad.
		    $threshold = $Thr;        # $ARGV[1]; # cutoff thresold for the prediction 
		    open (PR, $Rfpred) || die "Cant open the file:\n";
		    open (PRR, ">$Rfpred.response2") || die "Cant open the file:\n";

		    while (chomp (my $line = <PR>))
		    {
		#   EC5.1"	"EC5.2"	"EC5.3"	"EC5.4"	"EC5.5"
  #          "1"	0.418	0.104	0.314	0.108	0.056
		    
			    my @arr = split ("\t", $line);
			    
			    my %hash; $hash{$arr[1]}="EC5.1"; $hash{$arr[2]}="EC5.2"; $hash{$arr[3]}="EC5.3"; $hash{$arr[4]}="EC5.4"; $hash{$arr[5]}="EC5.5";
			    
			    my $x = 0; foreach my $y (@arr){ $x = $y if $y > $x; }
			    
			    if ($arr[1] > $threshold) { print PRR $arr[0]."\tEC5.1\n";}
			    elsif ($arr[2] > $threshold) { print PRR $arr[0]."\tEC5.2\n";}
			    elsif ($arr[3] > $threshold) { print PRR $arr[0]."\tEC5.3\n";}
			    elsif ($arr[4] > $threshold) { print PRR $arr[0]."\tEC5.4\n";}
			    elsif ($arr[5] > $threshold) { print PRR $arr[0]."\tEC5.5\n";}
			    else { print PRR $arr[0]."\tNP (Probability below thresold)\t$hash{$x}\n"; } # We have to decide wheteher this class should included or not
		    }
	  close PR; close PRR;
	  }
system ("rm *Out*"); 
}

#*-----------------------------------------------------------------------------------------------EC6----------------------------------------------------------------------------------------------------------------------------#

if ($ec =~ /ec6/)
{
	                                                     #--------------------------------- Hybrid-----------------------------------#
	  if ($model =~ /Hybrid/)
	  {
		    system ("cut -d\",\" -f2-10209 Hybrid_fingerprints.out > Hybrid_fingerprints.Out");		    
		    system ("cut -d\",\" -f17,18,87,92,122,128,159,193,200,212,221,239,243,254,258,266,288,315,356,359,394,412,464,477,495,508,512,514,525,529,549,562,598,633,658,729,744,801,802,832,837,853,884,929,937,949,974,997,1002,1005,1007 Hybrid_fingerprints.Out > 1");
		    system ("cut -d\",\" -f1032,1035,1042,1045,1048,1057,1058,1060 Hybrid_fingerprints.Out > 2");
		    system ("cut -d\",\" -f1133,1195,1196,1253,1259,1269,1331,1345,1367,1389,1409,1410,1412,1431,1436,1465,1801,1825,1826,1828,1834,1855,1937,1989,2023,2085,2108 Hybrid_fingerprints.Out > 3");
		    system ("cut -d\",\" -f2140,2153,2161,2170,2176,2177,2180,2184,2211,2217,2225,2233,2237,2250,2252,2253,2263,2264,2267,2269,2273,2282,2285  Hybrid_fingerprints.Out > 4");
		    system ("cut -d\",\" -f2307,2309,2314,2315,2475,2480,2549,2617,2627,2659,2674,2699,2726,2730,2734,2745,2756,2822,2908,2953,2956,2958,2980,3112 Hybrid_fingerprints.Out > 5");
		    system ("cut -d\",\" -f3190,3192,3221,3222,3223,3249,3262,3317,3378,3449,3455,3456,3470 Hybrid_fingerprints.Out > 6");
		    system ("cut -d\",\" -f3489,3497,3670,4158,4887,5047,6028,6045,6601,6641,6749,6863,6869,7135,7149,7180,7210,7225,7270,7287,7301,7302,7363,7438,7496,7603,7674,7718,7776,8176,8292,8298,8302,8307 Hybrid_fingerprints.Out > 7");
		    system ("cut -d\",\" -f8343,8433,8443,8677,8679,8744,8911,8989 Hybrid_fingerprints.Out > 8");
		    system ("cut -d\",\" -f9124,9126,9137,9168,9169,9170,9209,9264,9325,9395,9396,9402,9403,9416,9417,9423 Hybrid_fingerprints.Out > 9");
		    system ("cut -d\",\" -f9429,9430,9431,9508,9520,9530,9585,9608,9667,9764,9831,9842,9910,9998,10076 Hybrid_fingerprints.Out > 10");

		    system ("paste -d\",\" 1 2 3 4 5 6 7 8 9 10 > Hybrid_fingerprints.Out.Final");
		    system ("rm 1 2 3 4 5 6 7 8 9 10 Hybrid_fingerprints.Out");

		    if ($sample =~ /Upsampling/) { system ("sh ec6RscriptHybridUp.sh"); }
		    if ($sample =~ /Normal/)     { system ("sh ec6RscriptHybrid.sh");   }
		    system ("sed -i '1d' Hybrid_Pred2");

		    my $Rfpred = "Hybrid_Pred2"; # prediction file with probabilities obtained from the random forest commnad.
		    $threshold = $Thr;        # $ARGV[1]; # cutoff thresold for the prediction 
		    open (PR, $Rfpred) || die "Cant open the file:\n";
		    open (PRR, ">$Rfpred.response2") || die "Cant open the file:\n";

		    while (chomp (my $line = <PR>))
		    {
	  #	       "EC6.2"	"EC6.3"	"EC6.4"	"EC6.5"	"EC6.6"
#            "1"	0.114	0.734	0.044	0.006	0.102
		    
			    my @arr = split ("\t", $line);
			    
			    my %hash; $hash{$arr[1]}="EC6.2"; $hash{$arr[2]}="EC6.3"; $hash{$arr[3]}="EC6.4"; $hash{$arr[4]}="EC6.5"; $hash{$arr[5]}="EC6.6";
			    
			    my $x = 0; foreach my $y (@arr){ $x = $y if $y > $x; }
			    
			    if ($arr[1] > $threshold) { print PRR $arr[0]."\tEC6.2\n";}
			    elsif ($arr[2] > $threshold){ print PRR $arr[0]."\tEC6.3\n";}
			    elsif ($arr[3] > $threshold) { print PRR $arr[0]."\tEC6.4\n";}
			    elsif ($arr[4] > $threshold) { print PRR $arr[0]."\tEC6.5\n";}
			    elsif ($arr[5] > $threshold) { print PRR $arr[0]."\tEC6.6\n";}
			    else { print PRR $arr[0]."\tNP (Probability below thresold)\t$hash{$x}\n"; } # We have to decide wheteher this class should included or not
		    }
	  close PR; close PRR;
	  }

	                                                      #-----------------------------------------Fingerprinter-------------------------------------#
	  if ($model =~ /Fingerprinter/)
	  {
		    system ("cut -d\",\" -f2-1025 Hybrid_fingerprints.out  > Fingerprinter.Out.Final");

		    if ($sample =~ /Upsampling/) { system ("sh ec6RscriptFingerprinterUp.sh"); }
		    if ($sample =~ /Normal/)    { system ("sh ec6RscriptFingerprinter.sh"); }
		    system ("sed -i '1d' Fingerprinter_Pred2");

		    my $Rfpred = "Fingerprinter_Pred2"; # prediction file with probabilities obtained from the random forest commnad.
		    $threshold = $Thr;        # $ARGV[1]; # cutoff thresold for the prediction 
		    open (PR, $Rfpred) || die "Cant open the file:\n";
		    open (PRR, ">$Rfpred.response2") || die "Cant open the file:\n";

		    while (chomp (my $line = <PR>))
		    {
		#      "EC6.2"	"EC6.3"	"EC6.4"	"EC6.5"	"EC6.6"
#                "1"	0.114	0.734	0.044	0.006	0.102
		    
			    my @arr = split ("\t", $line);
			    
			    my %hash; $hash{$arr[1]}="EC6.2"; $hash{$arr[2]}="EC6.3"; $hash{$arr[3]}="EC6.4"; $hash{$arr[4]}="EC6.5"; $hash{$arr[5]}="EC6.6";
			    
			    my $x = 0; foreach my $y (@arr){ $x = $y if $y > $x; }
			    
			    if ($arr[1] > $threshold) { print PRR $arr[0]."\tEC6.2\n";}
			    elsif ($arr[2] > $threshold){ print PRR $arr[0]."\tEC6.3\n";}
			    elsif ($arr[3] > $threshold) { print PRR $arr[0]."\tEC6.4\n";}
			    elsif ($arr[4] > $threshold) { print PRR $arr[0]."\tEC6.5\n";}
			    elsif ($arr[5] > $threshold) { print PRR $arr[0]."\tEC6.6\n";}
			    else { print PRR $arr[0]."\tNP (Probability below thresold)\t$hash{$x}\n"; } # We have to decide wheteher this class should included or not
		    }
	  close PR; close PRR;
	  }
                                                             #--------------------------APC2D---------------------------------------------#
	  if ($model =~ /APC2D/)
	  {
		    system ("cut -d\",\" -f9430-10209 Hybrid_fingerprints.out > APC2D.Out.Final");

		    if ($sample =~ /Upsampling/) { system ("sh ec6RscriptAPC2DUp.sh");}
		    if ($sample =~ /Normal/)  { system ("sh ec6RscriptAPC2D.sh"); }
		    system ("sed -i '1d' APC2D_Pred2");

		    my $Rfpred = "APC2D_Pred2"; # prediction file with probabilities obtained from the random forest commnad.
		    $threshold = $Thr;        # $ARGV[1]; # cutoff thresold for the prediction 
		    open (PR, $Rfpred) || die "Cant open the file:\n";
		    open (PRR, ">$Rfpred.response2") || die "Cant open the file:\n";

		    while (chomp (my $line = <PR>))
		    {
		#   "EC6.2"	"EC6.3"	"EC6.4"	"EC6.5"	"EC6.6"
#            "1"	0.114	0.734	0.044	0.006	0.102
		    
			    my @arr = split ("\t", $line);
			    
			     my %hash; $hash{$arr[1]}="EC6.2"; $hash{$arr[2]}="EC6.3"; $hash{$arr[3]}="EC6.4"; $hash{$arr[4]}="EC6.5"; $hash{$arr[5]}="EC6.6";
			    
			    my $x = 0; foreach my $y (@arr){ $x = $y if $y > $x; }
			    
			    if ($arr[1] > $threshold) { print PRR $arr[0]."\tEC6.2\n";}
			    elsif ($arr[2] > $threshold){ print PRR $arr[0]."\tEC6.3\n";}
			    elsif ($arr[3] > $threshold) { print PRR $arr[0]."\tEC6.4\n";}
			    elsif ($arr[4] > $threshold) { print PRR $arr[0]."\tEC6.5\n";}
			    elsif ($arr[5] > $threshold) { print PRR $arr[0]."\tEC6.6\n";}
			    else { print PRR $arr[0]."\tNP (Probability below thresold)\t$hash{$x}\n"; } # We have to decide wheteher this class should included or not
		    }
	  close PR; close PRR;
	  }
system ("rm *Out*"); 
}
close IF;


