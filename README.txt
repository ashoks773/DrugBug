Download the complete program from
http://metagenomics.iiserb.ac.in/drugbug/download.php

Please note that this stand alone version of the tool is for the testing purpose and this will be assessible only to reviwers. The distribution of complete stand alone package for the users with all the suporting custom made databases and tools will be available soon once the manuscript is accepted. 
NOTE: Only input.sdf provided with this package can be used as a inputfile. (Only sample database for this molecule provided for reviewers to check the functionality)
Program requires OpenBabel as a prerequisite please make sure that it is accessible from working directory. For downloading and installing OpenBabel package please visit http://sourceforge.net/projects/openbabel/

For exceution run the following command on terminal. 
$tar -zxvf SampleRun.tar.gz 
$cd SampleRun 
$./run.exe input.sdf

Or to directly run the perl program use following command
perl XenoSampleRun.pl input.sdf

Please check "inputfile".Genus.Final.Result file for taxonomy results and "inputfile".proteins.Final.Result file for complete results 