#!/usr/bin/perl 
# This script was written by Josh Schmidt. Last Modified 13/2/12.
# usage: perl data_to_kaplan.pl input1 ...inputn outputfile
# input files are tab-delimited and have one header line: tray day timepoint1 timepoint2. The
# timepoints are the time elapsed since the start i.e relative to 0.
# the other lines contain the data, a single line per vial.
# their format is also tab-delimited, and contain: line_name #total_in_vial rep# #_affected_t1 #_affected_t2..... etc
# Output file. SPSS and R expect a specific data format for survival curve analysis. Essentially each individual
# has a data entry, one individual per data line: Strain Time Response	Tray_#	Day_#	Rep_#.
# the critical info is the time and response. An individual that responds i.e. dies, has the time the observation
# was first made, and a Response value of 1 (i.e it responded). Individuals that do not repsond, are still recorded, but
# given the final timepoint, and a response value of zero i.e. by the end of observation this invidual has still not responded.
# These are censored observations and are just as important for the analysis.


use strict;
my $outName = pop @ARGV; ### @ARGV is the command line input. Pop takes the last element off, which if the directions were followed is the user supplied outputfile name.
open (OUTPUT1, ">$outName") or die; ### create our outputfile.
print OUTPUT1 "strata\ttime\tstatus\treplica\tnFlies\ttimeHours\n"; ### print the header line in the output.

foreach my $fileName (@ARGV){ ##### in @ARGV, all the remaining elements are the names of input files. So foreach of these input files....
open (INPUT1, "<$fileName") or die "Couldn't open $fileName\n"; ### open that file, or report that it cant be opened.
my @inputLines = <INPUT1>; #### read all the lines into an array.
my @timePoints = split (/\t/, (shift @inputLines)); ### the first element in this new array is the header line of the input file. We get that via the shift
													### command. Splitting this on the tabs gives the @timePoints.
chomp @timePoints; ### chomp makes sure there are no cariage returns lingering.
my ($Strain, $Replicate, $Flies_total, $Time_on_Hours)= splice @timePoints, 0,4; #### cut off the first two elements in @timePoints. These should be the tray and day numbers respectively.

foreach my $dataLines (@inputLines){#### @inputLines is an array with each element being the datalines from the inputfile (rememebr we already got rid of the header line)
									#### so for each of the datalines......		
	my @lineData = split (/\t/, $dataLines); ### split on the tabs, create a new array, @lineData
	chomp @lineData;
	my ($strainNumber, $repNumber, $totalRepFlies, $hours) =  splice @lineData, 0, 4; ### cut of the first 3 elements, nd store them as the strain, toal fly, and rep numbers.
	my $numberObservations = scalar @lineData; #### how many timepoints do we have? scalar 'counts' this for us.
	my $maxObservationArrayPosition = $numberObservations -1; ### because array numbering starts at 0, the last position in the array is its length minus 1.
	my $finalTimePoint = $timePoints[$maxObservationArrayPosition]; ### create a variable $finalTimePoint, which is the last element in the @timePoints array.

### so the tricky part. we have two arrays that have our equired info in: 
### @timePoints, with our time of observation, and @lineData. They should contain the same number of elements
### which are addressed as scalars with an index number e.g $timePoints[0] is the FIRST element in the aray @timePoints.
	for (my $i=0;$i <= $maxObservationArrayPosition; $i++){#### the fun begins! The for loop takes 3 arguments. we use $i as a variable whose value changes
														#### each time we go through the loop. So the arguments are:
														#### set $i=0.
														#### keep going till $i equals $maxObservationArrayPosition
														#### $i++: increrase the value of $i by one each cycle.
															
		my $timeOfEvent = $timePoints[$i];
		
		if ($lineData[$i] == $lineData[$i-1]){;}
		elsif ($lineData[$i] > $lineData[$i-1]){
			
			my $numberResponded = $lineData[$i] - $lineData[$i-1];
			for (my $j=$numberResponded; $j>=1; $j--) {
					print OUTPUT1 "$strainNumber\t$timeOfEvent\t1\t$repNumber\t$totalRepFlies\t$hours\n";
					}													
		}
	}
	if ($totalRepFlies > $lineData[$maxObservationArrayPosition]){
				for (my $k=$totalRepFlies-$lineData[$maxObservationArrayPosition]; $k>=1; $k--){
					print OUTPUT1 "$strainNumber\t$finalTimePoint\t0\t$repNumber\t$totalRepFlies\t$hours\n";
					}
				}
}
}
close OUTPUT1;