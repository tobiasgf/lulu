#!/usr/bin/perl
my @files = glob("*.fas");
foreach my $fn (@files) {
	$fn =~ s/.fas//; 
	print "renaming: ", $fn, "\n";
	system("sed 's/>/>$fn;/g' $fn.fas > renamed_$fn.fas");
}

