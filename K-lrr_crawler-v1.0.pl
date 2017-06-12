#!/bin/perl
 
use strict;
use Getopt::Long;
use Data::Dumper;
use Bio::SeqIO;

my ($file_pep, $file_domains, $out_file, $motif);
my $before='';
my $after='';

GetOptions(
       'p|inseq_pep:s'        => \$file_pep,
       'd|inseq_domains:s'        => \$file_domains,
       'o|outfile:s'        => \$out_file,
       'm|motif:s'        => \$motif,
       'a|after:s'   => \$after,
       'b|before:s'   => \$before,

    );

my $usage="                                                                                                                                                                                                                                                                     
Usage: perl script.pl <options>                                                                                                                               Required options:                                                                                                                                                                                                                                                               
       '-p|--inseq_pep'  input fasta protein file
       '-d|--inseq_domains'   input domain annotation
       '-o|--outfile'    prefix for all output files
       '-m|--motif'      degenerate motif to search for
       '-b|--before'     consider motives only N' to this domain
       '-a|--after'      consider motives only C' to this domain                        
";

die $usage unless ( (defined $file_pep) and (defined $motif) and (defined $file_domains) );

# Load all fasta sequences

my %sequences;

my $inseq_pep    = Bio::SeqIO->new(
									-file => $file_pep,
									-format => 'Fasta',
									);

while (my $seq = $inseq_pep->next_seq){

	my ($header)=split(" ", $seq->id);
	$sequences{$header}{'pep'}=$seq->seq();

}

# Load domain annotations

#print $file_domains, "\n";
my %domains = load_domain_annotation($file_domains);

#define aa palette 
my $palette = {};
@$palette{'L','I','W','F','Y','V','P'}=('yellow') x 7;
@$palette{'E','D'}=('red') x 2;
@$palette{'K','R','H'}=('blue') x 3;
@$palette{'S','T','Q','N','C'} = ('green') x 5;
@$palette{'A','G'} =('grey') x 2;
@$palette{'M','*','X'} =('black') x 3;
#search for motif in each sequence and record aa composition

my %aa_composition;
my %obs_total;

open FILEOUT1, ">", $out_file . ".motifs.tsv";

foreach my $seq_id (keys %sequences){

    while ($sequences{$seq_id}{'pep'}=~ m/$motif/g ){
	
		if (defined($after)){
	
			my ($foo, $scan_after) = split ("-", $domains{$seq_id}{$after}[-1]);

			if 	($-[0] > $scan_after){
	
				print FILEOUT1 $seq_id, "\t", $-[0], "\t", substr($sequences{$seq_id}{'pep'}, $-[0], 8), "\n";	
				
				my $n=0;
				
				while ($n<8){
				
				my $aa = substr($sequences{$seq_id}{'pep'}, $-[0]+$n, 1);
				$aa_composition{$n}{$aa}++;
				$obs_total{$n}++;
				$n++;
				}
			}	
	    }
    }
}
close FILEOUT1;

open FILEOUT2, ">", $out_file . ".aa_composition.tsv";

foreach my $pos (sort keys %aa_composition){

	foreach my $aa (sort keys %{$aa_composition{$pos}}){
	
		print FILEOUT2 $pos, "\t", $aa, "\t", $aa_composition{$pos}{$aa}, "\t", 
		$aa_composition{$pos}{$aa}/$obs_total{$pos}, "\t", $palette->{$aa}, "\n";
		
	}

}

close FILEOUT2;

sub load_domain_annotation {
    
    my ($infile) = @_[0];
  
    open (FILE1, "<", $infile) or die "cannot open domain annotation file";
    
    my %domain_db;

    while (my $line = <FILE1> ){
    
	chomp $line; 

	my ($seqid, $domainstring)=split("\t", $line);
	    
	my @domains = split('~', $domainstring);

	foreach (@domains){

	    my $domain=$_;

	    my ($domainid, $attributes)=split(/\(/, $domain);
	    $attributes=~ s/\)//; #remove right hand bracket
	    $attributes=~ s/\s//g;
	    my ($start_string, $stop_string, $evalue_string)=split(",", $attributes);
		my ($foo, $start)=split("=", $start_string);
		my ($foo, $stop)=split("=", $stop_string);
	    my ($foo, $evalue)=split("=", $evalue_string);

#	    if ($evalue < $evalue_cutoff){
		my $record = $start . "-" . $stop;
	    push @{$domain_db{$seqid}{$domainid}}, $record;
#	}
	}
    
    }

    close FILE1;

    return %domain_db;

}


sub match_all_positions {
    my ($regex, $string) = @_;
    my @ret;
    while ($string =~ /$regex/g) {
        push @ret, [ $-[0], $+[0] ];
    }
    return @ret;
}

sub count_matches{

 my ($regex, $string) = @_;
 my @count = $string =~ /$regex/g;
 
 return scalar @count; 

}
