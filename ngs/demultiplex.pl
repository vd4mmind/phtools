#!/usr/bin/env perl

use strict;
use Getopt::Std;
use Text::Levenshtein qw(distance);
use FileHandle;

use vars qw($opt_h $opt_b $opt_m $opt_p $opt_o);

my $Usage="Usage:

   ... | $0 -b barcodes.txt [-m mismatches] [ -p outputprefix] [ -o outputdir ] 
";

if ( !getopts("b:p:o:m:h") || ! $opt_b ||  $opt_h ) {
    die $Usage; 
}

my  $mismatches_allowed = 1;
$mismatches_allowed = $opt_m if defined($opt_m);  # 0 also possible

die "only zero-mismatches allowed for now ..." if $mismatches_allowed >0;

my $special=0;

sub readbarcodes {
  my ($file)=@_;
  my $bc={};
  my $libs={};
  my $lengths={};
  my ($lib, $code);
  open(FILE, "$file") or die "Barcode '$file': $!";
  while(<FILE>) {
    s/[\n\r]*$//g;
    s/#.*//;
    next LINE unless $_;
    ($lib, $code)=split(' ');            # e.g. 'G7 \t CCAACAAT'
    die "Library '$lib' not unique" if $libs->{$lib}++;
    die "Barcode '$code' not unique" if $bc->{$code};
    $lengths->{length($code)}++;        # should be all same length
    $bc->{$code}=$lib;
  }
  close(FILE);
  my (@lengths)=keys %$lengths;
  if ( int(@lengths)  > 1 ) {
    die "lengths of the barcodes not all the same, found:" . 
        join(',', @lengths) . "\n";
  }
  ## also add the mononucleotides under their own sequence as name
  for my $nuc (qw(A C G T N)) {
    die "adding of mononucs untested";
    $code=$nuc x $lengths[0];
    $bc->{$code}=$code;
  }
  $bc;
}                                       # readbarcodes

sub open_infile {
  die "not used nor tested";
  my($file)=@_;
  my $fh=FileHandle->new();
  if ($file =~ /\.gz/) { 
    $fh->open("zcat $file | ", "r")  or die "'$file': $!";
  } else { 
    $fh->open("< $file")  or die "'$file': $!";
  }
  $fh;
}

sub open_outfiles { 
  my(@libs)=@_;
  my $fhs={};
  my $names={};

  for my $lib (@libs) { 
    my $name=sprintf("%s.fastq.gz", $lib);
    $name="$opt_p$name" if $opt_p;
    $name="$opt_o/$name" if $opt_o;
    my $fh = FileHandle->new("| gzip > $name") or die "library $lib, file $name: $!";
    warn "Creating/overwriting file $name ...\n";
    $fhs->{$lib}=$fh;
    $names->{$lib}=$name;
  }
  ($fhs, $names);
}                                       # open_outfiles

sub close_outfiles {
  my($fhs, $names)=@_;
  
  for my $lib (keys %$fhs) {
    my $name=$names->{$lib};
    $fhs->{$lib}->close() or die "could not close demultiplexed file $name for library $lib; investigate";
  }
}

sub rescue {
  die "to be written";
}

sub ambiguous {
  die "to be written";
}

my $codes = readbarcodes($opt_b);       # eg. $code->{'AGCGTT') => 'M3'
my @files=(values %$codes, 'AMBIGUOUS', 'UNKNOWN');
my ($filehandles, $filenames)=open_outfiles(@files);
my $filecounts={};
my $nexact=0;
my $nrescued=0;                         # having at most $mismatch mismatches
my $nunknown=0;
my $nambiguous=0;

RECORD:
while(1) { 
  my $record=<>;
  last RECORD if (eof(STDIN) || !$record);
  ### e.g.:  ^@NS500413:172:HVFHWBGXX:1:11101:4639:1062 1:N:0:CCGTCCAT$
  my ($code)=(split(':', $record))[-1];
  $code =~ s/[\n\r]*$//;
  $record .= <>; # sequence line
  $record .= <>; # '+'
  $record .= <>; # quality line
  
  my $lib;
 CASE:
  while(1) {
    $lib=$codes->{$code};
    if ($lib) {
      $nexact++;
      last CASE;
    }
    if ($mismatches_allowed == 0) {
      $nunknown++;
      $lib='UNKNOWN';
      last CASE;
    }
    $lib=rescue($code, $codes, $mismatches_allowed);
    if(!$lib) {
      $nunknown++;
      $lib='UNKNOWN';
      last CASE;
    }
    if($special) {
      # check if mismatch is in 7th bp; if so, call it ambiguous
      if (ambiguous() ) {
        $lib='AMBIGUOUS';
        $nambiguous++;
        last CASE;
      } else {
        $nrescued++;
        last CASE;
      }
    }
    die "should not reach this point";
  }                                     # CASE
  $filehandles->{$lib}->print($record);
  $filecounts->{$lib}++;
}                                       # RECORD
close_outfiles($filehandles);

for my $lib (keys %$filecounts) {
  die "closing of empty files untested";
  unlink $filenames->{$lib} unless $filecounts->{$lib};
}

sub commafy {
  # insert comma's to separate powers of 1000
  my($i)=@_;
  my $r = join('',reverse(split('',$i)));
  $r =~ s/(\d{3})/$1,/g;
  $r =~ s/,$//;
  join('',reverse(split('',$r)));
}

warn sprintf("exact: %s\nrescued:%s\nambiguous:%s\nunknown: %s\n",
             map { commafy $_ } ($nexact, $nrescued, $nambiguous, $nunknown ));
