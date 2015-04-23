#!/opt/bin/perl

# Author: Cyriac Kandoth
# Date: 02/07/2011
# Desc: Takes a FLAT file with gene names and updates them to HUGO names if possible
# Also provide the chromosome each gene is within to help resolve ambiguity


use strict;
use warnings;
use IO::File;
use LWP::Simple;
use FindBin qw($Bin); 
use lib "$Bin/lib";
use HugoGene;
use HugoGeneMethods;


if( scalar( @ARGV ) != 3 )
{
  print STDERR "\nUsage: perl $0 <tab_delimited_input_file> <gene_name_column_index> <chrom_name_column_index>";
  print STDERR "\nThe indexes must be zero-based\n\n";
  exit 1;
}

#Create hashes that resolve latest HUGO names
my $HugoObjRef = HugoGeneMethods::makeHugoGeneObjects();


my ( $PreviousSymbolToCurrentRef, undef ) = HugoGeneMethods::previousHugoSymbols( $HugoObjRef );
my ( $AliasToCurrentRef, undef ) = HugoGeneMethods::hugoAlias( $HugoObjRef );

my $flat_file = $ARGV[0];
my $gene_column = $ARGV[1];
my $chr_column = $ARGV[2];

my $annotFh = IO::File->new( $flat_file ) or die "Can't open $flat_file $!";
my $newFh = IO::File->new( "$flat_file\_hugoified", ">" ) or die "Can't open $flat_file\_hugoified $!";
my $ambigFh = IO::File->new( "$flat_file\_ambiguous", ">" ) or die "Can't open $flat_file\_ambiguous $!";
my ( %undefCnt, %defCnt, %sameCnt ) = ((), (), ());
while( my $line = $annotFh->getline )
{

  if( $line =~ m/^(#|chromosome_name|Hugo_Symbol|$)/ )
  {
    $newFh->print( $line ); #Copy any commented lines or empty lines to the new file
    next;
  }
  chomp( $line );

  my @segs = split( /\t/, $line );
  my ( $name, $chr ) = ( $segs[$gene_column], $segs[$chr_column] );

  #Searching for it as all uppercase, makes it case insensitive
  my $hugo = getHugoSymbol( uc( $name ), $chr );
  if( defined $hugo )
  {
    ++$defCnt{$hugo} if ( uc( $hugo ) ne uc( $name ));
    $segs[$gene_column] = $hugo;
  }
  else
  {
    $ambigFh->print( "Couldn't find a unique HUGO name for: $name in chrom $chr\n" );
    ++$undefCnt{$name};
  }

  #Write the modified (or unchanged) line to file
  $newFh->print( join( "\t", @segs ), "\n" );
}
$ambigFh->close;
$newFh->close;
$annotFh->close;

print scalar( keys %defCnt ), " gene names needed converting to HUGO\n";
print scalar( keys %undefCnt ), " gene names could not be resolved\n";
print "Unresolved names with more than two variants in the MAF:\n";
foreach my $key ( sort {$undefCnt{$b} <=> $undefCnt{$a}} ( keys %undefCnt ))
{
  if( $undefCnt{ $key } > 2) { print $undefCnt{ $key }, "\t$key\n"  };
}

sub getHugoSymbol
{
  my ( $oldSymbol, $chr ) = @_;
  #If it's already the latest Hugo name
  if( defined $$HugoObjRef{$oldSymbol} )
  {
    return $$HugoObjRef{$oldSymbol}->symbol();
  }

  #If it's a previous symbol of a unique Hugo name
  if( defined $$PreviousSymbolToCurrentRef{$oldSymbol} &&
         scalar( @{$$PreviousSymbolToCurrentRef{$oldSymbol}} ) == 1 &&
         $$HugoObjRef{uc( ${$$PreviousSymbolToCurrentRef{$oldSymbol}}[0] )}->sameChr( $chr ))
  {
    return ${$$PreviousSymbolToCurrentRef{$oldSymbol}}[0];
  }
  #If it's a previous symbol of multiple Hugo names, then disambiguate by chromosome
  elsif( defined $$PreviousSymbolToCurrentRef{$oldSymbol} &&
         scalar( @{$$PreviousSymbolToCurrentRef{$oldSymbol}} ) > 1 )
  {
    #Return a hugo name only if it is unique to this chromosome
    my %chrMatches = ();
    foreach my $hugo ( @{$$PreviousSymbolToCurrentRef{$oldSymbol}} )
    {
      ++$chrMatches{$hugo} if( $$HugoObjRef{uc( $hugo )}->sameChr( $chr ));
    }
    if( scalar( keys %chrMatches ) == 1 )
    {
      my @tmp = keys %chrMatches;
      return $tmp[0];
    }
  }

  #If it's an alias of a unique Hugo name
  if( defined $$AliasToCurrentRef{$oldSymbol} &&
         scalar( @{$$AliasToCurrentRef{$oldSymbol}} ) == 1 &&
         $$HugoObjRef{uc( ${$$AliasToCurrentRef{$oldSymbol}}[0] )}->sameChr( $chr ))
  {
    return ${$$AliasToCurrentRef{$oldSymbol}}[0];
  }
  #If it's an alias of multiple Hugo names, then disambiguate by chromosome
  elsif( defined $$AliasToCurrentRef{$oldSymbol} &&
         scalar( @{$$AliasToCurrentRef{$oldSymbol}} ) > 1 )
  {
    #Return a hugo name only if it is unique to this chromosome
    my %chrMatches = ();
    foreach my $hugo ( @{$$AliasToCurrentRef{$oldSymbol}} )
    {
      ++$chrMatches{$hugo} if( $$HugoObjRef{uc( $hugo )}->sameChr( $chr ));
    }
    if( scalar( keys %chrMatches ) == 1 )
    {
      my @tmp = keys %chrMatches;
      return $tmp[0];
    }
  }

  #Sorry mate. Couldn't figure it out :-|
  return undef;
}
