#!/opt/bin/perl

# Author: Cyriac Kandoth
# Date: 08/28/2010
# Desc: Takes a MAF file as input and updates the first two columns (HUGO symbol and Entrez ID)

use strict;
use warnings;
use IO::File;
use LWP::Simple;
use FindBin qw($Bin); 
use lib "$Bin/lib";
use HugoGene;
use HugoGeneMethods;

if( scalar( @ARGV ) != 1 )
{
  print STDERR "\nUsage: perl $0 <maf_file_to_update>\n\n";
  exit 1;
}
my $maf_file = $ARGV[0];
my ( $total_lines ) = ( `wc -l $maf_file` =~ m/^(\d+)/ );
print STDERR $total_lines, " lines in MAF. Countdown to finish: ";
my ( $line_cnt, $progress ) = ( 0, 10 );

# Create hashes that resolve latest HUGO names
my $HugoObjRef = HugoGeneMethods::makeHugoGeneObjects();
my ( $PreviousSymbolToCurrentRef, undef ) = HugoGeneMethods::previousHugoSymbols( $HugoObjRef );
my ( $AliasToCurrentRef, undef ) = HugoGeneMethods::hugoAlias( $HugoObjRef );

my %ambigCnt = ();
my $defCnt = 0;

# Use these hashes to avoid re-processing duplicate entries in the file
my ( %hugoName, %entrezID ) = ((), ());

my $mafFh = IO::File->new( $maf_file ) or die "Can't open $maf_file. $!";
my $newFh = IO::File->new( "$maf_file\_hugo_modified", ">" ) or die "Can't open $maf_file\_hugoified. $!";
my $ambigFh = IO::File->new( "$maf_file\_ambiguous", ">" ) or die "Can't open $maf_file\_ambiguous. $!";
while( my $line = $mafFh->getline )
{
  ++$line_cnt;
  if( $line_cnt / $total_lines > ( 10 - $progress ) / 10 )
  {
    print STDERR "$progress ";
    --$progress;
  }
  if( $line =~ m/^#/ or $line =~ m/Hugo_Symbol/ )
  {
    $newFh->print( $line ); #Copy header or comment lines to new file
    next;
  }
  chomp( $line );

  my @segs = split( /\t/, $line );
  my ( $name, $chr ) = ( $segs[0], $segs[4] );

  if( defined $hugoName{"$chr:$name"} )
  {
    ++$ambigCnt{"$chr:$name"} if( $hugoName{"$chr:$name"} eq '0' );
    $segs[0] = $name = $hugoName{"$chr:$name"} unless( $hugoName{"$chr:$name"} eq '0' );
  }
  else
  {
    #Searching for it as all uppercase, makes it case insensitive
    my $hugo = getHugoSymbol( uc( $name ), $chr );
    if( defined $hugo )
    {
      $hugoName{"$chr:$name"} = $hugo;
      ++$defCnt if ( uc( $hugo ) ne uc( $name ));
      $segs[0] = $name = $hugo;
    }
    else
    {
      $hugoName{"$chr:$name"} = '0';
      $ambigFh->print( "Couldn't find a unique HUGO name for: $name\tchr$chr:$segs[5]-$segs[6]\n" );
      ++$ambigCnt{"$chr:$name"};
    }
  }

  if( defined $entrezID{$name} )
  {
    $segs[1] = $entrezID{$name};
  }
  else
  {
    #Look for the Entrez ID for that gene name and change it if necessary
    my $entrez_id = getEntrezID( $name );
    if( defined $entrez_id )
    {
      $segs[1] = $entrez_id;
    }
    else
    {
      $segs[1] = '0';
      $ambigFh->print( "Couldn't find a unique Entrez ID for: $name\t$chr:$segs[5]-$segs[6]\n" );
    }
    $entrezID{$name} = $segs[1];
  }

  #Write the modified (or unchanged) line to file
  $newFh->print( join( "\t", @segs ), "\n" );
}

$ambigFh->close;
$newFh->close;
$mafFh->close;

print "\n";
print $defCnt, " gene names needed converting to HUGO\n";
print scalar( keys %ambigCnt ), " gene names could not be resolved\n";
print "Unresolved names with more than two variants in the MAF:\n";
foreach my $key ( sort {$ambigCnt{$b} <=> $ambigCnt{$a}} ( keys %ambigCnt ))
{
  print $ambigCnt{ $key }, "\tChr$key\n" if( $ambigCnt{ $key } > 2 );
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

sub getEntrezID
{
  #These are URLs to query some Entrez utilities that we'll need
  my $esearch_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=' .
                    '<gene>[Gene%20Name]%20AND%20%22homo%20sapiens%22[Organism]';
  my $esummary_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=';

  my $gene = shift;
  $esearch_url =~ s/<gene>/$gene/;
  my ( $search_result, $result_cnt ) = ( '', '' );
  do # If the site is non-responsive, keep bugging it till we get a properly formatted page
  {
    $search_result = get( $esearch_url );
    ( $result_cnt ) = ( $search_result =~ m/<Count>(\d+)<\/Count>/ );
  } while( $result_cnt !~ m/^\d+$/ );

  if( $result_cnt == 1 )
  {
    my ( $entrez_id ) = ( $search_result =~ m/<Id>(\d+)<\/Id>/ );
    return $entrez_id;
  }
  elsif( $result_cnt > 1 )
  {
    #Resolve ambiguity by querying each ID using Entrez's E-summary utility
    my ( @entrez_ids ) = ( $search_result =~ m/<Id>(\d+)<\/Id>/g );
    foreach my $id ( @entrez_ids )
    {
      my $summary_result = get( "$esummary_url$id" );
      my ( $curr_name ) = ( $summary_result =~ /<Item Name="Name" Type="String">(\S+)<\/Item>/ );
      if( defined $curr_name && $curr_name eq $gene )
      {
        my ( $curr_id ) = ( $summary_result =~ /<Item Name="CurrentID" Type="Integer">(\d+)<\/Item>/ );
        return $curr_id if( defined $curr_id && $curr_id != 0 ); #newer ID available
        return $id if( defined $curr_id && $curr_id == 0 ); #queried ID is the current ID
      }
    }
  }

  #If the gene name is in LOC format, then the postfix is the Entrez ID
  return $1 if( $gene =~ m/^LOC(\d+)$/ );

  #Sorry comrade. In soviet Russia, Entrez ID finds YOU!
  return undef;
}
