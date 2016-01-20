#!/opt/common/CentOS_6/perl/perl-5.20.2/bin/perl

#use lib qw(/ifs/data/zeng/dmp/resources/perllib);
use strict;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);
use Tie::IxHash;
use File::HomeDir;

my ($pre, $input_dir, $output_dir, $PERL, $hg19);

GetOptions ('pre=s' => \$pre,
		'input=s' => \$input_dir,
		'output=s' => \$output_dir,
		'hg19' => \$hg19,
            	'perl=s' => \$PERL) or exit(1);


if(!$pre || !$input_dir || !$output_dir || !$PERL){
    print <<HELP;

    USAGE: ./dmp_pipeline_wrapper_post_processing.pl [options]
    * -pre <string>: project prefix (project name)
    * -input <string>: input diretory
    * -output <string>: output directory
    * -perl <string>: path to perl executable

HELP
exit;
}

die "[ERROR]: perl executable does not exist: $PERL" if(!-e $PERL);
die "[ERROR]: input directory does not exist: $input_dir" if(!-d $input_dir);
die "[ERROR]: output directory does not exist: $output_dir" if(!-d $output_dir);



MergeCopyNumberSeg($input_dir, "$output_dir/${pre}_ALL_copynumber.seg");

### generate discrete copy number calls
GenerateCopyNumberCalls("$input_dir/${pre}_copynumber_segclusp.genes.txt", "$output_dir/${pre}_discrete_CNA.txt");

### copy over files to generate report for PI
GenerateReport();

### merge all copy number seg file, replace . with _ in header, and replace space with tab in the file
sub MergeCopyNumberSeg {
	my ($seg_dir, $seg_merged_file) = @_;
    print "[INFO]: Merging seg files in $seg_dir to $seg_merged_file\n";
	opendir(workDir, "$seg_dir");
    my @unsorted = readdir workDir;
    closedir workDir;
    my @files = sort @unsorted;

	open(SM, ">$seg_merged_file") or die "[ERROR]: Can't write to $seg_merged_file $!";
	my $print_header = 1;
    foreach my $file (@files){
            if($file =~ m/^.*copynumber\.seg$/){
            	print "[INFO]: Processing $file\n";
            	open(CS, "$seg_dir/$file") or die "[ERROR]: Can't open seg file $file $!";
            	my $cs_header = <CS>;
            	if($print_header)
            	{
        			chomp($cs_header);
        			#$cs_header =~ s/\./_/g;
        			$cs_header =~ s/\s+/\t/g;
            		print SM "$cs_header\n";
            		$print_header = 0;
            	}
 				while(<CS>)
 				{
                	my $data_line = $_;
                	chomp($data_line);
                	$data_line =~ s/\s+/\t/g;
			if($hg19)
			{
				my @data = split(/\t/, $data_line);
				$data[1] = "chr" . $data[1];
				$data_line = join("\t", @data);
			}
                	print SM "$data_line\n";      
        		}
            	close CS;
            }
        }
	close SM;
}



sub GenerateCopyNumberCalls{
	my ($input_file, $output_file) = @_;
	`$PERL $Bin/generate_copy_number_calls.pl $input_file >$output_file`;
}




### create directory and copy over files that need to be given to PI
sub GenerateReport
{
	`cp $input_dir/${pre}_copynumber_segclusp.pdf $output_dir/`;
	if(-e "$input_dir/${pre}_copynumber_segclusp.nvn.pdf")
	{
		`cp $input_dir/${pre}_copynumber_segclusp.nvn.pdf $output_dir/${pre}_normvsnorm_segclusp.pdf`;
	}
}



