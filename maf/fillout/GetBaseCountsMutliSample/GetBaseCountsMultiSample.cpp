//
//  GetBaseCountsMultiSample.cpp
//  GetBaseCountsMultiSample
//
//  Created by Zeng, Zheng/Sloan-Kettering Institute on 9/5/14.
//  Copyright (c) 2014 Zeng, Zheng/Sloan-Kettering Institute. All rights reserved.
//



#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sstream>
#include <fstream>
#include <getopt.h>
#include <vector>
#include <algorithm>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include "api/BamReader.h"
#include "omp.h"

//#define _DEBUG
//#define _PARSING_DEBUG
//#define _FASTA_DEBUG

using namespace std;
using namespace BamTools;


const string VERSION = "GetBaseCountsMultiSample 1.1.1";

string input_fasta_file;
map<string, string> input_bam_files;
vector<string> input_variant_files;
string output_file;

int mapping_quality_threshold = 20;
int base_quality_threshold = 20;
int quality_scale = 33;
int filter_duplicate = 1;
int filter_improper_pair = 1;
int filter_qc_failed = 0;
int filter_indel = 0;
int filter_non_primary = 1;
char bases[] = {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't', 'R', 'r', 'M', 'm', 'U', 'u'}; // M/m for matched insertion/deletion, u/U for unmatched insertion/deletion, R/r for matched reference
size_t base_size = 14;
int maximum_variant_block_size = 10000;
int maximum_variant_block_distance = 100000;
int num_thread = 1;
string count_method = "DMP";
bool input_variant_is_maf = false;
bool input_variant_is_vcf = false;
const size_t BIN_SIZE = 16*1024;


void printUsage(string msg = "")
{
    cout << endl;
    cout << VERSION << endl;
    cout << "Usage: " << endl;
    cout << "[REQUIRED ARGUMENTS]" << endl;
    cout << "\t--fasta                 <string>                        Input reference sequence file" << endl;
    cout << "\t--bam                   <string>                        Input bam file, in the format of SAMPLE_NAME:BAM_FILE. This paramter need to be specified at least once" << endl;
    cout << "\t                                                        e.g: --bam s_EV_crc_007:Proj_4495.eta_indelRealigned_recal_s_EV_crc_007_M3.bam." << endl;
    cout << "\t--maf                   <string>                        Input variant file in TCGA maf format. --maf or --vcf need to be specified at least once. But --maf and --vcf are mutually exclusive" << endl;
    cout << "\t--vcf                   <string>                        Input variant file in vcf-like format(the first 5 columns are used). --maf or --vcf need to be specified at least once. But --maf and --vcf are mutually exclusive" << endl;
    cout << "\t--output                <string>                        Output file" << endl;
    cout << endl;
    cout << "[OPTIONAL ARGUMENTS]" << endl;
    cout << "\t--thread                <int>                           Number of thread. Default " << num_thread << endl;
    cout << "\t--maq                   <int>                           Mapping quality threshold. Default " << mapping_quality_threshold << endl;
    cout << "\t--baq                   <int>                           Base quality threshold, Default " << base_quality_threshold << endl;
    cout << "\t--filter_duplicate      [0, 1]                          Whether to filter reads that are marked as duplicate. 0=off, 1=on. Default 1" << endl;
    cout << "\t--filter_improper_pair  [0, 1]                          Whether to filter reads that are marked as improperly paired. 0=off, 1=on. Default 1" << endl;
    cout << "\t--filter_qc_failed      [0, 1]                          Whether to filter reads that are marked as failed quality control. 0=off, 1=on. Default 0" << endl;
    cout << "\t--filter_indel          [0, 1]                          Whether to filter reads that have indels. 0=off, 1=on. Default 0" << endl;
    cout << "\t--filter_non_primary    [0, 1]                          Whether to filter reads that are marked as non primary alignment. Default 1" << endl;
    cout << "\t--help                                                  Print command line usage" << endl;
    cout << endl;
    cout << "[ADVANCED ARGUMENTS, CHANGING THESE ARGUMENTS WILL SIGNIFICANTLY AFFECT MEMORY USAGE AND RUNNING TIME. USE WITH CAUTION]" << endl;
    cout << "\t--max_block_size        <int>                           The maximum size of vcf chunks that can be processed at once per thread. Default " << maximum_variant_block_size << endl;
    cout << "\t--max_block_dist        <int>                           The longest spanning region (bp) of vcf chunks that can be processed at once per thread. Default " << maximum_variant_block_distance << endl;
    cout << endl;
    if(!msg.empty())
        cerr << msg << endl;
    exit(1);
}


static struct option long_options[] =
{
    {"fasta",                   required_argument,      0,     'f'},
    {"bam",                     required_argument,      0,     'b'},
    {"maf",                     required_argument,      0,     'v'},
    {"vcf",                     required_argument,      0,     'V'},
    {"output",                  required_argument,      0,     'o'},
    {"thread",                  required_argument,      0,     't'},
    {"maq",                     required_argument,      0,     'Q'},
    {"baq",                     required_argument,      0,     'q'},
    {"filter_duplicate",        required_argument,      0,     'd'},
    {"filter_improper_pair",    required_argument,      0,     'p'},
    {"filter_qc_failed",        required_argument,      0,     'l'},
    {"filter_indel",            required_argument,      0,     'i'},
    {"filter_non_primary",      required_argument,      0,     'n'},
    {"max_block_size",          required_argument,      0,     'M'},
    {"max_block_dist",          required_argument,      0,     'm'},
    {"help",                    no_argument,            0,     'h'},
    {0, 0, 0, 0}
};

void split(const string& line, char delim, vector<std::string>& parsed_item, bool ignore_empty_item = false) // split a string with specified delimiter
{
    stringstream my_ss(line);
    string item;
    while (getline(my_ss, item, delim))
    {
#ifdef _PARSING_DEBUG
        cout << "[DEBUG] Parsed item: " << item << endl;
#endif
        if (ignore_empty_item && item.empty())
            continue;    // use to skip empty item
        parsed_item.push_back(item);
    }
}

void addBamFile(string bam_string)  // check and add a bam file to input list
{
    vector<string> bam_items;
    split(bam_string, ':', bam_items, true);
    if(bam_items.size() != 2)
    {
        cerr << "[ERROR] Incorrect format of -bam parameter: " << bam_string << endl;
        exit(1);
    }
    string sample_name = bam_items[0];
    string bam_file_name = bam_items[1];
    if(input_bam_files.find(sample_name) != input_bam_files.end())
    {
        cerr << "[ERROR] Multiple bam files specified for sample: " << sample_name << endl;
        exit(1);
    }
    struct stat buffer;
    if(stat(bam_file_name.c_str(), &buffer) != 0)
    {
        cerr << "[ERROR] Unable to access bam file:" << bam_file_name << endl;
        exit(1);
    }
    input_bam_files.insert(make_pair(sample_name, bam_file_name));
}

void addVariantFile(string variant_string)  // add a variant file to input list
{
    input_variant_files.push_back(variant_string);
}


void parseOption(int argc, const char* argv[])
{
    if(argc == 1)
        printUsage();
    int next_option;
    int option_index = 0;
    do
    {
        next_option = getopt_long(argc, const_cast<char**>(argv), "f:b:v:V:o:t:Q:q:d:p:l:i:n:M:m:h", long_options, &option_index);
        switch(next_option)
        {
            case 'f':
                input_fasta_file = optarg;
                break;
            case 'b':
                addBamFile(optarg);
                break;
            case 'v':
                addVariantFile(optarg);
                input_variant_is_maf = true;
                break;
            case 'V':
                addVariantFile(optarg);
                input_variant_is_vcf = true;
            case 'o':
                output_file = optarg;
                break;
            case 't':
                num_thread = atoi(optarg);
                break;
            case 'Q':
                mapping_quality_threshold = atoi(optarg);
                break;
            case 'q':
                base_quality_threshold = atoi(optarg);
                break;
            case 'd':
                filter_duplicate = atoi(optarg);
                break;
            case 'p':
                filter_improper_pair = atoi(optarg);
                break;
            case 'l':
                filter_qc_failed = atoi(optarg);
                break;
            case 'i':
                filter_indel = atoi(optarg);
                break;
            case 'n':
                filter_non_primary = atoi(optarg);
                break;
            case 'M':
                maximum_variant_block_size = atoi(optarg);
                break;
            case 'm':
                maximum_variant_block_distance = atoi(optarg);
                break;
            case 'h':
                printUsage();
                break;
            case -1:
                break; //parsed all options
            default:
                printUsage(string("[ERROR] Argument error: ") + argv[optind - 1]);
        }
    } while(next_option != -1);
    if(input_fasta_file.empty())
        printUsage("[ERROR] Please specify input fasta file");
    if(input_bam_files.empty())
        printUsage("[ERROR] Please specify at least one input bam file");
    if(input_variant_files.empty())
        printUsage("[ERROR] Please specify at least one input variant file with --maf or --vcf");
    if(output_file.empty())
        printUsage("[ERROR] Please specify output file");
    if(num_thread <= 0)
        printUsage("[ERROR] Invalid number of threads");
    if(maximum_variant_block_size <= 0)
        printUsage("[ERROR] Invalid max_block_size");
    if(maximum_variant_block_distance <= 0)
        printUsage("[ERROR] Invalid max_block_dist");
    if(filter_duplicate != 0 && filter_duplicate != 1)
        printUsage("[ERROR] --filter_duplicate should be 0 or 1");
    if(filter_improper_pair != 0 && filter_improper_pair != 1)
        printUsage("[ERROR] --filter_improper_pair should be 0 or 1");
    if(filter_qc_failed != 0 && filter_qc_failed != 1)
        printUsage("[ERROR] --filter_qc_failed should be 0 or 1");
    if(filter_indel != 0 && filter_indel != 1)
        printUsage("[ERROR] --filter_indel should be 0 or 1");
    if(filter_non_primary != 0 && filter_non_primary != 1)
        printUsage("[ERROR] --filter_non_primary should be 0 or 1");
    if(input_variant_is_maf && input_variant_is_vcf)
        printUsage("[ERROR] --maf and --vcf are mutually exclusive");
    base_quality_threshold += quality_scale;
#ifdef _DEBUG
    cout << "[DEBUG] Parsing options complete." << endl;
#endif
}



void outputReferenceSequence(string output_fastafile, map<string, string>& ref_seq, vector<string>& orignal_header)
{
    size_t base_per_line = 80;
    ofstream out_fs(output_fastafile.c_str());
    for(size_t i = 0; i < orignal_header.size(); i++)
    {
        map<string, string>::iterator it = ref_seq.find(orignal_header[i]);
        size_t chrom_len = it->second.length();
#ifdef _FASTA_DEBUG
        cout << "[DEBUG] output reference sequence: " << it->first << ": " << chrom_len << endl;
#endif
        size_t current_len = 0;
        out_fs << ">" << it->first << endl;
        while(current_len < chrom_len)
        {
            size_t output_len = min(base_per_line, chrom_len - current_len);
            out_fs << it->second.substr(current_len, output_len) << endl;
            current_len += output_len;
        }
    }
    out_fs.close();
}


void loadReferenceSequenceSpeedup(string fasta_filename, map<string, string>& ref_seq)   // load refseq fasta
{
    cout << "[INFO] Loading reference sequence: " << fasta_filename << endl;
    string fasta_index_filename = fasta_filename + ".fai";
    ifstream index_fs(fasta_index_filename.c_str());
    if(!index_fs)
    {
        cerr << "[ERROR] fail to open reference fasta index file: " << fasta_index_filename << endl;
        exit(1);
    }
#pragma omp parallel num_threads(num_thread)
    {
        int thread_num = omp_get_thread_num();
	    string line;
	    ifstream ref_fs(fasta_filename.c_str());
	    if(!ref_fs)
	    {
	        cerr << "[ERROR] fail to open reference fasta file: " << fasta_filename << endl;
	        exit(1);
	    }
	    while(!index_fs.eof())
	    {
#pragma omp critical(read_index)
	        {
	        	getline(index_fs, line);
            }
            if(!index_fs.eof())
            {
		        vector<string> index_items;
		        split(line, '\t', index_items);
		        string chrom_name = index_items[0];
		        int chrom_len = atoi(index_items[1].c_str());
		        long long int chrom_offset = atoll(index_items[2].c_str());
		        int base_per_line = atoi(index_items[3].c_str());
		        int byte_per_line = atoi(index_items[4].c_str());
		        int byte_len = chrom_len + (chrom_len / base_per_line) * (byte_per_line - base_per_line);
		        string new_seq;
#pragma omp critical(access_ref_seq_map)
		        {
		        	ref_seq.insert(make_pair(chrom_name, new_seq));
		        	ref_seq[chrom_name].resize(chrom_len);
		      	}
		        ref_fs.seekg(chrom_offset);
		        char* seq_buff = new char[byte_len];
		        ref_fs.read(seq_buff, byte_len);
                string::iterator it_target;
#pragma omp critical(access_ref_seq_map)
		        {
                    it_target = ref_seq[chrom_name].begin();
                }
		        char* it_source = seq_buff;
		        for(int i = 0; i < byte_len; i++)
		        {
		            if(!isspace(*it_source))
		            {
		                *it_target = toupper(*it_source);
		                it_target ++;
		            }
		            it_source ++;
		        }
		        delete[] seq_buff;
	      	}
	    }
	    ref_fs.close();
  	}
    index_fs.close();
    cout << "[INFO] Finished loading reference sequence" << endl;
}

class VariantFile
{
public:
	
	VariantFile(string input_variant_file)
	{
        variant_fs.open(input_variant_file.c_str());
        if(!variant_fs)
        {
            cerr << "[ERROR] fail to open input variant file: " << input_variant_file << endl;
            exit(1);
        }
	}
    
	~VariantFile() {}
	
	bool get_next(string &line)
	{
		//if(cur_line.empty())
		//{
        if(!getline(variant_fs, line))
            return false;
		//}
		//else
		//{
		//	line = cur_line;
		//	cur_line.clear();
		//}
		return true;
	}
	
	//void roll_back(string &line)
	//{
	//	cur_line = line;
	//}
	
	void close()
	{
		if(variant_fs)
			variant_fs.close();
	}
    
	bool eof()
	{
		return variant_fs.eof();
	}
	
	ifstream variant_fs;
	//string cur_line;
};


class VariantEntry
{
public:
    
    VariantEntry()
    {
        for(map<string, string>::iterator it = input_bam_files.begin(); it != input_bam_files.end(); it ++)
        {
            for(size_t i = 0; i < base_size; i++)
            {
                base_count[it->first][bases[i]] = 0;
            }
        }
        pos = 0;
        end_pos = 0;
        insertion = false;
        deletion = false;
        t_ref_count = 0;
        t_alt_count = 0;
        n_ref_count = 0;
        n_alt_count = 0;
        duplicate_variant_ptr = NULL;
    }
    
    ~VariantEntry() {}
    
    string chrom;
    int pos;
    int end_pos;
    string ref;
    string alt;
    string tumor_sample;
    string normal_sample;
    bool insertion;
    bool deletion;
    map<string, map<char, int> > base_count;    //<sample_name, <allele, count> >
    string gene;                                // gene name from input maf
    string effect;                              // effect from input maf
    int t_ref_count;                            // tumor ref count from input maf
    int t_alt_count;                            // tumor alt count from input maf
    int n_ref_count;                            // normal ref count from input maf
    int n_alt_count;                            // normal alt count from input maf
    VariantEntry* duplicate_variant_ptr;        // pointer to the first identical variant, no need to do duplicate count for the same variant of different normal/tumor pair
};

bool isNumber(const string& s)  //check if a string(chrom name) is number
{
    string::const_iterator it = s.begin();
    while (it != s.end() && isdigit(*it))
    {
        it++;
    }
    return !s.empty() && it == s.end();
}

bool sortVariantByPos(const VariantEntry& lhs, const VariantEntry& rhs) // sort varaints by genomic position, 1-22, X, Y, etc
{
    if(lhs.chrom != rhs.chrom)
    {
        string chrom1 = lhs.chrom;
        string chrom2 = rhs.chrom;
        if(chrom1.length() > 3 && chrom1.substr(0, 3) == "chr")
            chrom1 = chrom1.substr(3);
        if(chrom2.length() > 3 && chrom2.substr(0, 3) == "chr")
            chrom2 = chrom2.substr(3);
        bool chrom1_num = isNumber(chrom1);
        bool chrom2_num = isNumber(chrom2);
        if(chrom1_num && !chrom2_num)
            return true;
        if(!chrom1_num && chrom2_num)
            return false;
        if(chrom1_num && chrom2_num)
            return atoi(chrom1.c_str()) < atoi(chrom2.c_str());
        else //both string
            return chrom1 < chrom2;
    }
    else
        return lhs.pos < rhs.pos;
}


bool alignmentHasIndel(BamAlignment& my_bam_alignment) // check if an alignment has indel
{
    vector<CigarOp>::const_iterator cigarIter;
    for(cigarIter = my_bam_alignment.CigarData.begin() ; cigarIter != my_bam_alignment.CigarData.end(); cigarIter++)
    {
        if(cigarIter->Type == 'I' || cigarIter->Type == 'D')
            return true;
    }
    return false;
}


void loadVariantFileDMPTXT(vector<string>& input_file_names, vector<VariantEntry>& variant_vec, vector<pair<size_t, size_t> >& variant_block_vec)  // load variant from dmp-impact annotated_exonic_variants.txt
{
    for(size_t i = 0; i < input_file_names.size(); i++)
    {
        int variant_count = 0;
        VariantFile my_vf(input_file_names[i]);
        string line;
        my_vf.get_next(line); //header
        while(my_vf.get_next(line))
        {
            vector<string> variant_items;
            split(line, '\t', variant_items);
            if(variant_items.size() < 6)
            {
                cerr << "[ERROR] Incorrectly formatted input variant entry:" << endl;
                cerr << "[ERROR] " << line << endl;
                exit(1);
            }
            VariantEntry new_variant;
            new_variant.tumor_sample = variant_items[0];
            new_variant.normal_sample = variant_items[1];
            new_variant.chrom = variant_items[2];
            new_variant.pos = atoi(variant_items[3].c_str()) - 1; //convert to 0-indexed, to be consistent with bam entry
            new_variant.ref = variant_items[4];
            new_variant.alt = variant_items[5];
            new_variant.insertion = (new_variant.alt.length() > new_variant.ref.length());
            new_variant.deletion = (new_variant.alt.length() < new_variant.ref.length());
            variant_vec.push_back(new_variant);
            variant_count ++;
        }
        my_vf.close();
        cout << "[INFO] " << variant_count << " variants has been loaded from file: " << input_file_names[i] << endl;
    }
}


void loadVariantFileVCF(vector<string>& input_file_names, vector<VariantEntry>& variant_vec)  // load variant from vcf-like format, only the first 5 columns are used
{
    for(size_t i = 0; i < input_file_names.size(); i++)
    {
        int variant_count = 0;
        VariantFile my_vf(input_file_names[i]);
        string line;
        while(my_vf.get_next(line))
        {
            if(line[0] == '#')
                continue; // header
            vector<string> variant_items;
            split(line, '\t', variant_items);
            if(variant_items.size() < 5)
            {
                cerr << "[ERROR] Incorrectly formatted input variant entry:" << endl;
                cerr << "[ERROR] " << line << endl;
                exit(1);
            }
            VariantEntry new_variant;
            new_variant.chrom = variant_items[0];
            new_variant.pos = atoi(variant_items[1].c_str()) - 1; //convert to 0-indexed, to be consistent with bam entry
            new_variant.ref = variant_items[3];
            new_variant.alt = variant_items[4];
            new_variant.insertion = (new_variant.alt.length() > new_variant.ref.length());
            new_variant.deletion = (new_variant.alt.length() < new_variant.ref.length());
            variant_vec.push_back(new_variant);
            variant_count ++;
        }
        my_vf.close();
        cout << "[INFO] " << variant_count << " variants has been loaded from file: " << input_file_names[i] << endl;
    }
}


void loadVariantFileMAF(vector<string>& input_file_names, vector<VariantEntry>& variant_vec, map<string, string>& reference_sequence)  // load from exome-pipeline tcga maf format
{
    for(size_t file_index = 0; file_index < input_file_names.size(); file_index++)
    {
        int variant_count = 0;
        map<string, size_t> header_index_hash;
        
        VariantFile my_vf(input_file_names[file_index]);
        string line;
        
        my_vf.get_next(line); //header
        vector<string> header_items;
        split(line, '\t', header_items);
        for(size_t header_index = 0; header_index < header_items.size(); header_index++)
        {
            header_index_hash[header_items[header_index]] = header_index;
        }
        if(header_index_hash.find("Hugo_Symbol") == header_index_hash.end() || header_index_hash.find("Chromosome") == header_index_hash.end() ||
           header_index_hash.find("Start_Position") == header_index_hash.end() || header_index_hash.find("End_Position") == header_index_hash.end() ||
           header_index_hash.find("Reference_Allele") == header_index_hash.end() || header_index_hash.find("Tumor_Seq_Allele2") == header_index_hash.end() ||
           header_index_hash.find("Tumor_Sample_Barcode") == header_index_hash.end() || header_index_hash.find("Matched_Norm_Sample_Barcode") == header_index_hash.end() ||
           header_index_hash.find("t_ref_count") == header_index_hash.end() || header_index_hash.find("t_alt_count") == header_index_hash.end() ||
           header_index_hash.find("n_ref_count") == header_index_hash.end() || header_index_hash.find("n_alt_count") == header_index_hash.end() ||
           header_index_hash.find("Variant_Classification") == header_index_hash.end())
        {
            cerr << "[ERROR] Incorrect TCGA MAF file header:" << endl;
            cerr << "[ERROR] " << line << endl;
            exit(1);
        }
        
        while(my_vf.get_next(line))
        {
            vector<string> variant_items;
            split(line, '\t', variant_items);
            VariantEntry new_variant;
            new_variant.gene = variant_items[header_index_hash["Hugo_Symbol"]];
            new_variant.chrom = variant_items[header_index_hash["Chromosome"]];
            new_variant.pos = atoi(variant_items[header_index_hash["Start_Position"]].c_str()) - 1; //convert to 0-indexed, to be consistent with bam entry
            new_variant.end_pos = atoi(variant_items[header_index_hash["End_Position"]].c_str()) - 1;
            new_variant.ref = variant_items[header_index_hash["Reference_Allele"]];
            new_variant.alt = variant_items[header_index_hash["Tumor_Seq_Allele2"]];
            new_variant.tumor_sample = variant_items[header_index_hash["Tumor_Sample_Barcode"]];
            new_variant.normal_sample = variant_items[header_index_hash["Matched_Norm_Sample_Barcode"]];
            new_variant.t_ref_count = atoi(variant_items[header_index_hash["t_ref_count"]].c_str());
            new_variant.t_alt_count = atoi(variant_items[header_index_hash["t_alt_count"]].c_str());
            new_variant.n_ref_count = atoi(variant_items[header_index_hash["n_ref_count"]].c_str());
            new_variant.n_alt_count = atoi(variant_items[header_index_hash["n_alt_count"]].c_str());
            new_variant.effect = variant_items[header_index_hash["Variant_Classification"]];
            if(new_variant.effect.empty())
            {
                if(header_index_hash.find("ONCOTATOR_VARIANT_CLASSIFICATION") != header_index_hash.end() && !(variant_items[header_index_hash["ONCOTATOR_VARIANT_CLASSIFICATION"]].empty()) )
                {
                    new_variant.effect = variant_items[header_index_hash["ONCOTATOR_VARIANT_CLASSIFICATION"]];
                }
                else
                {
                    cerr << "[ERROR] Could not find annotation information for variant:" << endl;
                    cerr << "[ERROR] " << line << endl;
                    exit(1);
                }
            }
            if(new_variant.ref == "-") // convert maf convention insertion to vcf convention         - A => G GA
            {
                char prev_ref = reference_sequence[new_variant.chrom][new_variant.pos]; // get allele before the current position
                new_variant.ref = prev_ref;
                new_variant.alt = prev_ref + new_variant.alt;
            }
            if(new_variant.alt == "-")  // convert maf convention deletion to vcf convention         G - => CG C
            {
                new_variant.pos --;
                new_variant.end_pos --;
                char prev_ref = reference_sequence[new_variant.chrom][new_variant.pos];  // get allele before the current position
                new_variant.ref = prev_ref + new_variant.ref;
                new_variant.alt = prev_ref;
            }
            new_variant.insertion = (new_variant.alt.length() > new_variant.ref.length());
            new_variant.deletion = (new_variant.alt.length() < new_variant.ref.length());
            variant_vec.push_back(new_variant);
            variant_count ++;
        }
        my_vf.close();
        cout << "[INFO] " << variant_count << " variants has been loaded from file: " << input_file_names[file_index] << endl;
    }
}


void sortAndIndexVariant(vector<VariantEntry>& variant_vec, vector<pair<size_t, size_t> >& variant_block_vec) // sort and bin the variants, each thread takes a bin to process each time, bin determined by maximum_variant_block_size and maximum_variant_block_distance,
{
    if(variant_vec.size() == 0)
    {
        cerr << "[ERROR] No variant need to be processed " << endl;
        exit(1);
    }
    sort(variant_vec.begin(), variant_vec.end(), sortVariantByPos);
    size_t start_index = 0;
    size_t end_index = 0;
    size_t cur_num_variant = 0;
    map<pair< pair<string, int>, pair<string, string> >, VariantEntry*> duplicate_variant_ptr_map;
    for(size_t i = 0; i < variant_vec.size(); i++)
    {
        pair< pair<string, int>, pair<string, string> > variant_key = make_pair(make_pair(variant_vec[i].chrom, variant_vec[i].pos), make_pair(variant_vec[i].ref, variant_vec[i].alt));
        if(duplicate_variant_ptr_map.find(variant_key) == duplicate_variant_ptr_map.end()) // there is no duplicate variant yet
        {
            duplicate_variant_ptr_map[variant_key] = &variant_vec[i];
        }
        else // there is already a duplicate variant
        {
            variant_vec[i].duplicate_variant_ptr = duplicate_variant_ptr_map[variant_key];
        }
        if((cur_num_variant >= maximum_variant_block_size) || (i != start_index && (variant_vec[i].chrom != variant_vec[start_index].chrom || variant_vec[i].pos - variant_vec[start_index].pos > maximum_variant_block_distance)))
        {
            end_index = i - 1;
            variant_block_vec.push_back(make_pair(start_index, end_index));
            start_index = i;
            cur_num_variant = 0;
        }
        cur_num_variant ++;
    }
    variant_block_vec.push_back(make_pair(start_index, variant_vec.size() - 1));
    cout << "[INFO] " <<  duplicate_variant_ptr_map.size() << " out of " << variant_vec.size() << " variants are unique" << endl;
}

void sortAndIndexVariant16K(vector<VariantEntry>& variant_vec, vector<pair<size_t, size_t> >& variant_block_vec) // sort and bin the variants, each thread takes a bin to process each time, bin determined by BIN_SIZE(smallest bam file block)
{
    if(variant_vec.size() == 0)
    {
        cerr << "[ERROR] No variant need to be processed " << endl;
        exit(1);
    }
    sort(variant_vec.begin(), variant_vec.end(), sortVariantByPos);
    size_t start_index = 0;
    size_t end_index = 0;
    size_t cur_num_variant = 0;
    map<pair< pair<string, int>, pair<string, string> >, VariantEntry*> duplicate_variant_ptr_map;
    for(size_t i = 0; i < variant_vec.size(); i++)
    {
        pair< pair<string, int>, pair<string, string> > variant_key = make_pair(make_pair(variant_vec[i].chrom, variant_vec[i].pos), make_pair(variant_vec[i].ref, variant_vec[i].alt));
        if(duplicate_variant_ptr_map.find(variant_key) == duplicate_variant_ptr_map.end()) // there is no duplicate variant yet
        {
            duplicate_variant_ptr_map[variant_key] = &variant_vec[i];
        }
        else // there is already a duplicate variant
        {
            variant_vec[i].duplicate_variant_ptr = duplicate_variant_ptr_map[variant_key];
        }
        if(i != start_index && (variant_vec[i].chrom != variant_vec[start_index].chrom || (variant_vec[i].pos / BIN_SIZE) != (variant_vec[start_index].pos / BIN_SIZE)) )
        {
            end_index = i - 1;
            variant_block_vec.push_back(make_pair(start_index, end_index));
            start_index = i;
            cur_num_variant = 0;
        }
        cur_num_variant ++;
    }
    variant_block_vec.push_back(make_pair(start_index, variant_vec.size() - 1));
    cout << "[INFO] " <<  duplicate_variant_ptr_map.size() << " out of " << variant_vec.size() << " variants are unique" << endl;
}


/*void printCounts(vector<VariantEntry>& variant_vec)
{
    ofstream output_fs(output_file.c_str());
    if(!output_fs)
    {
        cerr << "[ERROR] fail to open output file: " << output_file << endl;
        exit(1);
    }
    output_fs << "chrom" << "\t" << "pos" << "\t" << "ref" << "\t" << "alt";
    for(map<string, map<char, int> >::iterator it = variant_vec[0].base_count.begin(); it != variant_vec[0].base_count.end(); it++)
    {
        output_fs << "\t" << it->first;
    }
    output_fs << endl;
    
    for(size_t i = 0; i < variant_vec.size(); i++)
    {
        output_fs << variant_vec[i].chrom << "\t" << (variant_vec[i].pos + 1) << "\t" << variant_vec[i].ref << "\t" << variant_vec[i].alt;
        for(map<string, map<char, int> >::iterator it = variant_vec[i].base_count.begin(); it != variant_vec[i].base_count.end(); it++)
        {
            int dp_count = 0;
            int rd_count = 0;
            int ad_count = 0;
            for(map<char, int>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++)
            {
                dp_count += it2->second;
                if(variant_vec[i].insertion || variant_vec[i].deletion)
                {
                    if(toupper(it2->first) == 'M')
                        ad_count += it2->second;
                    else if(toupper(it2->first) == 'R')
                        rd_count += it2->second;
                }
                else if(toupper(it2->first) == variant_vec[i].ref[0]) // ref base
                {
                    rd_count += it2->second;
                }
                else if(toupper(it2->first) == variant_vec[i].alt[0])
                {
                    ad_count += it2->second;
                }
            }
            double vf_count = 0;
            if(dp_count > 0)
                vf_count = (double)ad_count / (double)dp_count;
            output_fs << "\t" << "DP=" << dp_count << ";RD=" << rd_count << ";AD=" << ad_count << ";VF=" << vf_count;
        }
        output_fs << endl;
    }
    output_fs.close();
}*/

void printCountsSIMPLE(vector<VariantEntry>& variant_vec) // print counts for vcf-like variants
{
    ofstream output_fs(output_file.c_str());
    if(!output_fs)
    {
        cerr << "[ERROR] fail to open output file: " << output_file << endl;
        exit(1);
    }
    output_fs << "Chrom" << "\t" << "Pos" << "\t" << "Ref" << "\t" << "Alt";
    for(map<string, map<char, int> >::iterator it = variant_vec[0].base_count.begin(); it != variant_vec[0].base_count.end(); it++)
    {
        output_fs << "\t" << it->first;
    }
    output_fs << endl;
    
    for(size_t i = 0; i < variant_vec.size(); i++)
    {
        output_fs << variant_vec[i].chrom << "\t" << (variant_vec[i].pos + 1 )<< "\t" << variant_vec[i].ref << "\t" << variant_vec[i].alt;
        map<string, map<char, int> >* base_count_ptr = NULL;
        if(variant_vec[i].duplicate_variant_ptr != NULL)
            base_count_ptr = &(variant_vec[i].duplicate_variant_ptr->base_count); // there is another duplicate variant that has the count
        else
            base_count_ptr = &(variant_vec[i].base_count);                        // current variant has the count
        for(map<string, map<char, int> >::iterator it = base_count_ptr->begin(); it != base_count_ptr->end(); it++)
        {
            int dp_count = 0;
            int rd_count = 0;
            int ad_count = 0;
            for(map<char, int>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++)
            {
                dp_count += it2->second;
                if(variant_vec[i].insertion || variant_vec[i].deletion) // insertion or deletion
                {
                    if(toupper(it2->first) == 'M')    // M for exact matched indel
                        ad_count += it2->second;
                    else if(toupper(it2->first) == 'R')  // R for matched ref
                        rd_count += it2->second;
                    //else if(toupper(it2->first) == 'U') // U for unmatched indel, snp
                }
                else if(toupper(it2->first) == variant_vec[i].ref[0]) // ref allele for snp
                {
                    rd_count += it2->second;
                }
                else if(toupper(it2->first) == variant_vec[i].alt[0]) // alt allele for snp
                {
                    ad_count += it2->second;
                }
                //else                                                // other unmatched snp, indel
            }
            double vf_count = 0;
            if(dp_count > 0)
                vf_count = (double)ad_count / (double)dp_count;
            output_fs << "\t" << "DP=" << dp_count << ";RD=" << rd_count << ";AD=" << ad_count << ";VF=" << vf_count;
        }
        output_fs << endl;
    }
    output_fs.close();
}

void printCountsFILLOUT(vector<VariantEntry>& variant_vec) // print counts for tcga maf variants
{
    ofstream output_fs(output_file.c_str());
    if(!output_fs)
    {
        cerr << "[ERROR] fail to open output file: " << output_file << endl;
        exit(1);
    }
    output_fs << "Sample\tNormalUsed\tChrom\tStart\tRef\tAlt\tVariantClass\tGene\tExon\tCall_Confidence\tComments\tTranscriptID\tcDNAchange\tAAchange\tdbSNP_ID\tCosmic_ID\t1000G_MAF\tFailureReason\tN_TotalDepth\tN_RefCount\tN_AltCount\tN_AltFreq\tT_TotalDepth\tT_RefCount\tT_AltCount\tT_AltFreq\tT_Ref+\tT_Ref-\tT_Alt+\tT_Alt-\tAll_N_Aggregate_AlleleDepth\tAll_N_Median_AlleleFreq\tT_freq/All_N_Freq\tOccurence_in_Normals";
    for(map<string, map<char, int> >::iterator it = variant_vec[0].base_count.begin(); it != variant_vec[0].base_count.end(); it++)
    {
        output_fs << "\t" << it->first;
    }
    output_fs << endl;
    
    for(size_t i = 0; i < variant_vec.size(); i++)
    {
        output_fs << variant_vec[i].tumor_sample << "\t" << variant_vec[i].normal_sample << "\t" << variant_vec[i].chrom << "\t" << (variant_vec[i].pos + 1)
        << "\t" << variant_vec[i].ref << "\t" << variant_vec[i].alt << "\t" << variant_vec[i].effect << "\t" << variant_vec[i].gene << "\t" << "\t" << "\t"
        << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << variant_vec[i].n_ref_count << "\t" << variant_vec[i].n_alt_count << "\t"
        << "\t" << "\t" << variant_vec[i].t_ref_count << "\t" << variant_vec[i].t_alt_count << "\t"  << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t"
        << "\t";
        map<string, map<char, int> >* base_count_ptr = NULL;
        if(variant_vec[i].duplicate_variant_ptr != NULL)
            base_count_ptr = &(variant_vec[i].duplicate_variant_ptr->base_count);
        else
            base_count_ptr = &(variant_vec[i].base_count);
        for(map<string, map<char, int> >::iterator it = base_count_ptr->begin(); it != base_count_ptr->end(); it++)
        {
            int dp_count = 0;
            int rd_count = 0;
            int ad_count = 0;
            for(map<char, int>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++)
            {
                dp_count += it2->second;
                if(variant_vec[i].insertion || variant_vec[i].deletion)
                {
                    if(toupper(it2->first) == 'M')                     // M for exact matched indel
                        ad_count += it2->second;
                    else if(toupper(it2->first) == 'R')                // R for matched ref
                        rd_count += it2->second;
                    //else if(toupper(it2->first) == 'U')              // U for unmatched indel, snp
                }
                else if(toupper(it2->first) == variant_vec[i].ref[0])  // ref allele for snp
                {
                    rd_count += it2->second;
                }
                else if(toupper(it2->first) == variant_vec[i].alt[0])  // alt allele for snp
                {
                    ad_count += it2->second;
                }
                //else                                                // other unmatched snp, indel
            }
            double vf_count = 0;
            if(dp_count > 0)
                vf_count = (double)ad_count / (double)dp_count;
            output_fs << "\t" << "DP=" << dp_count << ";RD=" << rd_count << ";AD=" << ad_count << ";VF=" << vf_count;
        }
        output_fs << endl;
    }
    output_fs.close();
}

void baseCountSNP(VariantEntry& my_variant_entry, vector<BamAlignment>& bam_vec, size_t& start_bam_index, string sample_name)  // get counts for snps
{
    bool start_bam_index_set = false;
    for(size_t bam_index = start_bam_index; bam_index < bam_vec.size(); bam_index++)
    {
        BamAlignment &my_bam_alignment = bam_vec[bam_index];
        if(my_bam_alignment.Position > my_variant_entry.pos)  // my_bam_alignment.Position is 0-indexed
            break;
        if(my_bam_alignment.GetEndPosition(false, true) < my_variant_entry.pos)  // GetEndPosition(bool padding, bool closed_interval)
            continue;
        if(!start_bam_index_set) // save roll back point
        {
            start_bam_index = bam_index;
            start_bam_index_set = true;
        }
        size_t ref_pos = my_bam_alignment.Position; // the position on the reference sequence for the next base to be processed
        size_t read_pos = 0; // the position on the read for the next base to be processed
        size_t variant_pos = string::npos; // target vcf position in the read
        vector<CigarOp>::const_iterator cigarIter, cigarIterNext;
        for(cigarIter = my_bam_alignment.CigarData.begin() ; cigarIter != my_bam_alignment.CigarData.end(); cigarIter++)         // parse cigar string
        {
            // check if the next postition to be processed is already > the target position of vcf
            if(ref_pos > my_variant_entry.pos)
                break;
            const CigarOp& op = (*cigarIter);
            switch(op.Type)
            {
                case 'M':
                {
                    if(ref_pos + op.Length > my_variant_entry.pos) // target vcf postion fall in this cigar region
                    {
                        variant_pos = my_variant_entry.pos - ref_pos + read_pos; //calculate the index of target vcf position in the read
                    }
                    ref_pos += op.Length;
                    read_pos += op.Length;
                    break;
                }
                case 'I':
                {
                    read_pos += op.Length;
                    break;
                }
                case 'S' :
                {
                    read_pos += op.Length;
                    break;
                }
                case 'D': case 'N':
                {
                    ref_pos += op.Length;
                    break;
                }
                case 'H':
                    break;
                default:
                    break;
            }
        }
        char cur_alt;
        int cur_bq;
        if(variant_pos == string::npos) // variant pos fall into the deletion region of the alignment
        {
            cur_alt = 'U';
            cur_bq = base_quality_threshold; // ignore base quality threshold
        }
        else // alt or ref
        {
            cur_alt = my_bam_alignment.QueryBases[variant_pos];
            cur_bq = my_bam_alignment.Qualities[variant_pos];
        }
        
        if(cur_bq >= base_quality_threshold) // possible base values for snp: A T G C a t g c U u
        {
            if(my_bam_alignment.IsReverseStrand())
                my_variant_entry.base_count[sample_name][tolower(cur_alt)] ++;  // lower case for - strand
            else
                my_variant_entry.base_count[sample_name][toupper(cur_alt)] ++;  // upper case for + strand
        }
    }

}

// DMP count method, for indels, count the ref depth and total read depth at 1 base ahead (vcf_pos + 1)
void baseCountIndelDMP(VariantEntry& my_variant_entry, vector<BamAlignment>& bam_vec, size_t& start_bam_index, string sample_name, map<string, string>& reference_sequence)
{
    bool start_bam_index_set = false;
    for(size_t bam_index = start_bam_index; bam_index < bam_vec.size(); bam_index++)
    {
        BamAlignment &my_bam_alignment = bam_vec[bam_index];
        if(my_bam_alignment.Position > my_variant_entry.pos + 1)  // my_bam_alignment.Position is 0-indexed
            break;
        if(my_bam_alignment.GetEndPosition(false, true) < my_variant_entry.pos)  // GetEndPosition(bool padding, bool closed_interval)
            continue;
        // any alignment that overlap with my_variant_entry.pos or (my_variant_entry.pos + 1) will pass through here
        if(!start_bam_index_set) // save roll back point
        {
            start_bam_index = bam_index;
            start_bam_index_set = true;
        }
        size_t ref_pos = my_bam_alignment.Position; // the position on the reference sequence for the next base to be processed
        size_t read_pos = 0; // the position on the read for the next base to be processed
        size_t variant_pos = string::npos; // target vcf position in the read
        bool matched_indel = false;
        bool unmatched_indel = false;
        bool do_not_count = false;
        bool break_loop = false;
        vector<CigarOp>::const_iterator cigarIter, cigarIterNext;
        for(cigarIter = my_bam_alignment.CigarData.begin() ; cigarIter != my_bam_alignment.CigarData.end(); cigarIter++)         // parse cigar string
        {
            // check if the next postition to be processed is already > the target position of vcf + 1
            if(break_loop || ref_pos > my_variant_entry.pos + 1)
                break;
            const CigarOp& op = (*cigarIter);
            switch(op.Type)
            {
                case 'M':
                {
                    if(ref_pos + op.Length - 1 == my_variant_entry.pos) // vcf postion is the last base of a M cigar string
                    {
                        // check if there is an indel right next to current position
                        cigarIterNext = cigarIter + 1;
                        if(cigarIterNext != my_bam_alignment.CigarData.end())
                        {
                            const CigarOp& opNext = (*cigarIterNext);
                            if(opNext.Type == 'I')
                            {
                                if(my_variant_entry.insertion && opNext.Length == (my_variant_entry.alt.length() - my_variant_entry.ref.length()) && my_bam_alignment.QueryBases.substr(read_pos + op.Length, opNext.Length) == my_variant_entry.alt.substr(my_variant_entry.ref.length()) )
                                {
                                    matched_indel = true;
                                    break_loop = true;
                                }
                                else //if(my_variant_entry.insertion)
                                {
                                    unmatched_indel = true;
                                    break_loop = true;
                                }
                            }
                            else if(opNext.Type == 'D')
                            {
                                if(my_variant_entry.deletion && opNext.Length == (my_variant_entry.ref.length() - my_variant_entry.alt.length()) )
                                {
                                    matched_indel = true;
                                    break_loop = true;
                                }
                                else
                                {
                                    unmatched_indel = true;
                                    break_loop = true;
                                }
                            }
                            else //if(opNext.Type == 'S' || opNext.Type == 'H') // alignment ends at vcf position
                            {
                                do_not_count = true;
                                break_loop = true;
                            }
                        }
                        else  // alignment ends at vcf position
                        {
                            do_not_count = true;
                            break_loop = true;
                        }
                        if(break_loop)
                            variant_pos = my_variant_entry.pos - ref_pos + read_pos; //calculate the index of target vcf position in the read
                    }
                    else if(ref_pos + op.Length - 1 > my_variant_entry.pos) //ref
                    {
                        variant_pos = my_variant_entry.pos + 1 - ref_pos + read_pos; //calculate the index of target vcf position + 1 in the read
                        break_loop = true;
                    }
                    ref_pos += op.Length;
                    read_pos += op.Length;
                    break;
                }
                case 'I':
                {
                    read_pos += op.Length;
                    break;
                }
                case 'S' :
                {
                    read_pos += op.Length;
                    break;
                }
                case 'D': case 'N':
                {
                    ref_pos += op.Length;
                    break;
                }
                case 'H':
                    break;
                default:
                    break;
            }
        }
        
        char cur_alt;
        int cur_bq;
        if(do_not_count)
            continue;
        if(variant_pos == string::npos || unmatched_indel) // noise
        {
            cur_alt = 'U';
            cur_bq = base_quality_threshold;   //ignore base quality threshold for this case
        }
        else if(matched_indel) //alt
        {
            cur_alt = 'M';
            cur_bq = my_bam_alignment.Qualities[variant_pos]; // apply base quality threshold to the base right before indel
        }
        else // ref
        {
            cur_alt = my_bam_alignment.QueryBases[variant_pos];
            cur_bq = my_bam_alignment.Qualities[variant_pos];
            if(cur_alt == reference_sequence[my_variant_entry.chrom][my_variant_entry.pos + 1])
                cur_alt = 'R';   // ref
            else
                cur_alt = 'U';   // noise
        }
        if(cur_bq >= base_quality_threshold) // possible base values for indel: U M R u m r
        {
            if(my_bam_alignment.IsReverseStrand())
                my_variant_entry.base_count[sample_name][tolower(cur_alt)] ++;  // lower case for - strand
            else
                my_variant_entry.base_count[sample_name][toupper(cur_alt)] ++;  // upper case for + strand
        }
    }

}


// to be implemented
void baseCountIndelBIC(VariantEntry& my_variant_entry, vector<BamAlignment>& bam_vec, size_t& start_bam_index, string sample_name)
{}




void getBaseCounts()
{
    vector<VariantEntry> variant_vec;
    vector<pair<size_t, size_t> >variant_block_vec;
    
    map<string, string> reference_sequence;
    loadReferenceSequenceSpeedup(input_fasta_file, reference_sequence);
    
    if(input_variant_is_maf)
        loadVariantFileMAF(input_variant_files, variant_vec, reference_sequence);
    else if(input_variant_is_vcf)
        loadVariantFileVCF(input_variant_files, variant_vec);
    sortAndIndexVariant(variant_vec, variant_block_vec);      //sortAndIndexVariant16K(variant_vec, variant_block_vec);

    for(map<string, string>::iterator it_bam = input_bam_files.begin(); it_bam != input_bam_files.end(); it_bam ++)
    {
        size_t variant_block_index = 0;
        cout << "[INFO] Processing bam file: " << it_bam->second << endl;
#pragma omp parallel num_threads(num_thread)
        {
            int thread_num = omp_get_thread_num();
            BamReader my_bam_reader;
            if(!my_bam_reader.Open(it_bam->second))
            {
#pragma omp critical(output_stderr)
                {
                    cerr << "[ERROR] fail to open input bam file: " << it_bam->second << endl;
                }
                exit(1);
            }
            string input_bam_index_file1 = it_bam->second.substr(0, it_bam->second.length() - 3) + "bai";
            string input_bam_index_file2 = it_bam->second + ".bai";
            if(!my_bam_reader.OpenIndex(input_bam_index_file1))
            {
                if(!my_bam_reader.OpenIndex(input_bam_index_file2))
                {
#pragma omp critical(output_stderr)
                    {
                        cerr << "[ERROR] fail to open input bam index file: " << input_bam_index_file1 << ", or " << input_bam_index_file2 << endl;
                    }
                    exit(1);
                }
            }
            while(variant_block_index < variant_block_vec.size())
            {
                size_t variant_start_index = -1;
                size_t variant_end_index = -1;
#pragma omp critical(get_next_variant_block)
                if(variant_block_index < variant_block_vec.size())
                {
                    variant_start_index = variant_block_vec[variant_block_index].first;
                    variant_end_index = variant_block_vec[variant_block_index].second;
                    variant_block_index ++;
                    //cout << "Thread " << thread_num << " getting assigned block: " << (variant_block_index - 1) << " of total " << variant_block_vec.size() << " blocks. Range from " << variant_start_index << " to " << variant_end_index << endl;
                }
                if(variant_start_index != -1 && variant_end_index != -1)
                {
                    int refid1 = my_bam_reader.GetReferenceID(variant_vec[variant_start_index].chrom);
                    int refid2 = my_bam_reader.GetReferenceID(variant_vec[variant_end_index].chrom);
                    if(refid1 == -1)
                    {
                        cerr << "[Error] Could not find variant chrom: " << variant_vec[variant_start_index].chrom << " in the bam file" << endl;
                        exit(1);
                    }
                    if(refid2 == -1)
                    {
                        cerr << "[Error] Could not find variant chrom: " << variant_vec[variant_end_index].chrom << " in the bam file" << endl;
                        exit(1);
                    }
                    my_bam_reader.SetRegion(refid1, variant_vec[variant_start_index].pos, refid2, variant_vec[variant_end_index].pos + 2); // 2bp buffer for indels DMP counting method, may need to extend more for BIC counting method
                    //my_bam_reader.SetRegion(refid1, int(variant_vec[variant_start_index].pos / BIN_SIZE) * BIN_SIZE, refid2, (int(variant_vec[variant_end_index].pos / BIN_SIZE) + 1) * BIN_SIZE);
                    
                    vector<BamAlignment> bam_vec;
                    BamAlignment new_bam_alignment;
                    while(my_bam_reader.GetNextAlignment(new_bam_alignment))
                    {
                        if((filter_duplicate && new_bam_alignment.IsDuplicate()) || (filter_improper_pair && !new_bam_alignment.IsProperPair()) || (filter_indel && alignmentHasIndel(new_bam_alignment)) || (filter_qc_failed && new_bam_alignment.IsFailedQC()) || (filter_non_primary && !(new_bam_alignment.IsPrimaryAlignment()) ) ||  new_bam_alignment.MapQuality < mapping_quality_threshold)
                            continue; //apply bam filtering
                        bam_vec.push_back(new_bam_alignment);
                    }
                    size_t start_bam_index = 0;
                    for(size_t variant_index = variant_start_index; variant_index <= variant_end_index ; variant_index++)
                    {
                        if(variant_vec[variant_index].duplicate_variant_ptr != NULL) // there is a duplicate variant has been counted already
                            continue;
                        if(!variant_vec[variant_index].insertion && !variant_vec[variant_index].deletion)
                            baseCountSNP(variant_vec[variant_index], bam_vec, start_bam_index, it_bam->first);
                        else if(count_method == "DMP")
                            baseCountIndelDMP(variant_vec[variant_index], bam_vec, start_bam_index, it_bam->first, reference_sequence);
                        else if(count_method == "BIC")
                            baseCountIndelBIC(variant_vec[variant_index], bam_vec, start_bam_index, it_bam->first);
                    }
                }
            }
            my_bam_reader.Close();
        }
    }
    if(input_variant_is_maf)
        printCountsFILLOUT(variant_vec);
    else if(input_variant_is_vcf)
        printCountsSIMPLE(variant_vec);
    cout << "[INFO] Finished processing all bam files" << endl;
}

int main(int argc, const char * argv[])
{
    parseOption(argc, argv);
    getBaseCounts();
    return 0;
}



