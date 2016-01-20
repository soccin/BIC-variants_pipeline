//
//  ParseRefseqGene.cpp
//  ParseRefseqGene
//
//  Created by Zeng, Zheng/Sloan-Kettering Institute on 7/9/14.
//  Copyright (c) 2014 Zeng, Zheng/Sloan-Kettering Institute. All rights reserved.
//
//
//  CountErrors.cpp
//  CountErrors
//
//  Created by Zeng, Zheng/Sloan-Kettering Institute on 1/10/14.
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
#include <map>
#include <iomanip>


//#define _DEBUG
//#define _PARSING_DEBUG
//#define _FASTA_DEBUG

using namespace std;

const string VERSION = "ParseRefseqGene 1.3.0";

string input_refseq_file;
string input_refseq_canonical_file;
string input_cytoband_file;
string input_dbsnp_file;
string input_custom_gene_file;
string input_custom_bait_file;
string input_custom_tiling_file;
string input_custom_bed_file;
string input_gene_interval_file;
string input_gene_interval_annotated_file;
string input_tiling_interval_file;
string input_tiling_interval_annotated_file;
string input_fp_interval_file;
string input_fp_genotype_file;
string input_gene_coord_file;
string input_gc_bias_file;
string input_target_ilist_file;
string input_bait_ilist_file;
string input_reference;
string input_aa_file;
string output_dir;
bool custom_only = false;
const double MINIMUM_OVERLAP_PERCENT = 0.5;

void printUsage(string msg = "")
{
    cerr << endl;
    cerr << VERSION << endl;
    cerr << "Usage: " << endl;
    cerr << "[REQUIRED ARGUMENTS]" << endl;
    cerr << "\t--reference                  <string>                       Input reference sequence" << endl;
    cerr << "\t--refseq                     <string>                       Input refseq file." << endl;
    cerr << "\t--refseq_canonical           <string>                       Input refseq canonical file." << endl;
    cerr << "\t--cytoband                   <string>                       Input cytoband txt file." << endl;
    cerr << "\t--ucsc_dnsnp                 <string>                       Input UCSC dbsnp common database txt file." << endl;
    cerr << "\t--custom_gene                <string>                       Input Impact+ custom gene interval file(exon+3bp)" << endl;
    cerr << "\t--custom_bait                <string>                       Input Impact+ custom bait interval file" << endl;
    cerr << "\t--custom_tiling              <string>                       Input Impact+ custom tiling interval file" << endl;
    cerr << "\t--custom_bed                 <string>                       Input merged(gene/bait/tiling) Impact+ custom interval bed file" << endl;
    cerr << "\t--gene_interval              <string>                       Standard DMP-IMPACT gene interval file" << endl;
    cerr << "\t--gene_interval_annotated    <string>                       Standard DMP-IMPACT gene interval file with annotation" << endl;
    cerr << "\t--gene_interval              <string>                       Standard DMP-IMPACT gene interval file" << endl;
    cerr << "\t--gene_interval_annotated    <string>                       Standard DMP-IMPACT gene interval file with annotation" << endl;
    cerr << "\t--tiling_interval            <string>                       Standard DMP-IMPACT tiling interval file" << endl;
    cerr << "\t--tiling_interval_annotated  <string>                       Standard DMP-IMPACT tiling interval file with annotation" << endl;
    cerr << "\t--fp_interval                <string>                       Standard DMP-IMPACT FP tiling interval file" << endl;
    cerr << "\t--fp_genotype                <string>                       Standard DMP-IMPACT FP tiling genotype file" << endl;
    cerr << "\t--gene_coord                 <string>                       Standard DMP-IMPACT gene coord file" << endl;
    cerr << "\t--gc_bias                    <string>                       Standard DMP-IMPACT gc bias file" << endl;
    cerr << "\t--canonical_aa               <string>                       Standard DMP-IMPACT caonical exon with amino acid file" << endl;
    cerr << "\t--target_ilist               <string>                       Standard DMP-IMPACT target ilist file" << endl;
    cerr << "\t--bait_ilist                 <string>                       Standard DMP-IMPACT bait ilist file" << endl;
    cerr << "\t--output                     <string>                       Output directory" << endl;
    cerr << "\t--custom_only                                               Use custom gene only, no standard DMP-IMPACT genes" << endl;
    cerr << "\t--help                                                      Print command line usage" << endl;
    cerr << endl;
    cerr << endl;
    if(!msg.empty())
        cerr << msg << endl;
    exit(1);
}


static struct option long_options[] =
{
    {"reference",                   required_argument,      0,     's'},
    {"refseq",                      required_argument,      0,     'r'},
    {"refseq_canonical",            required_argument,      0,     'a'},
    {"cytoband",                    required_argument,      0,     'c'},
    {"ucsc_dnsnp",                  required_argument,      0,     'd'},
    {"custom_gene",                 required_argument,      0,     'g'},
    {"custom_bait",                 required_argument,      0,     'w'},
    {"custom_tiling",               required_argument,      0,     'x'},
    {"custom_bed",                  required_argument,      0,     'p'},
    {"gene_interval",               required_argument,      0,     'i'},
    {"gene_interval_annotated",     required_argument,      0,     'n'},
    {"tiling_interval",             required_argument,      0,     'j'},
    {"tiling_interval_annotated",   required_argument,      0,     'k'},
    {"fp_interval",                 required_argument,      0,     'f'},
    {"fp_genotype",                 required_argument,      0,     'm'},
    {"gene_coord",                  required_argument,      0,     'e'},
    {"gc_bias",                     required_argument,      0,     'b'},
    {"canonical_aa",                required_argument,      0,     'q'},
    {"target_ilist",                required_argument,      0,     't'},
    {"bait_ilist",                  required_argument,      0,     'l'},
    {"output",                      required_argument,      0,     'o'},
    {"custom_only",                 no_argument,            0,     'y'},
    {"help",                        no_argument,            0,     'h'},
    {0, 0, 0, 0}
};


void parseOption(int argc, const char* argv[])
{
    if(argc == 1)
        printUsage();
    int next_option;
    int option_index = 0;
    do
    {
        next_option = getopt_long(argc, const_cast<char**>(argv), "s:r:a:c:d:g:w:x:p:i:n:j:k:f:m:e:b:q:t:l:o:yh", long_options, &option_index);
        switch(next_option)
        {
            case 's':
                input_reference = optarg;
                break;
            case 'r':
                input_refseq_file = optarg;
                break;
            case 'a':
                input_refseq_canonical_file = optarg;
                break;
            case 'c':
                input_cytoband_file = optarg;
                break;
            case 'd':
                input_dbsnp_file = optarg;
                break;
            case 'i':
                input_gene_interval_file = optarg;
                break;
            case 'n':
                input_gene_interval_annotated_file = optarg;
                break;
            case 'j':
                input_tiling_interval_file = optarg;
                break;
            case 'k':
                input_tiling_interval_annotated_file = optarg;
                break;
            case 'f':
                input_fp_interval_file = optarg;
                break;
            case 'm':
                input_fp_genotype_file = optarg;
                break;
            case 'e':
                input_gene_coord_file = optarg;
                break;
            case 'b':
                input_gc_bias_file = optarg;
                break;
            case 'q':
                input_aa_file = optarg;
                break;
            case 't':
                input_target_ilist_file = optarg;
                break;
            case 'l':
                input_bait_ilist_file = optarg;
                break;
            case 'g':
                input_custom_gene_file = optarg;
                break;
            case 'w':
                input_custom_bait_file = optarg;
                break;
            case 'x':
                input_custom_tiling_file = optarg;
                break;
            case 'p':
                input_custom_bed_file = optarg;
                break;
            case 'o':
                output_dir = optarg;
                break;
            case 'y':
                custom_only = true;
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
    if(input_refseq_file.empty())
        printUsage("[ERROR]: Please specify input refseq file");
    if(input_refseq_canonical_file.empty())
        printUsage("[ERROR]: Please specify input refseq canonical file");
    if(input_cytoband_file.empty())
        printUsage("[ERROR]: Please specify input cytoband file");
    if(input_gene_interval_file.empty())
        printUsage("[ERROR]: Please specify input gene interval file");
    if(input_gene_interval_annotated_file.empty())
        printUsage("[ERROR]: Please specify input gene interval annotated file");
    if(input_gc_bias_file.empty())
        printUsage("[ERROR]: Please specify input gc bias file");
    if(input_target_ilist_file.empty())
        printUsage("[ERROR]: Please specify input target ilist file");
    if(input_bait_ilist_file.empty())
        printUsage("[ERROR]: Please specify input bait ilist file");
    if(output_dir.empty())
        printUsage("[ERROR]: Please specify output file");
    if(input_custom_bed_file.empty() && input_custom_gene_file.empty())
        printUsage("[ERROR]: Please specify input custom gene file or custom bed file");
    if(input_custom_bed_file.empty() && input_custom_bait_file.empty())
        printUsage("[ERROR]: Please specify input custom bait file or custom bed file");
    if(!input_custom_bed_file.empty() && (!input_custom_gene_file.empty() || !input_custom_bait_file.empty() || !input_custom_tiling_file.empty()))
    {
        printUsage("[ERROR]: --custom_bed is mutually exclusive with --custom_gene/--custom_bait/--custom_tiling");
    }
#ifdef _DEBUG
    cerr << "[DEBUG]: Parsing options complete." << endl;
#endif
}


void split(const string& line, char delim, vector<std::string>& parsed_item, bool ignore_empty_item = false)
{
    stringstream my_ss(line);
    string item;
    while (getline(my_ss, item, delim))
    {
#ifdef _PARSING_DEBUG
        cerr << "[DEBUG]: Parsed item: " << item << endl;
#endif
        if (ignore_empty_item && item.empty())
            continue;    // use to skip empty item
        parsed_item.push_back(item);
    }
}


string IntToString(int numerical)
{
    stringstream ss;
    ss << setfill('0') << setw(2) << numerical;
    return ss.str();
}

string complement(string sequence)           //get the complement of DNA sequence
{
	for(size_t i = 0; i < sequence.length(); i++)
    {
        switch(sequence[i])
        {
            case 'G':
                sequence[i] = 'C';
                break;
            case 'C':
                sequence[i] = 'G';
                break;
            case 'T':
                sequence[i] = 'A';
                break;
            case 'A':
                sequence[i] = 'T';
                break;
        }
    }
    return sequence;
}

class Chromosome
{
public:
    
    Chromosome(): chrom_len(0), chrom_offset(0), base_per_line(0), byte_per_line(0) {}
    
    Chromosome(int _chrom_len, long long int _chrom_offset, int _base_per_line, int _byte_per_line): chrom_len(_chrom_len), chrom_offset(_chrom_offset), base_per_line(_base_per_line), byte_per_line(_byte_per_line) {}
    
    Chromosome(const Chromosome &other): chrom_len(other.chrom_len), chrom_offset(other.chrom_offset), base_per_line(other.base_per_line), byte_per_line(other.byte_per_line)  {}
    
    ~Chromosome() {}
    
    int chrom_len;
    long long int chrom_offset;
    int base_per_line;
    int byte_per_line;
};

class Interval
{
public:
    
    Interval(): start(0), end(0), custom(false), deleted(false) {}
    
    Interval(const Interval& other): chrom(other.chrom), start(other.start), end(other.end), cytoband(other.cytoband), gene_name(other.gene_name), gene_name_backup(other.gene_name_backup), transcript_name(other.transcript_name), target_ilist(other.target_ilist), strand(other.strand), exon_no(other.exon_no), aa_range(other.aa_range), annotation(other.annotation), dbsnp(other.dbsnp), genotype(other.genotype), custom(other.custom), deleted(other.deleted){}
    
    ~Interval() {}
    
    string chrom;
    int start;
    int end;
    string cytoband; // cytoband information
    string gene_name; // gene name
    string gene_name_backup; //gene name from original custom design file
    string transcript_name; // transcript name
    string target_ilist;  //target_ilist_name
    string strand;     // strand
    string exon_no; // exon number
    string aa_range;// amino acid range
    string annotation; // annotation
    string dbsnp;     // dbsnp id
    string genotype;  //observed genotype for dbsnp
    bool custom; // whether it's custom intervals
    bool deleted; //whether it's marked as deleted
};

bool operator<(const Interval& lhs, const Interval& rhs)
{
    if(lhs.chrom != rhs.chrom)
        return lhs.chrom < rhs.chrom;
    if(lhs.start != rhs.start)
        return lhs.start < rhs.start;
    return lhs.end < rhs.end;
}

bool compare_interval(const Interval& lhs, const Interval& rhs)
{
    if(lhs.chrom != rhs.chrom)
        return lhs.chrom < rhs.chrom;
    if(lhs.start != rhs.start)
        return lhs.start < rhs.start;
    return lhs.end < rhs.end;
}

int compare_overlap(const Interval& lhs, const Interval& rhs)
{
    if(lhs.chrom < rhs.chrom)
        return -1;
    if(lhs.chrom > rhs.chrom)
        return 1;
    if(lhs.end < rhs.start)
        return -1;
    if(lhs.start > rhs.end)
        return 1;
    return 0;
}

int compare_overlap(const int& start1, const int& end1, const int& start2, const int& end2)
{
    if(end1 < start2)
        return -1;
    if(start1 > end2)
        return 1;
    return 0;
}

bool is_number(const string& s)
{
    string::const_iterator it = s.begin();
    while (it != s.end() && isdigit(*it))
    {
        it++;
    }
    return !s.empty() && it == s.end();
}

bool sort_chrom_numerically(const Interval& lhs, const Interval& rhs)
{
    if(lhs.chrom != rhs.chrom)
    {
        bool chrom1_num = is_number(lhs.chrom);
        bool chrom2_num = is_number(rhs.chrom);
        if(chrom1_num && !chrom2_num)
            return true;
        if(!chrom1_num && chrom2_num)
            return false;
        if(chrom1_num && chrom2_num)
            return atoi(lhs.chrom.c_str()) < atoi(rhs.chrom.c_str());
        else //both string
            return lhs.chrom < rhs.chrom;
    }
    if(lhs.start != rhs.start)
        return lhs.start < rhs.start;
    return lhs.end < rhs.end;
}

bool sort_chrom_numerically_tiling(const Interval& lhs, const Interval& rhs)
{
    bool lhs_tiling_probe = (lhs.gene_name.length() >= 3 && lhs.gene_name.substr(0, 3) == "FP_") || (lhs.gene_name.length() >= 7 && lhs.gene_name.substr(0, 7) == "Tiling_");
    bool rhs_tiling_probe = (rhs.gene_name.length() >= 3 && rhs.gene_name.substr(0, 3) == "FP_") || (rhs.gene_name.length() >= 7 && rhs.gene_name.substr(0, 7) == "Tiling_");
    if(lhs_tiling_probe != rhs_tiling_probe)
        return !lhs_tiling_probe;
    if(lhs.chrom != rhs.chrom)
    {
        bool chrom1_num = is_number(lhs.chrom);
        bool chrom2_num = is_number(rhs.chrom);
        
        if(chrom1_num && !chrom2_num)
            return true;
        if(!chrom1_num && chrom2_num)
            return false;
        if(chrom1_num && chrom2_num)
            return atoi(lhs.chrom.c_str()) < atoi(rhs.chrom.c_str());
        else //both string
            return lhs.chrom < rhs.chrom;
    }
    if(lhs.start != rhs.start)
        return lhs.start < rhs.start;
    return lhs.end < rhs.end;
}

class Custom_Design
{
public:
    Custom_Design() {};
    
    ~Custom_Design()
    {
        if(ref_fs)
            ref_fs.close();
    };
    
    void LoadCytoband(string input_file) // Load ucsc cytoband informaiton
    {
        cerr << "[INFO]: Loading cytoband file: " << input_file << endl;
        ifstream cytoband_fs(input_file.c_str());
        if(!cytoband_fs)
        {
            cerr << "[ERROR]: Failed to open input cytoband file: " << input_file << endl;
            exit(1);
        }
        string line;
        while(getline(cytoband_fs, line))
        {
            vector<string> cytoband_item;
            split(line, '\t', cytoband_item);
            Interval new_cytoband;
            new_cytoband.chrom = cytoband_item[0];
            new_cytoband.start = atoi(cytoband_item[1].c_str());
            new_cytoband.end = atoi(cytoband_item[2].c_str());
            new_cytoband.cytoband = cytoband_item[0] + cytoband_item[3];
            cytoband_vec.push_back(new_cytoband);
        }
        sort(cytoband_vec.begin(), cytoband_vec.end(), compare_interval);
        cytoband_fs.close();
    }

    void LoadAllReferenceSequence(string fasta_filename, map<string, string>& reference_sequence)
    {
        cerr << "[INFO]: Loading reference sequence: " << fasta_filename << endl;
        ifstream ref_fs(fasta_filename.c_str());
        if(!ref_fs)
        {
            cerr << "[ERROR]: fail to open reference fasta file: " << fasta_filename << endl;
            exit(1);
        }
        string fasta_index_filename = fasta_filename + ".fai";
        ifstream index_fs(fasta_index_filename.c_str());
        if(!index_fs)
        {
            cerr << "[ERROR]: fail to open reference fasta index file: " << fasta_index_filename << endl;
            exit(1);
        }
        string line;
        while(getline(index_fs, line))
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
            reference_sequence.insert(make_pair(chrom_name, new_seq));
            reference_sequence[chrom_name].resize(chrom_len);
            ref_fs.seekg(chrom_offset);
            char* seq_buff = new char[byte_len];
            ref_fs.read(seq_buff, byte_len);
            string::iterator it_target = reference_sequence[chrom_name].begin();
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
        ref_fs.close();
        index_fs.close();
    }
    
    void OpenReferenceSequence(string fasta_filename)
    {
        cerr << "[INFO]: Openning reference sequence: " << fasta_filename << endl;
        ref_fs.open(fasta_filename.c_str());
        if(!ref_fs)
        {
            cerr << "[ERROR]: fail to open reference fasta file: " << fasta_filename << endl;
            exit(1);
        }
        string fasta_index_filename = fasta_filename + ".fai";
        ifstream index_fs(fasta_index_filename.c_str());
        if(!index_fs)
        {
            cerr << "[ERROR]: fail to open reference fasta index file: " << fasta_index_filename << endl;
            exit(1);
        }
        string line;
        while(getline(index_fs, line))
        {
            vector<string> index_items;
            split(line, '\t', index_items);
            string chrom_name = index_items[0];
            if(chrom_name.substr(0,3) == "chr")
		chrom_name = chrom_name.substr(3);  // for hg19 genome, remove chr from chromosome name
	    int chrom_len = atoi(index_items[1].c_str());
            long long int chrom_offset = atoll(index_items[2].c_str());
            int base_per_line = atoi(index_items[3].c_str());
            int byte_per_line = atoi(index_items[4].c_str());
            Chromosome new_chrom_info(chrom_len, chrom_offset, base_per_line, byte_per_line);
            chrom_info_vec.insert(make_pair(chrom_name, new_chrom_info));
        }
        index_fs.close();
    }
    
    void GetRefSequence(string chrom, int start, int end, string& fetched_seq)
    {
        int start_offset = start - 1 + ((start - 1)  / chrom_info_vec[chrom].base_per_line) * (chrom_info_vec[chrom].byte_per_line - chrom_info_vec[chrom].base_per_line);
        int end_offset = end - 1 + ((end - 1) / chrom_info_vec[chrom].base_per_line) * (chrom_info_vec[chrom].byte_per_line - chrom_info_vec[chrom].base_per_line);
        int seq_len = end - start + 1;
        int byte_len = end_offset - start_offset + 1 ;
        
        fetched_seq.resize(seq_len);
        char* seq_buff = new char[byte_len];
        
        ref_fs.seekg(chrom_info_vec[chrom].chrom_offset + start_offset);
        ref_fs.read(seq_buff, byte_len);
        
        string::iterator it_target = fetched_seq.begin();
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
        //cerr << "[DEBUG]: " << chrom << ":" << start << "-" << end << endl;
        //cerr << "[DEBUG]: " << fetched_seq << endl;
        //cerr << endl;
    }
    
    void LoadRefseqCanonicalExonAA(string input_file) // load refseq canonical transcript with exon and AA information
    {
        cerr << "[INFO]: Loading refseq canonical file: " << input_file << endl;
        ifstream canonical_fs(input_file.c_str());
        if(!canonical_fs)
        {
            cerr << "[ERROR]: Failed to open input refseq canonical file: " << input_file << endl;
            exit(1);
        }
        string line;
        while(getline(canonical_fs, line))
        {
            vector<string> canonical_item;
            split(line, '\t', canonical_item);
            vector<string> transcript_name_item;
            split(canonical_item[4], ':', transcript_name_item);
            if(transcript_name_item.size() != 5)
            {
                cerr << "[ERROR]: Incorrect format: " << input_file << endl;
                exit(1);
            }
            Interval new_canonical_transcript;
            new_canonical_transcript.chrom = canonical_item[0];
            new_canonical_transcript.start = atoi(canonical_item[1].c_str());
            new_canonical_transcript.end = atoi(canonical_item[2].c_str());
            new_canonical_transcript.strand = canonical_item[3];
            new_canonical_transcript.gene_name = transcript_name_item[0];
            new_canonical_transcript.transcript_name = transcript_name_item[1];
            new_canonical_transcript.exon_no = transcript_name_item[2];
            new_canonical_transcript.aa_range = new_canonical_transcript.gene_name + ":" + new_canonical_transcript.transcript_name + ":" + new_canonical_transcript.exon_no + ":" + transcript_name_item[4];
            refseq_exon_vec.push_back(new_canonical_transcript);  // save exon information
            
            pair<string, string> transcript_key = make_pair(new_canonical_transcript.transcript_name, new_canonical_transcript.chrom);
            pair<string, string> gene_key = make_pair(new_canonical_transcript.gene_name, new_canonical_transcript.chrom);
            refseq_canonical[transcript_key] = new_canonical_transcript.gene_name;  // create hashtable to map refseq transcript name to gene name
            if(gene_longest.find(gene_key) == gene_longest.end())  // keep track of the gene boundary
            {
                gene_longest[gene_key] = make_pair(new_canonical_transcript.start, new_canonical_transcript.end);
            }
            else
            {
                gene_longest[gene_key].first = min(gene_longest[gene_key].first, new_canonical_transcript.start);
                gene_longest[gene_key].second = max(gene_longest[gene_key].second, new_canonical_transcript.end);
            }
        }
        sort(refseq_exon_vec.begin(), refseq_exon_vec.end(), compare_interval);
        canonical_fs.close();
    }

    void LoadRefseqGene(string input_file)
    {
        cerr << "[INFO]: Loading refseq file: " << input_file << endl;
        ifstream refseq_fs(input_file.c_str());
        if(!refseq_fs)
        {
            cerr << "[ERROR]: Failed to open input refseq file: " << input_file << endl;
            exit(1);
        }
        string line;
        
        map<pair<string, string>, string> refseq_gene_to_canonical_gene;
        while(getline(refseq_fs, line))
        {
            vector<string> refseq_item;
            split(line, '\t', refseq_item);
            string gene_name = refseq_item[12];
            string refseq_name = refseq_item[1];
            string refseq_chrom = refseq_item[2];
            pair<string, string> transcript_key = make_pair(refseq_name, refseq_chrom);
            pair<string, string> gene_key = make_pair(gene_name, refseq_chrom);
            if(refseq_canonical.find(transcript_key) != refseq_canonical.end() && gene_name != refseq_canonical[transcript_key])
            {
                refseq_gene_to_canonical_gene[gene_key] = refseq_canonical[transcript_key]; // create hashmap to convert inconsistent gene name in refseq file to the gene name in refseq canonical transcript file
            }
        }
        
        refseq_fs.clear();
        refseq_fs.seekg(0, ios::beg);
        
        while(getline(refseq_fs, line))
        {
            vector<string> refseq_item;
            split(line, '\t', refseq_item);
            string gene_name = refseq_item[12];
            string refseq_name = refseq_item[1];
            string refseq_chrom = refseq_item[2];
            int refseq_txstart = atoi(refseq_item[4].c_str());
            int refseq_txend = atoi(refseq_item[5].c_str());
            pair<string, string> gene_key = make_pair(gene_name, refseq_chrom);
            
            if(refseq_gene_to_canonical_gene.find(gene_key) != refseq_gene_to_canonical_gene.end())
                gene_name = refseq_gene_to_canonical_gene[gene_key]; // convert inconsistent gene name in refseq file to the gene name in refseq canonical transcript file
            if(gene_longest.find(gene_key) != gene_longest.end()) // genes that have canonical transcript
            {
                if(compare_overlap(gene_longest[gene_key].first, gene_longest[gene_key].second, refseq_txstart, refseq_txend) == 0)
                {
                    gene_longest[gene_key].first = min(gene_longest[gene_key].first, refseq_txstart);
                    gene_longest[gene_key].second = max(gene_longest[gene_key].second, refseq_txend);
                }
            }
            else  // genes that don't have canonical transcript
            {
                if(gene_longest_non_canonical.find(gene_key) != gene_longest_non_canonical.end())
                {
                    if(compare_overlap(gene_longest_non_canonical[gene_key].first, gene_longest_non_canonical[gene_key].second, refseq_txstart, refseq_txend) == 0)
                    {
                        gene_longest_non_canonical[gene_key].first = min(gene_longest_non_canonical[gene_key].first, refseq_txstart);
                        gene_longest_non_canonical[gene_key].second = max(gene_longest_non_canonical[gene_key].second, refseq_txend);
                    }
                }
                else
                {
                    gene_longest_non_canonical[gene_key].first = refseq_txstart;
                    gene_longest_non_canonical[gene_key].second = refseq_txend;
                }
            }
        }
        /*map<pair<string, string>, pair<int, int> >::iterator test_it = gene_longest.find(make_pair("KAT7", "17"));
         if(test_it == gene_longest.end())
         cerr << "could not find" << endl;
         else
         cerr << test_it->second.first << "\t" << test_it->second.second << endl;
         */
        map<pair<string, string>, pair<int, int> >::iterator it;
        for(it = gene_longest.begin(); it != gene_longest.end(); it ++)
        {
            Interval new_refseq_gene;
            new_refseq_gene.gene_name = it->first.first;
            new_refseq_gene.chrom = it->first.second;
            new_refseq_gene.start = it->second.first;
            new_refseq_gene.end = it->second.second;
            refseq_gene_vec.push_back(new_refseq_gene);
        }
        sort(refseq_gene_vec.begin(), refseq_gene_vec.end(), compare_interval);
        for(it = gene_longest_non_canonical.begin(); it != gene_longest_non_canonical.end(); it ++)
        {
            Interval new_refseq_gene;
            new_refseq_gene.gene_name = it->first.first;
            new_refseq_gene.chrom = it->first.second;
            new_refseq_gene.start = it->second.first;
            new_refseq_gene.end = it->second.second;
            refseq_gene_non_canonical_vec.push_back(new_refseq_gene);
        }
        sort(refseq_gene_non_canonical_vec.begin(), refseq_gene_non_canonical_vec.end(), compare_interval);
        refseq_fs.close();
    }
    
    void LoadUcscDBsnp(string input_file)
    {
        cerr << "[INFO]: Loading input ucsc dbsnp file (this may take several minutes!): " << input_file << endl;
        ifstream input_fs(input_file.c_str());
        if(!input_fs)
        {
            cerr << "[ERROR]: Failed to open input ucsc dbsnp file: " << input_file << endl;
            exit(1);
        }
        string line;
        while(getline(input_fs, line))
        {
            if(line[0] == '#')
                continue;
            vector<string> dbsnp_items;
            split(line, '\t', dbsnp_items);
            pair<string, string> index_key = make_pair(dbsnp_items[0], dbsnp_items[1]);
            map<pair<string, string>, string >::iterator it = ucsc_dbsnp_map.find(index_key);
            if(it != ucsc_dbsnp_map.end())
            {
                if(dbsnp_items[2] == "+")
                    it->second = dbsnp_items[3];
                else
                    it->second = complement(dbsnp_items[3]);  // get complement of the common genotype
            }
        }
    }
    
    void LoadCustomInterval(string input_file, vector<Interval>& custom_interval)
    {
        cerr << "[INFO]: Loading custom interval file: " << input_file << endl;
        ifstream custom_interval_fs(input_file.c_str());
        if(!custom_interval_fs)
        {
            cerr << "[ERROR]: Failed to open custom interval file: " << input_file << endl;
            exit(1);
        }
        string line;
        while(getline(custom_interval_fs, line))
        {
            vector<string> chrom_item;
            vector<string> coordinate_item;
            split(line, ':', chrom_item);
            split(chrom_item[1], '-', coordinate_item);
            Interval new_custom_interval;
            new_custom_interval.chrom = chrom_item[0];
            new_custom_interval.start = atoi(coordinate_item[0].c_str());
            new_custom_interval.end = atoi(coordinate_item[1].c_str());
            new_custom_interval.custom = true;
            custom_interval.push_back(new_custom_interval);
        }
        custom_interval_fs.close();
        sort(custom_interval.begin(), custom_interval.end(), compare_interval);
    }
    
    void LoadCustomTilingInterval(string input_file, vector<Interval>& custom_interval)
    {
        cerr << "[INFO]: Loading custom interval file: " << input_file << endl;
        ifstream custom_interval_fs(input_file.c_str());
        if(!custom_interval_fs)
        {
            cerr << "[ERROR]: Failed to open custom interval file: " << input_file << endl;
            exit(1);
        }
        string line;
        while(getline(custom_interval_fs, line))
        {
            vector<string> dbsnp_item;
            split(line, '\t', dbsnp_item);
            vector<string> chrom_item;
            split(dbsnp_item[0], ':', chrom_item);
            vector<string> coordinate_item;
            split(chrom_item[1], '-', coordinate_item);
            Interval new_custom_interval;
            new_custom_interval.chrom = chrom_item[0];
            new_custom_interval.start = atoi(coordinate_item[0].c_str());
            new_custom_interval.end = atoi(coordinate_item[1].c_str());
            new_custom_interval.dbsnp = dbsnp_item[1];
            new_custom_interval.custom = true;
            custom_interval.push_back(new_custom_interval);
            ucsc_dbsnp_map.insert(make_pair(make_pair(new_custom_interval.chrom, new_custom_interval.dbsnp), ""));
        }
        custom_interval_fs.close();
        sort(custom_interval.begin(), custom_interval.end(), compare_interval);
    }
    
    void LoadCustomBed(string input_file)
    {
        cerr << "[INFO]: Loading custom bed file: " << input_file << endl;
        ifstream custom_bed_fs(input_file.c_str());
        if(!custom_bed_fs)
        {
            cerr << "[ERROR]: Failed to open custom bed file: " << input_file << endl;
            exit(1);
        }
        string line;
        while(getline(custom_bed_fs, line))
        {
            vector<string> custom_item;
            split(line, '\t', custom_item);
            if(custom_item.size() != 6)
            {
            	cerr << "incorrect number of columns in the custom design bed file" << endl;
            	cerr << line << endl;
            	exit(1);
            }
            vector<string> coordinate_item;
            split(custom_item[2], ':', coordinate_item);
            vector<string> position_item;
            split(coordinate_item[1], '-', position_item);

            Interval new_custom_interval;
            new_custom_interval.chrom = coordinate_item[0].substr(3);
            new_custom_interval.start = atoi(position_item[0].c_str());
            new_custom_interval.end = atoi(position_item[1].c_str());
            new_custom_interval.custom = true;
            new_custom_interval.gene_name_backup = custom_item[0];
            string type = custom_item[4];
            if(type == "Category")
            {
                continue;
            }
            else if(type == "Exon" || type == "Gene")
            {
                new_custom_interval.start = max(1, new_custom_interval.start - 3);
                new_custom_interval.end = min(chrom_info_vec[new_custom_interval.chrom].chrom_len, new_custom_interval.end + 3);
                custom_gene_vec.push_back(new_custom_interval);
            }
            else if(type == "PseudoExon")
            {
                custom_gene_vec.push_back(new_custom_interval);
            }
            else if(type == "Bait")
            {
                custom_bait_vec.push_back(new_custom_interval);
            }
            else if(type == "Tiling" || type == "Fingerprint")
            {
                new_custom_interval.dbsnp = custom_item[0];
                custom_tiling_vec.push_back(new_custom_interval);
                ucsc_dbsnp_map.insert(make_pair(make_pair(new_custom_interval.chrom, new_custom_interval.dbsnp), ""));
            }
            else if(type == "Intron")
            {
                custom_intron_vec.push_back(new_custom_interval);
            }
            else if(type == "Promoter")
            {
                custom_promoter_vec.push_back(new_custom_interval);
            }
            else
            {
                cerr << "[WARNING]: Ignoring unrecognized custom interval type: " << line << endl;
            }
        }
        custom_bed_fs.close();
        sort(custom_gene_vec.begin(), custom_gene_vec.end(), compare_interval);
        sort(custom_bait_vec.begin(), custom_bait_vec.end(), compare_interval);
        sort(custom_tiling_vec.begin(), custom_tiling_vec.end(), compare_interval);
    }
    
    string GetCytoband(string chrom, int start, int end)
    {
        for(size_t i = 0; i < cytoband_vec.size(); i++)
        {
            if(cytoband_vec[i].chrom == chrom && end >= cytoband_vec[i].start && start <= cytoband_vec[i].end)
                return cytoband_vec[i].cytoband;
        }
        return "";
    }
    
    double CalculateGC(string chrom, int start, int end)
    {
        if(chrom_info_vec.find(chrom) == chrom_info_vec.end())
        {
            cerr << "[ERROR]: Could not find chromosome " << chrom << " in reference sequence" << endl;
            exit(1);
        }
        if(end - start < 150)
        {
            int midpoint = (end + start) / 2;
            start = max(1, midpoint - 75);
            end = min((int)(chrom_info_vec[chrom].chrom_len), midpoint + 75);
        }
        string interval_seq;
        GetRefSequence(chrom, start, end, interval_seq);
        int gc_count = 0;
        for(size_t i = 0; i < interval_seq.length(); i++)
        {
            if(interval_seq[i] == 'G' || interval_seq[i] == 'C')
                gc_count ++;
        }
        return (double)(gc_count) / (double)(interval_seq.length());
    }
    
    double CalculateOverlapPercent(Interval &lhs, Interval &rhs)
    {
        int overlap_start = max(lhs.start, rhs.start);
        int overlap_end = min(lhs.end, rhs.end);
        return (double)(overlap_end - overlap_start + 1) / (double)(rhs.end - rhs.start + 1);
    }

    void AnnotateCustomGene(string output_file)
    {
        size_t custom_index = 0;
        size_t refseq_index = 0;
        size_t reserved_index = 0;
        bool index_reserved = false;
        map<string, int> genes_with_one_or_more_annotated_exons;
        while(custom_index < custom_gene_vec.size() && refseq_index < refseq_exon_vec.size())
        {
            int relation = compare_overlap(custom_gene_vec[custom_index], refseq_exon_vec[refseq_index]);
            if(relation == 1)
            {
                refseq_index ++;
            }
            else // relation == -1 or 0
            {
                if(!index_reserved)
                {
                    reserved_index = refseq_index;
                    index_reserved = true;
                }
                if(relation == 0)
                {
                    if(CalculateOverlapPercent(custom_gene_vec[custom_index], refseq_exon_vec[refseq_index]) >= MINIMUM_OVERLAP_PERCENT)
                    {
                        stringstream annotation_ss;
                        annotation_ss << refseq_exon_vec[refseq_index].gene_name << ":" << refseq_exon_vec[refseq_index].transcript_name << ":" << refseq_exon_vec[refseq_index].exon_no << ":" << refseq_exon_vec[refseq_index].chrom << ":" << refseq_exon_vec[refseq_index].start << ":" << refseq_exon_vec[refseq_index].end;
                        custom_gene_vec[custom_index].annotation = annotation_ss.str();
                        custom_gene_vec[custom_index].gene_name = refseq_exon_vec[refseq_index].gene_name;
                        //custom_gene_vec[custom_index].aa_range = refseq_exon_vec[refseq_index].aa_range;
                        Interval *refseq_exon_ptr = &(refseq_exon_vec[refseq_index]);
                        custom_aa_map[refseq_exon_ptr] = 1;
                        genes_with_one_or_more_annotated_exons[custom_gene_vec[custom_index].gene_name] = 1;
                    }
                    else
                    {
                        refseq_index ++;
                        continue;
                    }
                }
                custom_gene_vec[custom_index].cytoband = GetCytoband(custom_gene_vec[custom_index].chrom, custom_gene_vec[custom_index].start, custom_gene_vec[custom_index].end);
                custom_index ++;
                refseq_index = reserved_index;  // roll back to last reserved index
                index_reserved = false;
            }
        }
        
        custom_index = 0;
        refseq_index = 0;
        while(custom_index < custom_gene_vec.size() && refseq_index < refseq_gene_vec.size())
        {
            int relation = compare_overlap(custom_gene_vec[custom_index], refseq_gene_vec[refseq_index]);
            if(relation == -1)
                custom_index ++;
            else if(relation == 1)
                refseq_index ++;
            else
            {
                if(custom_gene_vec[custom_index].gene_name.empty())
                    custom_gene_vec[custom_index].gene_name = refseq_gene_vec[refseq_index].gene_name;
                custom_index ++;
            }
        }
        
        custom_index = 0;
        refseq_index = 0;
        while(custom_index < custom_gene_vec.size() && refseq_index < refseq_gene_non_canonical_vec.size())
        {
            int relation = compare_overlap(custom_gene_vec[custom_index], refseq_gene_non_canonical_vec[refseq_index]);
            if(relation == -1)
                custom_index ++;
            else if(relation == 1)
                refseq_index ++;
            else
            {
                if(custom_gene_vec[custom_index].gene_name.empty())
                    custom_gene_vec[custom_index].gene_name = refseq_gene_non_canonical_vec[refseq_index].gene_name;
                custom_index ++;
            }
        }
        
        int exon_count = 0;
        string previous_gene_name = "";
        //bool found_gene_interval_without_annotated_exon = false;
        for(size_t i = 0; i < custom_gene_vec.size(); i++)  // remove genes that does not have any annotated exons
        {
            if(custom_gene_vec[i].gene_name.empty())
            	custom_gene_vec[i].gene_name = custom_gene_vec[i].gene_name_backup;
						if(custom_gene_vec[i].gene_name != previous_gene_name)
            {
                previous_gene_name = custom_gene_vec[i].gene_name;
                exon_count = 0;
            }
            exon_count ++;
            custom_gene_vec[i].target_ilist = custom_gene_vec[i].gene_name + "_target_" + IntToString(exon_count);

            if(genes_with_one_or_more_annotated_exons.find(custom_gene_vec[i].gene_name) == genes_with_one_or_more_annotated_exons.end())
            {
                //custom_gene_vec[i].deleted = true;
                //found_gene_interval_without_annotated_exon = true;
                //cerr << "[WARNING]: Removing gene intervals that does not have any annotated exon: " << custom_gene_vec[i].chrom << ":" << custom_gene_vec[i].start << "-" << custom_gene_vec[i].end << " " << custom_gene_vec[i].gene_name << endl;
            }
        }
        /*if(found_gene_interval_without_annotated_exon)
        {
            cerr << "[ERROR]: Please contact Zheng or Nick for the above error" << endl;
            exit(1);
        }*/
        if(!output_file.empty())
        {
            ofstream output_fs(output_file.c_str());
            if(!output_fs)
            {
                cerr << "[ERROR]: Failed to open output gene interval file: " << output_file << endl;
                exit(1);
            }
            for(size_t i = 0; i < custom_gene_vec.size(); i++)
            {
                output_fs << custom_gene_vec[i].chrom << ":" << custom_gene_vec[i].start << "-" << custom_gene_vec[i].end << "\t" << custom_gene_vec[i].cytoband << "\t" << custom_gene_vec[i].annotation << "\t" << custom_gene_vec[i].gene_name << "\t" << custom_gene_vec[i].target_ilist << "\t" << custom_gene_vec[i].aa_range << endl;
            }
            output_fs.close();
        }
    }
    
    void AnnotateCustomTiling()
    {
        for(size_t custom_index = 0; custom_index < custom_tiling_vec.size(); custom_index++)
        {
            custom_tiling_vec[custom_index].cytoband = GetCytoband(custom_tiling_vec[custom_index].chrom, custom_tiling_vec[custom_index].start, custom_tiling_vec[custom_index].end);
            pair<string, string> dbsnp_idex = make_pair(custom_tiling_vec[custom_index].chrom, custom_tiling_vec[custom_index].dbsnp);
            if(ucsc_dbsnp_map.find(dbsnp_idex) != ucsc_dbsnp_map.end())
            {
                custom_tiling_vec[custom_index].genotype = ucsc_dbsnp_map[dbsnp_idex];
            }
            else
            {
                cerr << "[WARNING]: Could not find tiling probe in dbsnp database: " <<  dbsnp_idex.first << "\t" << dbsnp_idex.second << ", will be ignored"<< endl;
            }
            custom_tiling_vec[custom_index].target_ilist = "Tiling_" + custom_tiling_vec[custom_index].dbsnp;
            custom_tiling_vec[custom_index].gene_name = custom_tiling_vec[custom_index].target_ilist;
        }
    }
    
    bool CheckOverlap(vector<Interval> &lhs, vector<Interval> rhs, bool verbose)
    {
        size_t lhs_index = 0;
        size_t rhs_index = 0;
        size_t reserved_index = 0;
        bool index_reserved = false;
        bool overlap = false;
        while(lhs_index < lhs.size() && rhs_index < rhs.size())
        {
            int relation = compare_overlap(lhs[lhs_index], rhs[rhs_index]);
            if(relation == 1)
            {
                rhs_index ++;
            }
            else // relation == -1 or 0
            {
                if(!index_reserved)
                {
                    reserved_index = rhs_index;
                    index_reserved = true;
                }
                if(relation == 0)
                {
                    cerr << "[ERROR]: Found custom intervals overlapping with standard intervals: " << lhs[lhs_index].chrom << ":" << lhs[lhs_index].start << "-" << lhs[lhs_index].end << " and " << rhs[rhs_index].chrom << ":" << rhs[rhs_index].start << "-" << rhs[rhs_index].end << endl;
                    overlap = true;
                    if(!verbose)
                        break;;
                    rhs_index ++;
                }
                else
                {
                    lhs_index ++;
                    rhs_index = reserved_index;  // roll back to last reserved index
                    index_reserved = false;
                }
            }
        }
        return overlap;
    }
    
    
    bool CheckSelfOverlap(vector<Interval>& lhs, bool verbose)
    {
        bool overlap = false;
        if(lhs.size() != 0)
        {
            for(size_t lhs_index = 0; lhs_index < lhs.size() - 1; lhs_index++)
            {
                if(compare_overlap(lhs[lhs_index], lhs[lhs_index + 1]) == 0)
                {
                    cerr << "[ERROR]: Found two custom intervals overlapping with each other: " << lhs[lhs_index].chrom << ":" << lhs[lhs_index].start << "-" << lhs[lhs_index].end << " and " << lhs[lhs_index + 1].chrom << ":" << lhs[lhs_index + 1].start << "-" << lhs[lhs_index + 1].end << endl;
                    overlap = true;
                    if(!verbose)
                        break;
                }
            }
        }
        return overlap;
    }
    
    void MergeOverlapInterval(Interval &lhs, Interval &rhs, bool verbose)
    {
        if(verbose)
            cerr << "[WARNING]: Merging overlapping intervals: " << lhs.chrom << ":" << lhs.start << "-" << lhs.end << " and " << rhs.chrom << ":" << rhs.start << "-" << rhs.end << endl;
        lhs.end = rhs.end;
        if(!rhs.annotation.empty())
        {
            lhs.annotation += ",";
            lhs.annotation += rhs.annotation;
        }
        rhs.deleted = true;
    }
    
    void GenerateCombinedGeneInterval(string input_file, string output_file)
    {
        cerr << "[INFO]: Loading gene interval file: " << input_file << endl;
        ifstream gene_fs(input_file.c_str());
        if(!gene_fs)
        {
            cerr << "[ERROR]: Failed to open input gene interval file: " << input_file << endl;
            exit(1);
        }

        ofstream output_fs(output_file.c_str());
        if(!output_fs)
        {
            cerr << "[ERROR]: Failed to open output gene interval file: " << output_file << endl;
            exit(1);
        }

        string output_annotated_file = output_file + ".annotated";
        ofstream output_annotated_fs(output_annotated_file.c_str());
        if(!output_annotated_fs)
        {
            cerr << "[ERROR]: Failed to open output gene interval annotated file: " << output_annotated_file << endl;
            exit(1);
        }

        string output_file_chr = output_file + ".chr.list";
        ofstream output_fs_chr(output_file_chr.c_str());
        if(!output_fs_chr)
        {
            cerr << "[ERROR]: Failed to open output gene interval file: " << output_file_chr << endl;
            exit(1);
        }

        string output_annotated_file_chr = output_annotated_file + ".chr.list";
        ofstream output_annotated_fs_chr(output_annotated_file_chr.c_str());
        if(!output_annotated_fs_chr)
        {
            cerr << "[ERROR]: Failed to open output gene interval annotated file: " << output_annotated_file_chr << endl;
            exit(1);
        }

        vector<Interval> gene_vec; //add custom gene interval
        string line;
        if(!custom_only)
        {
            while(getline(gene_fs, line)) // add standard gene interval
            {
                vector<string> gene_item;
                vector<string> chrom_item;
                vector<string> coordinate_item;
                split(line, '\t', gene_item);
                split(gene_item[0], ':', chrom_item);
                split(chrom_item[1], '-', coordinate_item);
                string cytoband = (gene_item.size()) > 1 ? gene_item[1] : "";
                string annotation = (gene_item.size()) > 2 ? gene_item[2] : "";
                Interval new_custom_exon;
                new_custom_exon.chrom = chrom_item[0];
                new_custom_exon.start = atoi(coordinate_item[0].c_str());
                new_custom_exon.end = atoi(coordinate_item[1].c_str());
                new_custom_exon.annotation = annotation;
                new_custom_exon.cytoband = cytoband;
                gene_vec.push_back(new_custom_exon);
            }
        }
        gene_vec.insert(gene_vec.end(), custom_gene_vec.begin(), custom_gene_vec.end());
        sort(gene_vec.begin(), gene_vec.end(), compare_interval);
        
        for(size_t dest_index = 0; dest_index < gene_vec.size(); dest_index ++)
        {
            if(gene_vec[dest_index].deleted)
                continue;
            for(size_t src_index = dest_index + 1; src_index < gene_vec.size(); src_index ++)
            {
                if(gene_vec[src_index].deleted)
                    continue;
                if(compare_overlap(gene_vec[dest_index], gene_vec[src_index]) == 0) //overlap
                {
                    MergeOverlapInterval(gene_vec[dest_index], gene_vec[src_index], true);
                }
                else
                {
                    break;
                }
            }
        }
        sort(gene_vec.begin(), gene_vec.end(), sort_chrom_numerically);
        for(size_t i = 0; i < gene_vec.size(); i++)
        {
            if(!gene_vec[i].deleted)
            {
                output_fs << gene_vec[i].chrom << ":" << gene_vec[i].start << "-" << gene_vec[i].end << endl;
                output_annotated_fs << gene_vec[i].chrom << ":" << gene_vec[i].start << "-" << gene_vec[i].end << "\t" << gene_vec[i].cytoband << "\t" << gene_vec[i].annotation << endl;
                output_fs_chr << "chr" << gene_vec[i].chrom << ":" << gene_vec[i].start << "-" << gene_vec[i].end << endl;
                output_annotated_fs_chr << "chr" << gene_vec[i].chrom << ":" << gene_vec[i].start << "-" << gene_vec[i].end << "\t" << gene_vec[i].cytoband << "\t" << gene_vec[i].annotation << endl;
            }
        }
        gene_fs.close();
        output_fs.close();
        output_annotated_fs.close();
        output_fs_chr.close();
        output_annotated_fs_chr.close();
    }
    
    void GenerateCombinedGeneCoord(string input_file, string output_file)
    {
        cerr << "[INFO]: Loading gene coord file: " << input_file << endl;
        ifstream gene_fs(input_file.c_str());
        if(!gene_fs)
        {
            cerr << "Error: Failed to open input gene coord file: " << input_file << endl;
            exit(1);
        }
        ofstream output_fs(output_file.c_str());
        if(!output_fs)
        {
            cerr << "Error: Failed to open output gene coord file: " << input_file << endl;
            exit(1);
        }
        
        vector<Interval> gene_vec; //add custom gene interval
        string line;
        if(!custom_only)
        {
            while(getline(gene_fs, line)) // add standard gene interval
            {
                vector<string> gene_item;
                split(line, '\t', gene_item);
                Interval new_custom_exon;
                new_custom_exon.chrom = gene_item[2];
                new_custom_exon.start = atoi(gene_item[3].c_str());
                new_custom_exon.end = atoi(gene_item[4].c_str());
                new_custom_exon.gene_name = gene_item[1];
                gene_vec.push_back(new_custom_exon);
            }
        }
        bool gene_info_not_found = false;
        gene_vec.insert(gene_vec.end(), custom_gene_vec.begin(), custom_gene_vec.end());
        sort(gene_vec.begin(), gene_vec.end(), sort_chrom_numerically);
        for(size_t i = 0; i < gene_vec.size(); i++)
        {
            if(gene_vec[i].gene_name.empty())
            {
                cerr << "[ERROR]: Could not find gene information for: " << gene_vec[i].chrom << "\t" << gene_vec[i].start << "\t" << gene_vec[i].end << endl;
                gene_info_not_found = true;
            }
            if(!gene_vec[i].deleted)
            {
                output_fs << i << "\t" << gene_vec[i].gene_name << "\t" << gene_vec[i].chrom << "\t" << gene_vec[i].start << "\t" << gene_vec[i].end << endl;
            }
        }
        gene_fs.close();
        output_fs.close();
        if(gene_info_not_found)
            exit(1);
    }
    
    void GenerateCombinedGCBias(string input_file, string output_file, int tiling_left_buffer, int tiling_right_buffer)
    {
        cerr << "[INFO]: Loading gc bias file: " << input_file << endl;
        ifstream input_fs(input_file.c_str());
        if(!input_fs)
        {
            cerr << "[ERROR]: Failed to open input gc bias file: " << input_file << endl;
            exit(1);
        }
        ofstream output_fs(output_file.c_str());
        if(!output_fs)
        {
            cerr << "[ERROR]: Failed to open output gc bias file: " << input_file << endl;
            exit(1);
        }
        vector<Interval> gene_vec; //add custom gene interval
        string header;
        getline(input_fs, header);
        string line;
        if(!custom_only)
        {
            while(getline(input_fs, line)) // add standard gene interval
            {
                vector<string> gc_item;
                vector<string> chrom_item;
                vector<string> coordinate_item;
                split(line, '\t', gc_item);
                split(gc_item[3], ':', chrom_item);
                split(chrom_item[1], '-', coordinate_item);
                Interval new_custom_exon;
                new_custom_exon.chrom = chrom_item[0];
                new_custom_exon.start = atoi(coordinate_item[0].c_str());
                new_custom_exon.end = atoi(coordinate_item[1].c_str());
                new_custom_exon.gene_name = gc_item[2];
                //new_custom_exon.gc_content = atof(gc_item[6].c_str());
                gene_vec.push_back(new_custom_exon);
            }
        }
        gene_vec.insert(gene_vec.end(), custom_gene_vec.begin(), custom_gene_vec.end());
        sort(gene_vec.begin(), gene_vec.end(), compare_interval);
        for(size_t dest_index = 0; dest_index < gene_vec.size(); dest_index ++)
        {
            if(gene_vec[dest_index].deleted)
                continue;
            for(size_t src_index = dest_index + 1; src_index < gene_vec.size(); src_index ++)
            {
                if(gene_vec[src_index].deleted)
                    continue;
                if(compare_overlap(gene_vec[dest_index], gene_vec[src_index]) == 0) //overlap
                {
                    MergeOverlapInterval(gene_vec[dest_index], gene_vec[src_index], false);
                }
                else
                {
                    break;
                }
            }
        }
        
        for(size_t i = 0; i < custom_tiling_vec.size(); i ++)
        {
            gene_vec.push_back(custom_tiling_vec[i]);
            gene_vec[gene_vec.size() - 1].start = max(1, gene_vec[gene_vec.size() - 1].start - tiling_left_buffer);
            gene_vec[gene_vec.size() - 1].end = min(chrom_info_vec[gene_vec[gene_vec.size() - 1].chrom].chrom_len, gene_vec[gene_vec.size() - 1].end + tiling_right_buffer);
        }
        sort(gene_vec.begin(), gene_vec.end(), sort_chrom_numerically_tiling);
        output_fs << header << endl;
        for(size_t i = 0; i < gene_vec.size(); i++)
        {
            if(!gene_vec[i].deleted)
            {
                double calculated_gc_content = CalculateGC(gene_vec[i].chrom, gene_vec[i].start, gene_vec[i].end);
                output_fs << (i+1) << "\t" << gene_vec[i].chrom << "\t" << gene_vec[i].gene_name << "\t" << gene_vec[i].chrom << ":" << gene_vec[i].start << "-" << gene_vec[i].end << "\t" << "chr" << gene_vec[i].chrom << "\t" << gene_vec[i].start << "\t" << fixed << setprecision(16) << calculated_gc_content << endl;
            }
        }
        input_fs.close();
        output_fs.close();
    }
    
    void GenerateCombinedTargetIlist(string input_file, string output_file)
    {
        
        cerr << "[INFO]: Loading target ilist file: " << input_file << endl;
        ifstream input_fs(input_file.c_str());
        if(!input_fs)
        {
            cerr << "[ERROR]: Failed to open target ilist file: " << input_file << endl;
            exit(1);
        }
        ofstream output_fs(output_file.c_str());
        if(!output_fs)
        {
            cerr << "[ERROR]: Failed to open output target ilist file: " << input_file << endl;
            exit(1);
        }
        
        vector<Interval> gene_vec;
        string line;
        while(getline(input_fs, line)) // add standard gene interval
        {
            if(line[0] == '@')
            {
                output_fs << line << endl;
                continue;
            }
            if(!custom_only)
            {
                vector<string> target_item;
                split(line, '\t', target_item);
                Interval new_custom_exon;
                new_custom_exon.chrom = target_item[0];
                new_custom_exon.start = atoi(target_item[1].c_str());
                new_custom_exon.end = atoi(target_item[2].c_str());
                new_custom_exon.target_ilist = target_item[4];
                gene_vec.push_back(new_custom_exon);
            }
        }
        
        size_t start_index_custom_gene = gene_vec.size();
        gene_vec.insert(gene_vec.end(), custom_gene_vec.begin(), custom_gene_vec.end());   // add custom gene interval
        for(size_t i = start_index_custom_gene; i < gene_vec.size(); i++) // remove 3bp buffer
        {
            gene_vec[i].start += 3;
            gene_vec[i].end -= 3;
        }
        
        size_t start_index_custom_tiling = gene_vec.size();
        gene_vec.insert(gene_vec.end(), custom_tiling_vec.begin(), custom_tiling_vec.end()); // add custom tiling interval
        for(size_t i = start_index_custom_tiling; i < gene_vec.size(); i++) // adjust start position
        {
            gene_vec[i].start -= 1;
        }

        size_t start_index_custom_intron = gene_vec.size();
        gene_vec.insert(gene_vec.end(), custom_intron_vec.begin(), custom_intron_vec.end()); // add custom intron interval
        int intron_count = 0;
        for(size_t i = start_index_custom_intron; i < gene_vec.size(); i++)
        {
            intron_count++;
            gene_vec[i].target_ilist = "INTRON_target_" + IntToString(intron_count);
        }
        
        size_t start_index_custom_promoter = gene_vec.size();
        gene_vec.insert(gene_vec.end(), custom_promoter_vec.begin(), custom_promoter_vec.end()); // add custom promoter interval
        int promoter_count = 0;
        for(size_t i = start_index_custom_promoter; i < gene_vec.size(); i++)
        {
            promoter_count++;
            gene_vec[i].target_ilist = "PROMOTER_target_" + IntToString(promoter_count);
        }
        
        sort(gene_vec.begin(), gene_vec.end(), sort_chrom_numerically);
        for(size_t i = 0; i < gene_vec.size(); i++)
        {
            output_fs << gene_vec[i].chrom << "\t" << gene_vec[i].start << "\t" << gene_vec[i].end << "\t" << "+" << "\t" << gene_vec[i].target_ilist << endl;
        }
        input_fs.close();
        output_fs.close();
    }
    
    void GenerateCombinedBaitIlist(string input_file, string output_file)
    {
        
        cerr << "[INFO]: Loading bait ilist file: " << input_file << endl;
        ifstream input_fs(input_file.c_str());
        if(!input_fs)
        {
            cerr << "[ERROR]: Failed to open bait ilist file: " << input_file << endl;
            exit(1);
        }
        ofstream output_fs(output_file.c_str());
        if(!output_fs)
        {
            cerr << "[ERROR]: Failed to open output bait ilist file: " << input_file << endl;
            exit(1);
        }
        
        vector<Interval> gene_vec; //add custom bait interval
        string line;
        while(getline(input_fs, line)) // add standard bait interval
        {
            if(line[0] == '@')
            {
                output_fs << line << endl;
                continue;
            }
            if(!custom_only)
            {
                vector<string> bait_item;
                split(line, '\t', bait_item);
                Interval new_custom_bait;
                new_custom_bait.chrom = bait_item[0];
                new_custom_bait.start = atoi(bait_item[1].c_str());
                new_custom_bait.end = atoi(bait_item[2].c_str());
                gene_vec.push_back(new_custom_bait);
            }
        }
        
        /*if(!custom_only)
        {
            sort(gene_vec.begin(), gene_vec.end(), compare_interval);
            if( CheckOverlap(custom_bait_vec, gene_vec, true) )
                exit(1);
        }*/
        
        gene_vec.insert(gene_vec.end(), custom_bait_vec.begin(), custom_bait_vec.end());
        sort(gene_vec.begin(), gene_vec.end(), sort_chrom_numerically);
        for(size_t i = 0; i < gene_vec.size(); i++)
        {
            output_fs << gene_vec[i].chrom << "\t" << gene_vec[i].start << "\t" << gene_vec[i].end << "\t" << "+" << "\t" << "tiled_interval_" << setfill('0') << setw(4) << (i+1) << endl;
        }
        input_fs.close();
        output_fs.close();
    }
    
    void GenerateCombinedAminoAcid(string input_file, string output_file1, string output_file2)
    {
        cerr << "[INFO]: Loading canonical AA file: " << input_file << endl;
        ifstream input_fs(input_file.c_str());
        if(!input_fs)
        {
            cerr << "[ERROR]: Failed to open canonical exon with AA file: " << input_file << endl;
            exit(1);
        }
        ofstream output_fs1(output_file1.c_str());
        if(!output_fs1)
        {
            cerr << "[ERROR]: Failed to open output canonical exon with AA file: " << input_file << endl;
            exit(1);
        }
        ofstream output_fs2(output_file2.c_str());
        if(!output_fs2)
        {
            cerr << "[ERROR]: Failed to open output canonical exon file: " << input_file << endl;
            exit(1);
        }
        
        vector<Interval> gene_vec;
        string line;
        if(!custom_only)
        {
            while(getline(input_fs, line)) // add standard bait interval
            {
                vector<string> aa_item;
                split(line, '\t', aa_item);
                Interval new_custom_aa;
                new_custom_aa.chrom = aa_item[0];
                new_custom_aa.start = atoi(aa_item[1].c_str());
                new_custom_aa.end = atoi(aa_item[2].c_str());
                new_custom_aa.strand = aa_item[3];
                new_custom_aa.aa_range = aa_item[4];
                gene_vec.push_back(new_custom_aa);
            }
        }
        
        map<Interval*, int>::iterator it;
        for(it = custom_aa_map.begin(); it != custom_aa_map.end(); it++)
        {
            gene_vec.push_back(*(it->first));
        }
        sort(gene_vec.begin(), gene_vec.end(), sort_chrom_numerically);
        for(size_t i = 0; i < gene_vec.size(); i++)
        {
            output_fs1 << gene_vec[i].chrom << "\t" << gene_vec[i].start << "\t" << gene_vec[i].end << "\t" << gene_vec[i].strand << "\t" << gene_vec[i].aa_range << endl;
            output_fs2 << gene_vec[i].chrom << ":" << gene_vec[i].start << "-" << gene_vec[i].end << endl;
        }
        input_fs.close();
        output_fs1.close();
        output_fs2.close();
    }
    
    void GenerateCombinedTilingInterval(string input_file, string output_file, int left_buffer, int right_buffer, bool output_annotation)
    {
        cerr << "[INFO]: Loading tiling interval file: " << input_file << endl;
        ifstream input_fs(input_file.c_str());
        if(!input_fs)
        {
            cerr << "[ERROR]: Failed to open input tiling interval file: " << input_file << endl;
            exit(1);
        }
        ofstream output_fs(output_file.c_str());
        if(!output_fs)
        {
            cerr << "[ERROR]: Failed to open output tiling interval file: " << output_file << endl;
            exit(1);
        }

        string output_file_chr = output_file + ".chr.list";
        ofstream output_fs_chr(output_file_chr.c_str());
        if(!output_fs_chr)
        {
            cerr << "[ERROR]: Failed to open output tiling interval file: " << output_file_chr << endl;
            exit(1);
        }

        vector<Interval> gene_vec; //add custom tiling interval
        string line;
        if(!custom_only)
        {
            while(getline(input_fs, line)) // add standard tiling interval
            {
                vector<string> gene_item;
                vector<string> chrom_item;
                vector<string> coordinate_item;
                split(line, '\t', gene_item);
                split(gene_item[0], ':', chrom_item);
                split(chrom_item[1], '-', coordinate_item);
                string cytoband = (gene_item.size()) > 1 ? gene_item[1] : "";
                Interval new_custom_tiling;
                new_custom_tiling.chrom = chrom_item[0];
                new_custom_tiling.start = atoi(coordinate_item[0].c_str());
                new_custom_tiling.end = atoi(coordinate_item[1].c_str());
                new_custom_tiling.cytoband = cytoband;
                gene_vec.push_back(new_custom_tiling);
            }
            //sort(gene_vec.begin(), gene_vec.end(), compare_interval);
            //if( CheckOverlap(custom_gene_vec, gene_vec, true) )
            //    exit(1);
        }
        
        for(size_t i = 0; i < custom_tiling_vec.size(); i ++)
        {
            gene_vec.push_back(custom_tiling_vec[i]);
            gene_vec[gene_vec.size() - 1].start = max(1, gene_vec[gene_vec.size() - 1].start - left_buffer);
            gene_vec[gene_vec.size() - 1].end = min(chrom_info_vec[gene_vec[gene_vec.size() - 1].chrom].chrom_len, gene_vec[gene_vec.size() - 1].end + right_buffer);
        }
        sort(gene_vec.begin(), gene_vec.end(), sort_chrom_numerically);
        for(size_t i = 0; i < gene_vec.size(); i++)
        {
            output_fs << gene_vec[i].chrom << ":" << gene_vec[i].start << "-" << gene_vec[i].end;
            if(output_annotation)
                output_fs << "\t" << gene_vec[i].cytoband;
            output_fs << endl;

            output_fs_chr << "chr" << gene_vec[i].chrom << ":" << gene_vec[i].start << "-" << gene_vec[i].end;
            if(output_annotation)
                output_fs_chr << "\t" << gene_vec[i].cytoband;
            output_fs_chr << endl;

        }
        input_fs.close();
        output_fs.close();
        output_fs_chr.close();
    }
    
    void GenerateCombinedFpTilingGenotype(string input_file, string output_file)
    {
        cerr << "[INFO]: Loading fp tiling genotypes file: " << input_file << endl;
        ifstream input_fs(input_file.c_str());
        if(!input_fs)
        {
            cerr << "[ERROR]: Failed to open input fp tiling genotypes file: " << input_file << endl;
            exit(1);
        }
        ofstream output_fs(output_file.c_str());
        if(!output_fs)
        {
            cerr << "[ERROR]: Failed to open output fp tiling genotypes file: " << output_file << endl;
            exit(1);
        }
        
        vector<Interval> gene_vec; //add custom tiling interval
        string line;
        if(!custom_only)
        {
            while(getline(input_fs, line)) // add standard tiling interval
            {
                vector<string> gene_item;
                vector<string> chrom_item;
                split(line, '\t', gene_item);
                split(gene_item[0], ':', chrom_item);
                string genotype = (gene_item.size()) > 1 ? gene_item[1] : "";
                Interval new_custom_tiling;
                new_custom_tiling.chrom = chrom_item[0];
                new_custom_tiling.start = atoi(chrom_item[1].c_str());
                new_custom_tiling.end = atoi(chrom_item[1].c_str());
                new_custom_tiling.genotype = genotype;
                gene_vec.push_back(new_custom_tiling);
            }
            //sort(gene_vec.begin(), gene_vec.end(), compare_interval);
            //if( CheckOverlap(custom_gene_vec, gene_vec, true) )
            //    exit(1);
        }
        
        gene_vec.insert(gene_vec.end(), custom_tiling_vec.begin(), custom_tiling_vec.end());
        sort(gene_vec.begin(), gene_vec.end(), sort_chrom_numerically);
        for(size_t i = 0; i < gene_vec.size(); i++)
        {
            if(!gene_vec[i].genotype.empty())
                output_fs << gene_vec[i].chrom << ":" << gene_vec[i].start << "\t" << gene_vec[i].genotype << endl;
        }
        input_fs.close();
        output_fs.close();
    }
    
    void CreateEmptyFile(string output_file, string header)
    {
        ofstream output_fs(output_file.c_str());
        if(!output_fs)
        {
            cerr << "[ERROR]: Failed to write to finish file: " << output_file << endl;
            exit(1);
        }
        if(!header.empty())
            output_fs << header << endl;
        output_fs.close();
    }
    
    
    vector<Interval> cytoband_vec;
    map<pair<string, string>, string> refseq_canonical;
    map<pair<string, string>, pair<int, int> > gene_longest;
    map<pair<string, string>, pair<int, int> > gene_longest_non_canonical;
    vector<Interval> refseq_exon_vec;
    vector<Interval> refseq_gene_vec;
    vector<Interval> refseq_gene_non_canonical_vec;
    map<pair<string, string>, string> ucsc_dbsnp_map;
    vector<Interval> custom_gene_vec;
    vector<Interval> custom_bait_vec;
    vector<Interval> custom_tiling_vec;
    vector<Interval> custom_intron_vec;
    vector<Interval> custom_promoter_vec;
    map<Interval*, int> custom_aa_map;
    ifstream ref_fs;
    map<string, Chromosome> chrom_info_vec;
    
};



int main(int argc, const char * argv[])
{
    cerr << "[INFO]: Generating custom gene interval files, this may take up to several minutes" << endl;
    parseOption(argc, argv);
    Custom_Design cd;
    cd.OpenReferenceSequence(input_reference);
    cd.LoadCytoband(input_cytoband_file);
    cd.LoadRefseqCanonicalExonAA(input_refseq_canonical_file);
    cd.LoadRefseqGene(input_refseq_file);
    if(!input_custom_bed_file.empty())
    {
        cd.LoadCustomBed(input_custom_bed_file);
    }
    else
    {
        cd.LoadCustomInterval(input_custom_gene_file, cd.custom_gene_vec);
        cd.LoadCustomInterval(input_custom_bait_file, cd.custom_bait_vec);
        if(!input_custom_tiling_file.empty())
        {
            cd.LoadCustomTilingInterval(input_custom_tiling_file, cd.custom_tiling_vec);
        }
    }
    cd.AnnotateCustomGene(output_dir + "/custom_gene_intervals.list.annotated");
    if(!cd.custom_tiling_vec.empty())
    {
        if(input_dbsnp_file.empty())
            printUsage("[ERROR]: Please specify input ucsc dbsnp file");
        if(input_tiling_interval_file.empty())
            printUsage("[ERROR]: Please specify input tiling interval file");
        if(input_tiling_interval_annotated_file.empty())
            printUsage("[ERROR]: Please specify input tiling interval with annotation file");
        if(input_fp_interval_file.empty())
            printUsage("[ERROR]: Please specify input fp tiling interval file");
        if(input_fp_genotype_file.empty())
            printUsage("[ERROR]: Please specify input fp tiling genotype file");
        cd.LoadUcscDBsnp(input_dbsnp_file);
        cd.AnnotateCustomTiling();
    }
    cd.GenerateCombinedGeneInterval(input_gene_interval_annotated_file, output_dir + "/combined_gene_intervals.list");
    //cd.GenerateCombinedGeneInterval(input_gene_interval_annotated_file, output_dir + "/combined_gene_intervals.list.annotated", true);
    cd.GenerateCombinedGeneCoord(input_gene_coord_file, output_dir + "/combined_gene_coords.txt");
    cd.GenerateCombinedGCBias(input_gc_bias_file, output_dir + "/combined_gc_bias.txt", 60, 59);
    cd.GenerateCombinedTargetIlist(input_target_ilist_file, output_dir + "/combined_target.ilist");
    cd.GenerateCombinedBaitIlist(input_bait_ilist_file, output_dir + "/combined_bait.ilist");
    cd.GenerateCombinedAminoAcid(input_aa_file, output_dir + "/combined_canonical_exon_with_aa.list", output_dir + "/combined_canonical_exon.list");
    if(!cd.custom_tiling_vec.empty() || custom_only)
    {
        cd.GenerateCombinedTilingInterval(input_tiling_interval_file, output_dir + "/combined_tiling_interval.list", 60, 59, false);
        cd.GenerateCombinedTilingInterval(input_tiling_interval_annotated_file, output_dir + "/combined_tiling_interval.list.annotated", 60, 59, true);
        cd.GenerateCombinedTilingInterval(input_fp_interval_file, output_dir + "/combined_fp_tiling_interval.list", 60, 58, false);
        cd.GenerateCombinedFpTilingGenotype(input_fp_genotype_file, output_dir + "combined_fp_tiling_genotypes.txt");
    }
    if(custom_only)
    {
        cd.CreateEmptyFile(output_dir + "/Combined_FFPE_CtrlPool_BAM_title.txt", "Barcode\tPool\tSample_ID\tCollab_ID\tPatient_ID\tClass\tSample_type\tInput_ng\tLibrary_yield\tPool_input\tBait_version");
        cd.CreateEmptyFile(output_dir + "/Combined_FFPE_CtrlPool_BAM_ALL_intervalcoverage_loess.txt", "\"Order\"\t\"Target\"\t\"genes\"");
        cd.CreateEmptyFile(output_dir + "/Combined_FFPE_CtrlPool_BAM_ALL_intervalnomapqcoverage_loess.txt", "\"Order\"\t\"Target\"\t\"genes\"");
    }
    cd.CreateEmptyFile(output_dir + "/done.txt", "");
    cerr << "[INFO]: Generating IMPACT+ gene interval files finished" << endl;
    cerr << endl;
    return 0;
}



















