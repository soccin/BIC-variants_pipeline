//
//  TrimProblematicQuality.cpp
//  TrimProblematicQuality
//
//  Created by Zeng, Zheng/Sloan-Kettering Institute on 1/6/14.
//  Copyright (c) 2014 Zeng, Zheng/Sloan-Kettering Institute. All rights reserved.
//
//

#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <string.h>
#include <getopt.h>
#include <stdlib.h>

//#define _DEBUG

using namespace std;


const string version = "TrimProblematicQuality 1.0.1";

string input_file;
string output_file;
bool checking_only = false;
int check_num_read = 0;
int quality_threshold = 999;
int encoding = 33;
int min_read_len = 50;

int total_reads_processed = 0;
int total_reads_trimmed = 0;
int total_trimmed_base = 0;


void printUsage(string msg = "")
{
    cout << endl;
    cout << version << endl;
    cout << "This program checks if the quality score of input Fastq file is exceed the threshold and trim it" << endl;
    cout << "Usage: " << endl;
    cout << "[REQUIRED ARGUMENTS]" << endl;
    cout << "\t--input              <string>                        Input Fastq file" << endl;
    cout << "\t--output             <string>                        Output Fastq file" << endl;
    cout << "\t--quality            <int>                           Quality threshold, any base with quality score higher than this and all the following bases will be trimmed" << endl;
    cout << endl;
    
    
    cout << "[OPTIONAL ARGUMENTS]" << endl;
    cout << "\t--encoding           <int>                           Quality score encoding. Default [33] for Phred+33. Change this to 64 for Solexa+64 and Phred+64" << endl;
    cout << "\t--min_len           <int>                            Reads that are shorter than this after trimming will be discarded. Default [50]" << endl;
    cout << "\t--num_read           <int>                           Only use the first <int> reads. By default all reads will be used" << endl;
    cout << "\t--checking_only                                      Checking only, no trimming" << endl;
    if(!msg.empty())
        cerr << msg << endl;
    exit(1);
}

static struct option long_options[] =
{
    {"input",                   required_argument,      0,     'i'},
    {"output",                  required_argument,      0,     'o'},
    {"quality",                 required_argument,      0,     'q'},
    {"encoding",                required_argument,      0,     'e'},
    {"min_len",                 required_argument,      0,     'm'},
    {"num_read",                required_argument,      0,     'n'},
    {"checking_only",           no_argument,            0,     'c'},
    {0, 0, 0, 0}
};



void parseOption(int argc, char* argv[])
{
    if(argc == 1)
        printUsage();
    int next_option;
    int option_index = 0;
    do
    {
        next_option = getopt_long(argc, const_cast<char**>(argv), "i:o:q:e:m:n:c", long_options, &option_index);
        switch(next_option)
        {
            case 'i':
                input_file = optarg;
                break;
            case 'o':
                output_file = optarg;
                break;
            case 'q':
                quality_threshold = atoi(optarg);
                if(quality_threshold < -10 || quality_threshold > 100)
                    printUsage("--quality must between [-10, 100]");
                break;
            case 'e':
                encoding = atoi(optarg);
                if(encoding != 33 || encoding != 64)
                    printUsage("--encoding must be 33 or 64");
                break;
            case 'm':
                min_read_len = atoi(optarg);
                if(min_read_len <= 0)
                    printUsage("--min_len must > 0");
                break;
            case 'n':
                check_num_read = atoi(optarg);
                if(check_num_read <= 0)
                    printUsage("--check_num_read must > 0");
                break;
            case 'c':
                checking_only = true;
                break;
            case -1:
                break; //parsed all options
            default:
                printUsage(string("[ERROR]: Argument error: ") + argv[optind - 1]);
        }
    } while(next_option != -1);
    if(input_file.empty())
        printUsage("[ERROR]: Please specify input file");
    if(!checking_only && output_file.empty())
        printUsage("[ERROR]: Please specify output file");
    if(quality_threshold == 999)
        printUsage("[ERROR]: Please specify quality threshold");
    quality_threshold += encoding;
#ifdef _DEBUG
        cout << "[DEBUG]: Parsing options complete." << endl;
#endif
}


void trim_quality(const char* input_filename, const char* output_filename)
{
    cout << "Checking quality for file: " << input_filename << endl;
    ifstream input_fs(input_filename);
    if(!input_fs)
    {
        cerr << "open input file error: " << input_filename << endl;
        exit(1);
    }
    ofstream output_fs;
    if(!checking_only)
    {
        output_fs.open(output_filename, ofstream::out);
        if(!output_fs)
        {
            cerr << "open output file error: " << output_filename << endl;
            exit(1);
        }
    }
    
    string line1, line2, line3, line4;
    //bool report_warning = true;
    while(true)
    {
        getline(input_fs, line1);
        getline(input_fs, line2);
        getline(input_fs, line3);
        getline(input_fs, line4);
        if(input_fs.eof() || (check_num_read != 0 && total_reads_processed >= check_num_read))
            break;
        total_reads_processed ++;
        
        size_t trim_index;
        for (trim_index = 0; trim_index < line4.size(); trim_index ++)
        {
            char qual = line4[trim_index];
            if(qual > quality_threshold)
            {
#ifdef _DEBUG
                //if(report_warning)
                {
                    cerr << "[DEBUG]: warning, quality score at position " << (trim_index+1) << " above debug threhold: " << int(qual) << ">" << quality_threshold << endl;
                    cerr << "[DEBUG]: " << line1 << endl;
                    cerr << "[DEBUG]: " << line2 << endl;
                    cerr << "[DEBUG]: " << line3 << endl;
                    cerr << "[DEBUG]: " << line4 << endl;
                    cerr << "[DEBUG]: " << endl;
                    //report_warning = false;
                }
#endif
                break;
            }
        }
        size_t kept_len = trim_index;
        size_t trimmed_len = line4.size() - trim_index;
        if(trimmed_len > 0)
        {
                total_reads_trimmed ++;
                total_trimmed_base += trimmed_len;
        }
        if(!checking_only && kept_len >= min_read_len)
        {
            output_fs << line1 << endl;
            output_fs << line2.substr(0, kept_len) << endl;
            output_fs << line3 << endl;
            output_fs << line4.substr(0, kept_len) << endl;
        }
    }
    double average_trimmed_base = 0;
    if(total_reads_trimmed != 0)
        average_trimmed_base = ((double)total_trimmed_base / (double)total_reads_trimmed);
    
    cout << "total reads processed: " << total_reads_processed << endl;
    cout << "total reads trimmed: " << total_reads_trimmed << "(" << ((double)total_reads_trimmed / (double)total_reads_processed * 100 ) << "%)" << endl;
    cout << "average trimmed bases per trimmed read: " << average_trimmed_base << endl;
    cout << endl;
            
    input_fs.close();
    output_fs.close();
}

int main(int argc, char * argv[])
{
    parseOption(argc, argv);
    trim_quality(input_file.c_str(), output_file.c_str());
    return 0;
}


