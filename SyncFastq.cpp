//
//  SyncFastq.cpp
//  SyncFastq
//
//  Created by Zeng, Zheng/Sloan-Kettering Institute on 1/7/14.
//  Copyright (c) 2014 Zeng, Zheng/Sloan-Kettering Institute. All rights reserved.
//

#include <string.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <map>
#include <stdlib.h>
#include <getopt.h>


using namespace std;
#pragma warning(disable:4996)

const string version = "SyncFastq 1.0.0";

string input_file1;
string input_file2;
string output_file1;
string output_file2;
string output_file1_single;
string output_file2_single;
string file_format;
bool fastq_format = false;




void printUsage(string msg = "")
{
    cout << endl;
    cout << version << endl;
    cout << "This program sychronize the two ends of pair end read FASTA/FASTQ files" << endl;
    cout << "Usage: " << endl;
    cout << "[REQUIRED ARGUMENTS]" << endl;
    cout << "\t--s1                 <string>                        Input file for end1" << endl;
    cout << "\t--s2                 <string>                        Input file for end2" << endl;
    cout << "\t--t1                 <string>                        Output file for end1" << endl;
    cout << "\t--t2                 <string>                        Output file for end2" << endl;
    cout << "\t--f                  <string>                        Input file format. FASTA or FASTQ" << endl;
    cout << endl;
    if(!msg.empty())
        cerr << msg << endl;
    exit(1);
}

static struct option long_options[] =
{
    {"s1",                   required_argument,      0,     'i'},
    {"s2",                   required_argument,      0,     'I'},
    {"t1",                   required_argument,      0,     't'},
    {"t2",                   required_argument,      0,     'T'},
    {"f",                    required_argument,      0,     'f'},
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
        next_option = getopt_long(argc, const_cast<char**>(argv), "i:I:t:T:f:", long_options, &option_index);
        switch(next_option)
        {
            case 'i':
                input_file1 = optarg;
                break;
            case 'I':
                input_file2 = optarg;
                break;
            case 't':
                output_file1 = optarg;
                output_file1_single = output_file1 + ".single";
                break;
            case 'T':
                output_file2 = optarg;
                output_file2_single = output_file2 + ".single";
                break;
            case 'f':
                file_format = optarg;
                if(file_format != "FASTA" && file_format != "FASTQ")
                {
                    printUsage("unrecognized file format");
                    exit(1);
                }
                break;
            case -1:
                break; //parsed all options
            default:
                printUsage(string("[ERROR]: Argument error: ") + argv[optind - 1]);
        }
    } while(next_option != -1);
    if(input_file1.empty())
        printUsage("[ERROR]: Please specify input file1");
    if(input_file2.empty())
        printUsage("[ERROR]: Please specify input file2");
    if(input_file1 == input_file2)
        printUsage("[ERROR]: input file1 and input file2 are the same");
    if(output_file1.empty())
        printUsage("[ERROR]: Please specify output file1");
    if(output_file2.empty())
        printUsage("[ERROR]: Please specify output file2");
    if(output_file1 == output_file2)
        printUsage("[ERROR]: output file1 and output file2 are the same");
    if(file_format.empty())
        printUsage("[ERROR]: Please specify input file format");
    if(file_format == "FASTQ")
        fastq_format = true;
#ifdef _DEBUG
    cout << "[DEBUG]: Parsing options complete." << endl;
#endif
}


void sync_fastq()
{

	ifstream end1_input_fs(input_file1.c_str());
    if(!end1_input_fs)
    {
        cerr << "unable to open input file: " << input_file1 << endl;
        exit(1);
    }
	ifstream end2_input_fs(input_file2.c_str());
    if(!end2_input_fs)
    {
        cerr << "unable to open input file: " << input_file2 << endl;
        exit(1);
    }
    ofstream end1_paired_fs(output_file1.c_str());
    if(!end1_paired_fs)
    {
        cerr << "unable to open output file: " << output_file1 << endl;
        exit(1);
    }
 	ofstream end2_paired_fs(output_file2.c_str());
    if(!end2_paired_fs)
    {
        cerr << "unable to open output file: " << output_file2 << endl;
        exit(1);
    }
	ofstream end1_single_fs(output_file1_single.c_str());
    if(!end1_single_fs)
    {
        cerr << "unable to open output file: " << output_file1_single << endl;
        exit(1);
    }
	ofstream end2_single_fs(output_file2_single.c_str());
    if(!end2_single_fs)
    {
        cerr << "unable to open output file: " << output_file2_single << endl;
        exit(1);
    }
    map<string, int> read_table;
    string line;
    char name_tmp[1000];
    ofstream* current_fs;
	while(getline(end1_input_fs, line))
	{
		sscanf(line.c_str(), "%*c%s", name_tmp);
		string read_name = name_tmp;
        if(read_name.length() > 2 && read_name[read_name.size() - 2] == '/')
            read_name = read_name.substr(0, read_name.length() - 2);
        read_table.insert(make_pair(read_name, 1));
		getline(end1_input_fs, line);
        if(fastq_format)
        {
            getline(end1_input_fs, line);
            getline(end1_input_fs, line);
        }
	}
	while(getline(end2_input_fs, line))
	{
		sscanf(line.c_str(), "%*c%s", name_tmp);
		string read_name = name_tmp;
        if(read_name.length() > 2 && read_name[read_name.size() - 2] == '/')
            read_name = read_name.substr(0, read_name.length() - 2);
		if(read_table.find(read_name) != read_table.end())
		{
			current_fs = &end2_paired_fs;
			read_table[read_name] ++;
		}
		else
			current_fs = &end2_single_fs;
		(*current_fs) << line << endl;
		getline(end2_input_fs, line);
		(*current_fs) << line << endl;
		if(fastq_format)
        {
            getline(end2_input_fs, line);
            (*current_fs) << line << endl;
            getline(end2_input_fs, line);
            (*current_fs) << line << endl;
        }
	}
	end1_input_fs.clear();
	end1_input_fs.seekg(0, ios_base::beg);
	while(getline(end1_input_fs, line))
	{
		sscanf(line.c_str(), "%*c%s", name_tmp);
		string read_name = name_tmp;
        if(read_name.length() > 2 && read_name[read_name.size() - 2] == '/')
            read_name = read_name.substr(0, read_name.length() - 2);
		if(read_table[read_name] == 2)
			current_fs = &end1_paired_fs;
		else
			current_fs = &end1_single_fs;
		(*current_fs) << line << endl;
		getline(end1_input_fs, line);
		(*current_fs) << line << endl;
        if(fastq_format)
        {
            getline(end1_input_fs, line);
            (*current_fs) << line << endl;
            getline(end1_input_fs, line);
            (*current_fs) << line << endl;
        }
	}
	end1_input_fs.close();
	end2_input_fs.close();
	end1_paired_fs.close();
	end1_single_fs.close();
	end2_paired_fs.close();
	end2_single_fs.close();
}


int main(int argc, const char * argv[])
{
    parseOption(argc, argv);
    sync_fastq();
    return 0;
}

