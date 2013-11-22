/*

Copyright (c) 2013, Daniel S. Standage <daniel.standage@gmail.com>

Permission to use, copy, modify, and/or distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.


SMR: SAM mapped reads

The SAM file format encodes the alignment (mapping) of short sequence reads to
longer molecular sequences. This program reads each entry in the SAM file to
compute a tally of the number of short reads mapped to each molecule.

*/

#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>


/**
 * @type SmrOptions
 *
 * Container and parser for handling command-line options and arguments.
 */
typedef struct SmrOptions SmrOptions;
struct SmrOptions
{
  char delim;
  const char *outfile;
  FILE *outstream;
  unsigned numfiles;
  std::vector<const char *> infiles;

  SmrOptions(int argc, char **argv)
  {
    delim = ',';
    outfile = "stdout";

    char opt;
    const char *arg;
    while((opt = getopt(argc, argv, "d:ho:")) != -1)
    {
      switch(opt)
      {
        case 'd':
          arg = optarg;
          if(strcmp(optarg, "\\t") == 0)
            arg = "\t";
          else if(strlen(optarg) > 1)
          {
            fprintf(stderr, "warning: string '%s' provided for delimiter, using "
                    "only '%c'\n", optarg, optarg[0]);
          }
          delim = arg[0];
          break;
        case 'h':
          usage(stderr);
          exit(0);
          break;
        case 'o':
          outfile = optarg;
          break;
        default:
          fprintf(stderr, "error: unknown option '%c'\n", opt);
          usage(stderr);
          break;
      }
    }

    outstream = stdout;
    if(strcmp(outfile, "stdout") != 0)
    {
      outstream = fopen(outfile, "w");
      if(outstream == NULL)
      {
        fprintf(stderr, "error opening file %s", outfile);
        exit(1);
      }
    }

    numfiles = argc - optind;
    if(numfiles < 1)
    {
      fprintf(stderr, "error: expected 1 or more input files\n");
      usage(stderr);
      exit(1);
    }

    for(unsigned i = 0; i < numfiles; i++)
      infiles.push_back(argv[optind+i]);
  }

  ~SmrOptions()
  {
    fclose(outstream);
  }

  void usage(FILE *outstream)
  {
    fprintf(outstream, "\nSMR: SAM mapped reads\n\n"
"The input to SMR is 1 or more SAM files. The output is a table (1 column for\n"
"each input file) showing the number of reads that map to each molecule.\n\n"
"Usage: smr [options] sample-1.sam sample-2.sam ... sample-n.sam\n"
"  Options:\n"
"    -d|--delim: CHAR         delimiter for output data; default is comma\n"
"    -h|--help                print this help message and exit\n"
"    -o|--outfile: FILE       name of file to which read counts will be\n"
"                             written; default is terminal (stdout)\n\n");
  }
};


/**
 * @type ReadTally
 *
 * This class is an instance of an unordered map. Each key is a unique ID
 * corresponding to a molecule, and the value is the number of reads mapped to
 * that molecule.
 */
#define MAX_LINE_LENGTH 2048
typedef struct ReadTally ReadTally;
struct ReadTally : public std::unordered_map<std::string, unsigned>
{
  FILE *instream;
  ReadTally(const char *infilename)
  {
    instream = fopen(infilename, "r");
    if(instream == NULL)
    {
      fprintf(stderr, "error opening file %s\n", infilename);
      exit(1);
    }

    char buffer[MAX_LINE_LENGTH];
    while(fgets(buffer, MAX_LINE_LENGTH, instream) != NULL)
    {
      if(buffer[0] == '@')
        continue;

      char *tok = strtok(buffer, "\t\n");
      tok = strtok(NULL, "\t\n");
      int bflag = atoi(tok);
      if(bflag & 0x4)
        continue;

      tok = strtok(NULL, "\t\n");
      std::string molid (tok);

      ReadTally::iterator kvpair = this->find(molid);
      if(kvpair == this->end())
        this->emplace(molid, 1);
      else
        kvpair->second += 1;
    }
  }

  ~ReadTally() { fclose(instream); }
};


/**
 * @class ReadTallyMatrix
 *
 * This class is an instance of a vector, where the vector elements are
 * ReadTally objects (described above). Each row in the matrix corresponds to a
 * molecule, and each column corresponds to one of the input files. The order of
 * the columns is the same as the order of the input files.
 */
typedef struct ReadTallyMatrix ReadTallyMatrix;
struct ReadTallyMatrix : public std::vector<ReadTally>
{
  ReadTallyMatrix(std::vector<const char *>& infiles)
  {
    for(auto& infilename : infiles)
      this->emplace_back(infilename);
  }

  void print(FILE *outstream, char delim)
  {
    std::unordered_set<std::string> molids;
    for(auto& readTally : *this)
    {
      for(auto& kvpair : readTally)
        molids.emplace(kvpair.first);
    }
    
    for(auto& molid : molids)
    {
      fprintf(outstream, "%s%c", molid.c_str(), delim);
    
      bool printdelim = false;
      for(auto& readTally : *this)
      {
        if(printdelim)
          fputc(delim, outstream);
        else
          printdelim = true;
    
        ReadTally::const_iterator kvpair = readTally.find(molid);
        if(kvpair == readTally.end())
          fputc('0', outstream);
        else
          fprintf(outstream, "%u", kvpair->second);
      }
      fprintf(outstream, "\n");
    }
  }
};


// Main method
int main(int argc, char **argv)
{
  SmrOptions options(argc, argv);
  ReadTallyMatrix readTalliesPerSample(options.infiles);
  readTalliesPerSample.print(options.outstream, options.delim);
  return 0;
}
