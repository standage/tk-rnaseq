/*

Copyright (c) 2014, Daniel S. Standage <daniel.standage@gmail.com>

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


--------------------------------------------------------------------------------
lenpick: pick the minimum peak length from a set of read length distributions

Some applications of RNA-seq require reads to be exactly the same length.
However, adapter and quality trimming result in read sets with uneven lengths.
This program was written to find the peak length in each read set's length
distribution and then report the minimum peak length for all read sets. This
number serves to guide additional post-QC trimming to a uniform read length.

Compile like so:    g++ -Wall -O3 -std=c++11 -o lenpick lenpick.cpp
--------------------------------------------------------------------------------

*/

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

int main (int argc, char **argv)
{
  std::vector<unsigned> peaks;

  bool noarg = argc == 1;
  bool hflag = argc > 1 &&
               (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0);
  if (noarg || hflag)
  {
    puts("Usage: lenpick seq1.fq [seq2.fq seq3.fq ...]");
    return 0;
  }

  for (int i = 1; i < argc; i++)
  {
    std::time_t starttime = std::time(NULL);
    std::unordered_map<unsigned, unsigned> lendist;
    std::string buffer;
    std::ifstream instream (argv[i]);

    unsigned count = 0;
    while (std::getline(instream, buffer))
    {
      if(count++ % 4 == 1)
      {
        unsigned length = buffer.length();
        auto kvpair = lendist.find(length);
        if (kvpair == lendist.end())
          lendist.emplace(length, 1);
        else
          kvpair->second += 1;
      }
    }
    instream.close();

    auto lengthpeak = lendist.end();
    for (auto iter = lendist.begin(); iter != lendist.end(); iter++)
    {
      if (lengthpeak == lendist.end() || iter->second > lengthpeak->second)
        lengthpeak = iter;
    }
    peaks.push_back(lengthpeak->first);

    std::time_t elapsedtime = std::time(NULL) - starttime;
    fprintf(stderr, "[lenpick] Peak length for '%s': %u (%ld seconds)\n",
            argv[i], lengthpeak->first, elapsedtime);
  }

  auto minpeak = std::min_element(std::begin(peaks), std::end(peaks));
  printf("%u\n", *minpeak);
}

