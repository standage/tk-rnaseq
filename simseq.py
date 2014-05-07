#!/usr/bin/env python

# Copyright (c) 2014, Daniel S. Standage <daniel.standage@gmail.com>
#
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

import getopt
import os
import re
import shutil
import subprocess
import sys

class Molecule:
  """
  Represents the relative within-sample abundance of the given molecule
  """
  def __init__(self, seqid, abundance):
    assert abundance >= 0
    self.seqid = seqid
    self.abundance = int(abundance)

  def __str__(self):
    return "[%s:%d]" % (self.seqid, self.abundance)


class Sample:
  """
  Represents a sample of RNA, with a baseline abundance and a list of molecules
  with their relative abundances
  """
  def __init__(self, sampleid):
    self.id = sampleid
    self.baseline = 100
    self.molecules = []

  def __str__(self):
    return "[%s:%d %d]" % (self.id, self.baseline, len(self.molecules))

def id_safe(seqid):
  """
  Sequence IDs and sample IDs are both used in creating intermediate files. This
  function is used to screen for characters in IDs that may cause problems.
  """
  if re.search("[^A-Za-z0-9._-]", seqid):
    return False
  return True

def parse_fasta(fp):
  """
  Stolen shamelessly from http://stackoverflow.com/a/7655072/459780.
  """
  name, seq = None, []
  for line in fp:
    line = line.rstrip()
    if line.startswith(">"):
      if name: yield (name, ''.join(seq))
      name, seq = line, []
    else:
      seq.append(line)
  if name: yield (name, ''.join(seq))

def split_sequences(instream, tmpdir):
  if not os.path.exists(tmpdir):
    os.makedirs(tmpdir)
  with instream as fp:
    for defline, seq in parse_fasta(fp):
      seqid = re.match(">(\S+)", defline).group(1)
      assert id_safe(seqid)
      tmpfile = "%s/%s.fa" % (tmpdir, seqid)
      outstream = open(tmpfile, 'w')
      print >> outstream, "%s\n%s" % (defline, seq)
      outstream.close()

def print_usage(outstream):
  print >> outstream, (
"\nSimSeq: simulate a paired-end RNA-seq experiment\n"+
"Usage: python simseq.py [options] sampling-config\n\n"+
"'sampling-config' refers to a list of 3-tuples specifying a sample, a sequence\n"+
"ID, and the relative abundance of that sequence in the sample. This is\n"+
"demonstrated in the following example.\n\n"+
"  python simseq.py samp1,seq1,1 samp1,seq2,4 samp2,seq1,1 samp2,seq1,1 < seqs.fa\n\n"+
"Options:\n"+
"  -b|--baseline: STR,INT    define the baseline per-molecule abundance for the\n"+
"                            given sample, formatted like 'sample1,10000'; if a\n"+
"                            sample's baseline abundance is not explicitly\n"+
"                            stated, it defaults to 100\n"+
"  -c|--config-file: FILE    provide 'sampling-config' through a file rather than\n"+
"                            through command-line arguments; each line of the\n"+
"                            file should contain a single 3-tuple\n"+
"  -h|--help                 print this help message and exit\n"+
"  -i|--input: FILE          specify Fasta file containing sequences to be\n"+
"                            sequenced in silico; default is standard input\n"+
"  -t|--tmpdir: DIR          specify a directory for placing intermediate files;\n"+
"                            default is /tmpdir/simseq-$pid\n"+
"  -w|--wgsim-opts: STR      options to be passed to wgsim; default is\n"+
"                            '-d 270 -1 100 -2 100'\n")

def load_config(argv):
  baseline_args = []
  sampling_args = []
  instream = sys.stdin
  params = {
    "tmpdir": "/tmp/simseq-%d" % os.getpid(), 
    "wgsim-opts": "-d 270 -1 100 -2 100",
  }

  optstr = "b:c:hi:p:t:w:"
  longopts = ["baseline=", "config-file=", "help", "input=", "tmpdir=", "wgsim-opts="]
  (options, args) = getopt.getopt(argv[1:], optstr, longopts)
  for key, value in options:
    if key in ("-b", "--baseline"):
      baseline_args.append(value)
    elif key in ("-c", "--config-file"):
      cfg = open(value, 'r')
      for line in cfg:
        line = line.rstrip()
        if line == "" or line.startswith("#"):
          continue
        elif line.startswith("--baseline="):
          baseline_args.append(line[11:])
        else:
          sampling_args.append(line)
    elif key in ("-h", "--help"):
      print_usage(sys.stdout)
      exit(0)
    elif key in ("-i", "--input"):
      instream = open(value, 'r')
    elif key in ("-t", "--tmpdir"):
      params["tmpdir"] = value
    elif key in ("-w", "--wgsim-opts"):
      params["wgsim-opts"] = value
    else:
      print >> sys.stderr, "error: unknown option '%s'" % key
      exit(1)
  for arg in args:
    sampling_args.append(arg)  

  split_sequences(instream, params["tmpdir"])
  instream.close()

  baselines = {}
  for arg in baseline_args:
    sampleid, baseline = arg.split(",")
    baselines[sampleid] = baseline

  samples = {}
  for arg in sampling_args:
    sampleid, seqid, seqabund = arg.split(",")
    assert id_safe(sampleid) and id_safe(seqid)
    if sampleid not in samples:
      sample = Sample(sampleid)
      if sampleid in baselines:
        sample.baseline = int(baselines[sampleid])
      samples[sample.id] = sample
    mol = Molecule(seqid, int(seqabund))
    samples[sampleid].molecules.append(mol)

  return params, samples

if __name__ == "__main__":
  params, samples = load_config(sys.argv)
  for sampleid in samples:
    sample = samples[sampleid]
    for mol in sample.molecules:
      abundance = sample.baseline * mol.abundance
      cmd = "wgsim %s -N %d %s/%s.fa %s/%s.%s.1.fq %s/%s.%s.2.fq" % (params["wgsim-opts"], abundance, params["tmpdir"], mol.seqid, params["tmpdir"], sampleid, mol.seqid, params["tmpdir"], sampleid, mol.seqid)
      subprocess.call(cmd, shell=True)

  for sampleid in samples:
    subprocess.call("cat %s/%s.*.1.fq > %s.1.fq" % (params["tmpdir"], sampleid, sampleid), shell=True)
    subprocess.call("cat %s/%s.*.2.fq > %s.2.fq" % (params["tmpdir"], sampleid, sampleid), shell=True)
    
  shutil.rmtree(params["tmpdir"])

