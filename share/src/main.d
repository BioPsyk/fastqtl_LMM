/* The GPL v3 License

   Copyright (C) 2016 University of Geneva.
   #
   # Author: Andrew Brown <andrew.brown@unige.ch>

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

import arg_parse : Opts;
import calculation : genPerms;
import read_data : makeOut, readBed, readEigen, readSquare;
import run_analysis : analyseData;
import std.conv : to;
import std.range : enumerate;
import std.stdio : File, stderr, writeln;

extern (C)
{
  void gslKinship(double* kinship, size_t nInd, double* Sinv, double* eigen);
}

version (GNU)
{
  pragma(msg, "Compiled using gdc");
}
else
{
  pragma(lib, "m");
  pragma(lib, "gsl");
  pragma(lib, "gslcblas");
}

version (unittest)
  void main()
{
  writeln("All unit tests completed successfully.");
}
else
  void main(string[] args)
{
  pragma(msg, "Alfonso");

  const auto opts = new Opts(args.to!(string[]));

  auto phenotype = readBed(opts);

  auto nInd = phenotype[0].values.length;

  auto permutations = genPerms(opts, nInd);

  auto outFile = makeOut(opts);

  double[] Sinv = new double[](nInd * nInd);
  double[] eigen = new double[](nInd);

  if (opts.kinship != "")
  {
    if (opts.verbose)
    {
      stderr.writeln("Calculating Sinv matrix from kinship matrix.");
    }

    double[] kinship = readSquare(opts, nInd);
    if (opts.calc)
    {
      foreach (ref e; phenotype)
      {
        version (EMBEDR)
        {
          import call_r_heritability : h2CalcEmbedR;

          e.h2 = h2CalcEmbedR(e.values, kinship, e.geneName, opts.verbose);
        }
        else
        {
          import heritability : h2Calc;

          e.h2 = h2Calc(e.values, kinship, e.geneName, opts.verbose);
        }
      }
    }

    gslKinship(kinship.ptr, nInd, Sinv.ptr, eigen.ptr);

  }
  else
  {
    if (opts.verbose)
    {
      stderr.writeln("Reading Sinv matrix and eigenvalues");
    }

    Sinv = readSquare(opts, nInd);
    eigen = readEigen(opts, nInd);
  }

  foreach (ref e; phenotype.enumerate)
  {
    if (opts.verbose)
    {
      stderr.writeln("Analysing gene ", e[1].geneName, " (", e[0] + 1,
          " out of ", phenotype.length, ").");
    }
    analyseData(e[1], permutations, outFile, opts, Sinv, eigen);
  }
}

@system unittest
{
  import core.stdc.stdlib : exit;
  import std.digest.sha : SHA1, toHexString;
  import std.file : exists, remove;
  import std.range : put;
  import std.uuid : randomUUID;

  auto testFile = randomUUID.toString;

  while (testFile.exists)
  {
    testFile = randomUUID.toString;
  }

  scope (exit)
  {
    if (testFile.exists)
      testFile.remove;
  }

  // ./bin/alfonso --bed data/phenotype.bed --job-number 1 --genes 10 --vcf data/genotype.vcf.gz --perm 100000,4 --Sinv data/Sinv --eigen data/eigen
  immutable auto args = [
    "prog", "--bed", "data/phenotype.bed", "--job-number", "1", "--genes", "10", "--vcf", "data/genotype.vcf.gz",
    "--perm", "100000,4", "--Sinv", "data/Sinv", "--eigen", "data/eigen", "--out", testFile
  ];
  auto opts = new Opts(args.to!(string[]));

  auto phenotype = readBed(opts);
  auto nInd = phenotype[0].values.length;
  auto permutations = genPerms(opts, nInd);

  auto outFile = makeOut(opts);

  double[] Sinv = readSquare(opts, nInd);
  double[] eigen = readEigen(opts, nInd);

  foreach (ref e; phenotype)
  {
    analyseData(e, permutations, outFile, opts, Sinv, eigen);
  }

  outFile.close;

  SHA1 hash;
  hash.start;
  put(hash, File(testFile).byChunk(1024));

  assert(toHexString(hash.finish) == "80745CCF282B5069B75B99470F26DB8466B15ABE");

  stderr.writeln("Passed test with Sinv and eigen.");

  // ./bin/alfonso --bed data/phenotype.bed --job-number 1 --genes 10 --vcf data/genotype.vcf.gz --perm 100000,4 --kinship data/kinship

  opts.kinship = "data/kinship";

  double[] kinship = readSquare(opts, nInd);

  gslKinship(kinship.ptr, nInd, Sinv.ptr, eigen.ptr);

  outFile = makeOut(opts);

  foreach (ref e; phenotype)
  {
    analyseData(e, permutations, outFile, opts, Sinv, eigen);
  }

  outFile.close;

  hash.start;
  put(hash, File(testFile).byChunk(1024));

  assert(toHexString(hash.finish) == "0E50C4075C347F5E596EE120A68B7154C64625C2");

  stderr.writeln("Passed test with kinship matrix.");

  // ./bin/alfonso --bed data/phenotype.bed --job-number 1 --genes 10 --vcf data/genotype.vcf.gz --perm 100000,4 --kinship data/kinship --calc-h2

  opts.calc = true;

  phenotype = readBed(opts);
  kinship = readSquare(opts, nInd);

  foreach (ref e; phenotype)
  {
    version (EMBEDR)
    {
      import call_r_heritability : h2CalcEmbedR;

      e.h2 = h2CalcEmbedR(e.values, kinship, e.geneName, opts.verbose);
    }
    else
    {
      import heritability : h2Calc;

      e.h2 = h2Calc(e.values, kinship, e.geneName, opts.verbose);
    }
  }

  gslKinship(kinship.ptr, nInd, Sinv.ptr, eigen.ptr);

  outFile = makeOut(opts);

  foreach (ref e; phenotype)
  {
    analyseData(e, permutations, outFile, opts, Sinv, eigen);
  }

  outFile.close;

  hash.start;
  put(hash, File(testFile).byChunk(1024));

  //  assert(toHexString(hash.finish) == "C04C005DC615A9212F20D8BD304731626D9FDD61"); //heritability estimated using regress
  // assert(toHexString(hash.finish) == "F2C6522C50D9084D26ABD865C8AA3285DF9A1304"); //heritability estimated using genabel
  assert(toHexString(hash.finish) == "A6D2F414FC237BAB7D6E73D1C98D634EDACC9C87"); //heritability estimated using genabel_hglm

  stderr.writeln("Passed test with heritability calculation.");

}
