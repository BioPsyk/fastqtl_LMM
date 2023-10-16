module read_data;

import arg_parse : Opts;
import core.exception : RangeError;
import core.stdc.stdlib : exit;
import std.algorithm : countUntil, map, max;
import std.array : array, split;
import std.conv : to;
import std.exception : enforce;
import std.format : format;
import std.process : pipeShell, Redirect, wait;
import std.range : indexed, iota;
import std.stdio : File, readln, stderr, stdout, writeln;

class InputException : Exception
{
  pure this(string s)
  {
    super(s);
  }
}

struct Phenotype
{
  string geneName;
  string chromosome;
  double h2;
  size_t location;
  double[] values;

  this(char[] line, const size_t[] indices)
  {
    auto splitLine = line.split;
    geneName = splitLine[3].to!string;
    h2 = splitLine[1].to!double;
    chromosome = splitLine[0].to!string;
    location = splitLine[2].to!size_t;
    values = splitLine[4 .. $].indexed(indices).array.to!(double[]);
    if (countUntil!"a != b"(values, values[0]) == -1)
      throw new InputException("");
  }
}

struct Genotype
{
  string snpId;
  double[] values;
  double cor;

  this(char[] line, const size_t[] indices, const long loc, const bool gt)
  {
    auto splitLine = line.split;
    snpId = format("\t%-(%s\t%)\t", splitLine[0 .. 4]);
    values = splitLine[4 .. $].indexed(indices).map!(a => getDosage(a, loc, gt)).array;
    if (countUntil!"a != b"(values, values[0]) == -1)
      throw new InputException("");
    cor = 0;
  }
}

double getDosage(char[] field, long loc, bool gt)
{
  auto fieldSplit = field.split(':');
  enforce(fieldSplit.length > loc, new InputException(""));
  return gt ? cast(ubyte) fieldSplit[loc][0] + cast(ubyte) fieldSplit[loc][2] - 96
    : fieldSplit[loc].to!double;
}

auto readBed(const Opts opts)
{
  File bedFile;
  try
  {
    bedFile = File(opts.bed);
  }
  catch (Exception e)
  {
    stderr.writeln(e.msg);
    exit(1);
  }

  bedFile.readln;

  foreach (i; iota((opts.jobNumber - 1) * opts.genes))
  {
    try
    {
      bedFile.readln;
    }
    catch (Exception e)
    {
      stderr.writeln("Too few phenotypes in bed file.");
      exit(1);
    }
  }

  if (bedFile.eof)
  {
    stderr.writeln("Too few phenotypes in bed file.");
    exit(1);
  }

  Phenotype[] phenotype;

  size_t geneCount = opts.genes == 0 ? uint.max : opts.genes;

  foreach (line; bedFile.byLine)
  {
    if (geneCount == 0)
      break;

    try
    {
      phenotype ~= Phenotype(line, opts.phenotypeLocations);
    }
    catch (Exception e)
    {
    }
    geneCount--;
  }

  if (phenotype.length == 0)
  {
    stderr.writeln("No phenotypes read from file.");
    exit(1);
  }

  return phenotype;
}

auto readGenotype(const Opts opts, string chrom, size_t location, string geneName)
{

  immutable auto start = location < opts.window ? 0 : location - opts.window;

  auto tabixCommand = "tabix " ~ opts.vcf ~ " " ~ chrom ~ ":" ~ start.to!string ~ "-" ~ (
      location + opts.window).to!string ~ " | grep -v '#' | cut -f1,2,4,5,10-";

  auto pipes = pipeShell(tabixCommand, Redirect.stdout);
  scope (exit)
    wait(pipes.pid);

  Genotype[] genotype;

  foreach (line; pipes.stdout.byLine)
  {
    try
    {
      genotype ~= Genotype(line, opts.genotypeLocations, opts.loc, opts.gt);
    }
    catch (Exception e)
    {
    }
  }

  if (genotype.length == 0)
  {
    stderr.writeln("Failed to extract any useful SNPs for ", geneName,
        ". Run: \n\n\"", tabixCommand, "\"\n\nfor more information.\n");
  }

  return genotype;
}

auto makeOut(const Opts opts)
{
  File outFile;

  try
  {
    if (opts.output == "")
      outFile = stdout;
    else
      outFile = File(opts.output, "w");

    if (!opts.noheader)
      outFile.writeln("GENE\tCHROM\tPOS\tREF\tALT\tCOR\tP\tPerm\tFWER\tBETA");
  }
  catch (Exception e)
  {
    stderr.writeln("Failed to write to output file.");
    exit(1);
  }
  return outFile;

}

auto readSquare(const Opts opts, size_t nInd)
{
  auto matrix = new double[](nInd * nInd);

  auto fileName = opts.kinship != "" ? opts.kinship : opts.Sinv;

  try
  {
    auto kinFile = File(fileName);
    kinFile.readln;
    size_t i = 0;

    foreach (line; kinFile.byLine)
    {
      auto splitLine = line.split.array.to!(double[]);
      enforce(splitLine.length == nInd && i < nInd,
          new InputException(" Incorrect number of columns."));
      foreach (j; iota(nInd))
        matrix[nInd * opts.kinshipLocations[i] + opts.kinshipLocations[j]] = splitLine[j];
      i++;
    }
    enforce(i == nInd, new InputException(" Incorrect number of rows."));
  }
  catch (Exception e)
  {
    stderr.writeln("Failed to read , ", fileName, " file.", e.msg);
    exit(1);
  }
  catch (RangeError e)
  {
    stderr.writeln("Failed to read ", fileName, ". File contains too few individuals.");
    exit(1);
  }

  return matrix;
}

auto readEigen(const Opts opts, size_t nInd)
{
  double[] eigen;
  try
  {
    eigen = File(opts.eigen).byLine.map!(a => to!double(a)).array;
    eigen = eigen.indexed(opts.kinshipLocations).array;
    enforce(eigen.length == nInd, new InputException(" Incorrect number of eigenvalues."));
  }
  catch (Exception e)
  {
    stderr.writeln("Failed to read eigenvalues.");
    exit(1);
  }
  catch (RangeError e)
  {
    stderr.writeln("Failed to read eigenvalues. Too few values in file.");
    exit(1);
  }

  return eigen;
}
