module arg_parse;

import core.stdc.stdlib : exit;
import std.algorithm : canFind, countUntil, filter, joiner, map, setDifference,
  sort;
import std.array : array;
import std.array : split;
import std.conv : to;
import std.file : exists;
import std.getopt : arraySep, defaultGetoptFormatter, defaultGetoptPrinter,
  getopt;
import std.process : executeShell, pipeShell, Redirect, wait;
import std.range : indexed, iota;
import std.stdio : File, stderr, writeln;
import std.string : chomp;

class Opts
{

  //write appropriate string and quit
  bool version_ = false;
  bool verbose = false;
  //phenotype and genotype ids are given
  bool noheader = false;
  //number of genotype columns to skip, and phenotype column
  int genes = 0;
  int jobNumber = 0;
  //permutation numbers and seeds
  uint[] perms;
  //file names
  string vcf = "";
  string bed = "";
  string output = "";
  size_t window = 1_000_000;
  size_t[] genotypeLocations;
  size_t[] phenotypeLocations;
  size_t[] kinshipLocations;
  bool gt = false;
  bool calc = false;
  long loc = 0;

  bool nocheck = false;
  string kinship = "";
  string Sinv = "";
  string eigen = "";

  private auto parseOptions(string[] args)
  {
    arraySep = ",";
    // dfmt off
    auto options = getopt(args,
			  "bed", "Phenotype file [last argument].\n", &bed,
			  "vcf", "Genotype file.\n", &vcf,
			  "out|o", "Output file [stdout].\n", &output,
			  "kinship", "Specify kinship matrix for the sample.\n", &kinship,
			  "Sinv", "Specify the inverse kinship matrix for the sample.\n", &Sinv,
			  "eigen", "Specify the eigenvalues of the kinship matrix.\n", &eigen,
			  "calc-h2", "Calculate the heritability using the regress package in R.\n", &calc,
			  "job-number", "Split the analysis into a number of smaller runs which can run in parallel on a cluster. This option specifies which of the sub-analyses should be run.\n", &jobNumber,
			  "genes", "This specifies the number of genes to be analysed in each job.\n", &genes,
			  "perm", "Number of permutations to use. One following number indicates the number of permutations, two comma separated numbers gives the number of permutations and the seed.\n", &perms,
			  "window", "The size in base pairs of the cis window around the transcription start site of the gene [1,000,000].\n", &window,
			  "noheader", "Suppress writing of header line.\n", &noheader,
			  "nocheck", "Do not check IDs, assume correct order.\n", &nocheck,
			  "verbose", "Print additional information.\n", &verbose,
			  "version", "Display version information.\n", &version_,
			  );
    return options;
  }

  this(string[] args)
  {
    immutable bool noArgs = args.length == 1;
    auto options = parseOptions(args);
    try
    {
      // dfmt on
    if (options.helpWanted || noArgs)
    {
      defaultGetoptPrinter("Alfonso's method for eQTL mapping with related individuals.
",
          options.options);
      exit(0);
    }

    if (version_)
      giveHelp(versionString);

    immutable auto checkTabix = executeShell("command -v tabix");

    if (checkTabix.status != 0)
    {
      stderr.writeln("Error: tabix is not installed.");
      exit(1);
    }

    if (!vcf.exists)
    {
      stderr.writeln("Error: genotype file ", vcf, " does not exist.");
      exit(1);
    }

    if (!(vcf ~ ".tbi").exists && !(vcf ~ ".csi").exists)
    {
      stderr.writeln("Error: Neither ", vcf, ".tbi nor ", vcf,
          ".csi files are present, meaning genotype file hasn't been indexed with tabix or bcftools.");
      exit(1);
    }

    if (verbose)
    {
      stderr.writeln("Alfonso version, git commit: ", commitString, ".");
    }

    matchIds();

    if (bed == "" && args.length > 1)
      bed = args[$ - 1];

    if (perms.length == 0)
      perms = [10_000];

  }
  catch (Exception e)
  {
    stderr.writeln("Error with command: ", e.msg, "\n");

    auto placeholder = ["dummy"];
    auto helpOptions = parseOptions(placeholder);
    defaultGetoptFormatter(stderr.lockingTextWriter(),
        "Alfonso's method for eQTL mapping with related individuals.
", helpOptions.options);
    exit(1);
  }

}

private void matchIds()
{
  string[] phenotypeIds;

  try
  {
    auto bedFile = File(bed);

    phenotypeIds = bedFile.readln.chomp.split[4 .. $];

    if (verbose)
    {
      stderr.writeln(phenotypeIds.length, " individuals present in phenotype file.");
    }
  }
  catch (Exception e)
  {
    stderr.writeln("Failed to read phenotype IDs. ", e.msg);
    exit(1);
  }

  if (nocheck)
  {
    genotypeLocations = iota(phenotypeIds.length).array;
    phenotypeLocations = iota(phenotypeIds.length).array;
    kinshipLocations = iota(phenotypeIds.length).array;

    if (verbose)
    {
      stderr.writeln("Assuming same individuals in genotype and kinship files");
    }
  }
  else
  {
    string[] genotypeIds;
    try
    {
      auto pipes = pipeShell("zcat " ~ vcf ~ " | grep -v '##' | head -2", Redirect.stdout);
      scope (exit)
        wait(pipes.pid);

      auto line = pipes.stdout.readln.chomp;

      genotypeIds = line.split[9 .. $].to!(string[]);

      if (verbose)
      {
        stderr.writeln(genotypeIds.length, " individuals present in genotype file.");
      }

      auto formatField = pipes.stdout.readln.chomp.split[8].split(':');
      loc = countUntil(formatField, "DS");
      if (loc == -1)
      {
        loc = countUntil(formatField, "GT");
        gt = true;
        if (loc == -1)
        {
          stderr.writeln("DS and GT fields are both missing from the vcf file.");
          exit(1);
        }
      }
    }
    catch (Exception e)
    {
      stderr.writeln("Failed to read genotype IDs. ", e.msg);
      exit(1);
    }

    phenotypeLocations = genotypeIds.map!(a => phenotypeIds.countUntil(a))
      .filter!(a => a != -1).array.to!(size_t[]);

    genotypeLocations = iota(genotypeIds.length).filter!(
        a => phenotypeIds.canFind(genotypeIds[a])).array;

    if (genotypeLocations.length == 0 || phenotypeLocations.length == 0)
    {
      stderr.writeln("No individuals to analyse.");
      exit(1);
    }

    if (phenotypeIds.indexed(phenotypeLocations)
        .array != genotypeIds.indexed(genotypeLocations).array)
    {
      stderr.writeln("Failed to match genotype and phenotype IDs. THIS SHOULD NEVER HAPPEN.");
      exit(1);
    }

    if (verbose && genotypeLocations.length != genotypeIds.length)
    {
      stderr.writeln(genotypeIds.indexed(setDifference(iota(genotypeIds.length),
          genotypeLocations)).joiner(", "), " dropped from genotype file.");
    }

    if (verbose && phenotypeLocations.length != phenotypeIds.length)
    {
      stderr.writeln(phenotypeIds.indexed(setDifference(iota(phenotypeIds.length),
          phenotypeLocations.dup.sort!())).joiner(", "), " dropped from phenotype file.");
    }

    string[] kinshipIds;

    try
    {
      if (kinship != "")
        kinshipIds = File(kinship).readln.chomp.split;
      else
        kinshipIds = File(Sinv).readln.chomp.split;

      if (verbose)
      {
        stderr.writeln(kinshipIds.length, " individuals present in genotype file.");
      }

      auto temp = genotypeIds.indexed(genotypeLocations).map!(a => kinshipIds.countUntil(a)).array;
      if (temp.canFind(-1))
      {
        auto missing = genotypeIds.indexed(genotypeLocations)
          .filter!(a => !kinshipIds.canFind(a)).joiner(", ");
        auto fileName = kinship != "" ? "kinship" : "Sinv";
        stderr.writeln("Individuals in genotype and phenotype file, but not in ",
            fileName, " matrix. Missing individuals are ", missing, ".");
        exit(1);
      }

      kinshipLocations = temp.to!(size_t[]);

      if (verbose && kinshipLocations.length != kinshipIds.length)
      {
        stderr.writeln(kinshipIds.indexed(setDifference(iota(kinshipIds.length),
            kinshipLocations.dup.sort!())).joiner(", "), " dropped from kinship file.");
      }
    }
    catch (Exception e)
    {
      stderr.writeln("Failed to read kinship IDs. ", e.msg);
      exit(1);
    }

    if (kinshipIds.indexed(kinshipLocations).array != genotypeIds.indexed(genotypeLocations).array)
    {
      stderr.writeln("Failed to match genotype and Sinv IDs. THIS SHOULD NEVER HAPPEN.");
      exit(1);
    }
  }
}
}

static immutable string versionString = "Alfonso's method for eQTL mapping with related individuals, version 1.0-";
static immutable string commitString = chomp(cast(string) import("commit"));

void giveHelp(string versionString)
{
  import std.compiler : name, version_major, version_minor;

  static string[] dateString = __DATE__.split;

  writeln(versionString, commitString, ".");
  writeln("Compiled with ", name, " ", version_major, ".", version_minor,
      " at ", __TIME__, ", ", dateString[1], " ", dateString[0], " ", dateString[2], ".");

  exit(0);
}
