module run_analysis;

import calculation : betaParameters, corPvalue, correlation, Opts, transform,
  VarianceException;
import read_data : Genotype, Phenotype, readGenotype;
import std.algorithm : count, map, max, sort;
import std.array : array;
import std.conv : to;
import std.format : format;
import std.math : fabs, sqrt;
import std.numeric : dotProduct;
import std.range : chunks, indexed, iota, SearchPolicy, zip;
import std.stdio : File, stderr, writeln;

enum double EPSILON = 0.00000001;

class InputException : Exception
{
  pure this(string s)
  {
    super(s);
  }
}

pure nothrow extern (C)
{
  //call GSL to calculate P values from beta parameters
  double gsl_cdf_beta_P(double x, double alpha, double beta);
}

void analyseData(Phenotype phenotype, size_t[] perms, File outFile, const Opts opts,
    double[] Sinv, double[] eigen)
{
  auto genotype = readGenotype(opts, phenotype.chromosome, phenotype.location, phenotype.geneName);

  if (opts.verbose)
  {
    stderr.writeln("Extracted ", genotype.length, " usable genotypes.");
  }

  mixedModelsRun(opts, perms, phenotype, genotype, outFile, Sinv, eigen);
}

void mixedModelsRun(const Opts opts, size_t[] perms, Phenotype phenotype,
    Genotype[] genotypes, File outFile, double[] Sinv, double[] eigen)
{
  immutable size_t nInd = phenotype.values.length;

  double[] multSinv(const double[] vec)
  {
    return Sinv.chunks(nInd).map!(a => dotProduct(a, vec)).array;
  }

  immutable size_t nPerm = opts.perms[0];
  double[2] cor;

  try
  {
    transform(phenotype.values);
  }
  catch (VarianceException)
  {
    stderr.writeln("Phenotype is constant");
    return;
  }
  const auto yInv = multSinv(phenotype.values);

  auto maxCor = new double[](nPerm);
  //we need to store greatest statistic across all SNPs in maxCor
  maxCor[] = 0.0;
  foreach (ref genotype; genotypes)
  {
    try
    {
      transform(genotype.values);

      auto effect = dotProduct(genotype.values, phenotype.values);
      effect = effect * effect;
      auto myh2 = (phenotype.h2 - effect) / (1 - effect);

      myh2 = myh2 > 0 ? myh2 : 0;
      const auto Dp = eigen.map!(a => 1 / sqrt((myh2 * a + 1 - myh2) * (1 - effect))).array;

      auto xInv = multSinv(genotype.values);
      auto yTemp = yInv.dup;

      foreach (i; iota(nInd))
      {
        xInv[i] *= Dp[i];
        yTemp[i] *= Dp[i];
      }

      transform(xInv);
      transform(yTemp);
      cor = correlation(yTemp, xInv);

      genotype.cor = fabs(cor[0]) - EPSILON;

      // perm P value as before,and store maximum correlation for each permutation
      auto simplePerm = map!(a => fabs(dotProduct(yTemp, xInv.indexed(a))))(chunks(perms, nInd))
        .array;

      genotype.snpId ~= format("%g\t%g\t%g", cor[0], cor[1],
          (1.0 + simplePerm.count!(a => a > genotype.cor)) / (nPerm + 1));
      foreach (e; zip(iota(nPerm), simplePerm))
        maxCor[e[0]] = max(maxCor[e[0]], e[1]);

    }
    catch (VarianceException e)
    {
      genotype.cor = 2;
    }
  }

  //sort stored maximum statistics
  auto sortMax = sort!()(maxCor);
  immutable auto len = (sortMax.length + 1).to!double;

  auto minPvalues = sortMax.map!(a => corPvalue(a, nInd)).array;
  auto betaParam = betaParameters(minPvalues);

  if (opts.verbose)
  {
    outFile.writefln("##Beta parameters for %s: %.8g and %.8g.", phenotype.geneName, betaParam[0], betaParam[1]);
    outFile.writefln("##Heritability for %s: %g.", phenotype.geneName, phenotype.h2);
  }

  // //read through old file and compare correlations to maxCor to calculate FWER

  foreach (genotype; genotypes)
  {
    if (genotype.cor == 2)
    {
      writeln(phenotype.geneName, genotype.snpId, "\tNA\tNA\tNA\tNA\tNA");
    }
    else
    {
      auto adjusted = (sortMax.upperBound!(SearchPolicy.gallop)(genotype.cor).length + 1) / len;
      auto betaP = gsl_cdf_beta_P(corPvalue(genotype.cor, nInd), betaParam[0], betaParam[1]);
      outFile.writeln(phenotype.geneName, genotype.snpId, "\t", adjusted, "\t", betaP);
    }
  }

}
