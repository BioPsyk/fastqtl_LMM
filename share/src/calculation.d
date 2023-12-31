module calculation;

import arg_parse : Opts;
import std.algorithm : sum;
import std.conv : to;
import std.exception : enforce;
import std.math : fabs, sqrt;

version (unittest)
{
  import std.math : approxEqual;
}

class VarianceException : Exception
{
  //thrown if variable is constant
  pure nothrow this(string s)
  {
    super(s);
  }
}

pure nothrow extern (C)
{
  //call GSL to calculate P values from T statistics
  double gsl_cdf_tdist_P(double x, double nu);
}

extern (C)
{
  int mleBeta(double* pval, size_t nPerm, double* alpha, double* beta);
}

@system unittest
{
  // Checks GSL gives right p value for t statistic
  assert(approxEqual(gsl_cdf_tdist_P(-1.6, 7), 0.07681585));
}

pure void transform(ref double[] vector)
{
  //transforms array so mean =0 sum of squares = 1
  int n = 0;
  double mean = 0;
  double M2 = 0;
  double delta;

  foreach (ref e; vector)
  {
    n++;
    delta = e - mean;
    mean += delta / n;
    M2 += delta * (e - mean);
  }

  enforce(M2 != 0, new VarianceException(""));

  M2 = sqrt(M2);

  foreach (ref e; vector)
    e = (e - mean) / M2;
}

@system unittest
{
  //Checks that transform works on randomly generated vector
  import std.algorithm : reduce;
  import std.random : uniform;

  double[] x = new double[10];
  foreach (ref e; x)
    e = uniform(0.0, 10.0);

  transform(x);
  auto mean = 0.0.reduce!((a, b) => a + b)(x);

  assert(approxEqual(mean, 0.0));
  assert(approxEqual(0.0.reduce!((a, b) => a + (b - mean) * (b - mean))(x), 1));
}

pure nothrow double[2] correlation(ref double[] vector1, ref double[] vector2)
{
  //calculates correlation, t stat and p value for two arrays
  import std.numeric : dotProduct;

  double[2] results;
  results[0] = dotProduct(vector1, vector2);
  results[1] = gsl_cdf_tdist_P(
      -fabs(results[0] * sqrt((vector1.length - 2) / (1 - results[0] * results[0]))),
      vector1.length - 2) * 2;
  return results;
}

@system unittest
{
  //Check correlation of phenotype with 3rd row genotype against estimates from R
  double[2] corFromR = [-0.3060383, 0.3897973];

  double[] genotype = [0.115, 2, 0.0964, 1, 1, 1, 0, 1, 0, 0.0563];
  double[] phen = [
    -1.3853088072, -0.785797093643, 1.14540423638, -0.785797093643, 1.03820492508,
    -1.25652676836, -0.787662180447, -2.05355237841, -0.245457234103, 1.14277217712
  ];

  transform(phen);
  transform(genotype);
  double[2] cor = correlation(genotype, phen);

  assert(approxEqual(cor[0], corFromR[0]));
  assert(approxEqual(cor[1], corFromR[1]));
}

size_t[] genPerms(const Opts opts, size_t nInd)
{
  import std.array : array;
  import std.random : randomShuffle, rndGen;
  import std.range : chunks, cycle, iota, take;

  if (opts.perms.length > 1)
    rndGen.seed(opts.perms[1]);

  size_t[] outPerm = iota(nInd).cycle.take(opts.perms[0] * nInd).array;

  foreach (ref perm; chunks(outPerm, nInd))
    randomShuffle(perm);

  return outPerm;

}

double corPvalue(double cor, size_t nInd)
{
  return gsl_cdf_tdist_P(-fabs(cor * sqrt((nInd - 2) / (1 - cor * cor))), nInd - 2) * 2;
}

double[2] betaParameters(ref double[] minPval)
{
  immutable double mean = minPval.sum / minPval.length;
  double variance = 0;
  foreach (ref e; minPval)
    variance += (e - mean) * (e - mean);
  variance /= minPval.length;
  //Estimate shape1 & shape2
  double alpha = mean * (mean * (1 - mean) / variance - 1);
  double beta = alpha * (1 / mean - 1);

  immutable auto alphaCopy = alpha;
  immutable auto betaCopy = beta;

  immutable auto success = mleBeta(minPval.ptr, minPval.length, &alpha, &beta);

  if (success == 0)
    return [alphaCopy, betaCopy];

  return [alpha, beta];
}
