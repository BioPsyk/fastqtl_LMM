double h2Calc(ref double[] phenotype, ref double[] kinship, string geneName, bool verbose)
{
  import std.array : join, split;
  import std.conv : to;
  import std.format : format;
  import std.process : pipeProcess, Redirect, wait;
  import std.stdio : File, readln, stderr, writeln;
  import std.string : strip;

  //   auto rString = "library(regress)
  // phenotype <- c(%-(%a, %))
  // kinship <- matrix(c(%-(%a, %)), %d, %d, byrow = TRUE)
  // model <- regress(phenotype ~ 1, ~ kinship)$sigma
  // write(sprintf('%s', model[1] / (model[1] + model[2])), '')
  // ".format(phenotype, kinship, phenotype.length,
  //       kinship.length / phenotype.length, "%a").split("\n").join("; ");

  auto rString = "library(GenABEL, quietly = TRUE)
  suppressMessages(library(hglm, quietly = TRUE))
  phenotype <- c(%-(%a, %))
  kinship <- matrix(c(%-(%a, %)), %d, %d, byrow = TRUE)
  data <- data.frame(ids = 1:length(phenotype), phenotype = phenotype, kinship = kinship)
  h2 <- polygenic_hglm(phenotype ~ 1, kinship, data = data, quiet = TRUE)$esth2
  write(sprintf('%s', h2), '')
  ".format(phenotype, kinship, phenotype.length,
      kinship.length / phenotype.length, "%a");

  import std.uuid : randomUUID;
  import std.file : exists, remove;

  auto rFile = randomUUID.toString;

  while (rFile.exists)
  {
    rFile = randomUUID.toString;
  }

  scope (exit)
  {
    if (rFile.exists)
      rFile.remove;
  }

  auto scriptFile = File(rFile, "w");

  scriptFile.writeln(rString);

  scriptFile.flush;

  try
  {
    auto rProcess = pipeProcess(["Rscript", "--vanilla", rFile], Redirect.stdout);

    scope (exit)
      wait(rProcess.pid);

    double h2 = rProcess.stdout.readln.strip.to!double;

    if (h2 < 0)
    {
      stderr.writeln("Heritability estimate for ", geneName, " is ", h2,
          " (less than zero). Setting equal to 0.");
      return 0;
    }
    if (h2 > 1)
    {
      stderr.writeln("Heritability estimate for ", geneName, " is ", h2,
          " (greater than one). Setting equal to 1.");
      return 1;
    }
    if (verbose)
    {
      stderr.writeln("Heritability estimate for ", geneName, " is ", h2, ".");
    }
    return h2;
  }
  catch (Exception e)
  {
    stderr.writeln("Warning, failed to estimate heritability for ", geneName,
        ". Setting h2 to 0.");
    return 0;
  }

}

version (EMBEDR)
{
  double h2CalcEmbedR(ref double[] phenotype, ref double[] kinship, string geneName, bool verbose)
  {
    import embedr.r : evalRQ, RMatrix, RVector, scalar, toR;
    import std.range : chunks;
    import std.stdio : stderr, writeln;

    auto nInd = cast(int) phenotype.length;
    auto v = RVector(phenotype);

    v.toR("rV");

    auto m = RMatrix(nInd, nInd);

    int i = 0;
    foreach (e; kinship.chunks(nInd))
    {
      int j = 0;
      foreach (f; e)
      {
        m[i, j] = f;
        j++;
      }
      i++;
    }

    m.toR("rMat");

    try
    {
      evalRQ(`library(regress)`);
      evalRQ(`model <- regress(rV ~ 1, ~ rMat)$sigma`);
      evalRQ(`sig <- model[1] / (model[1] + model[2])`);

      auto h2 = scalar("sig");

      if (h2 < 0)
      {
        stderr.writeln("Heritability estimate for ", geneName, " is ", h2,
            " (less than zero). Setting equal to 0.");
        return 0;
      }
      if (h2 > 1)
      {
        stderr.writeln("Heritability estimate for ", geneName, " is ", h2,
            " (greater than one). Setting equal to 1.");
        return 1;
      }
      if (verbose)
      {
        stderr.writeln("Heritability estimate for ", geneName, " is ", h2, ".");
      }

      return h2;

    }
    catch (Exception e)
    {
      stderr.writeln("Warning, failed to estimate heritability for ", geneName,
          ". Setting h2 to 0.");
      return 0;
    }

  }

}
