using BenchmarkDotNet.Attributes;
using Biotech.Net.Chemistry;

namespace Biotech.Net.Benchmarks;

/// <summary>SMILES parsing, Morgan (ECFP4) fingerprints and Tanimoto similarity search.</summary>
public class ChemistryBenchmarks
{
    public static readonly string[] DrugSmiles =
    [
        "CC(=O)Oc1ccccc1C(=O)O",                                  // aspirin
        "CC(=O)Nc1ccc(O)cc1",                                     // paracetamol
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",                             // ibuprofen
        "Cn1cnc2c1c(=O)n(C)c(=O)n2C",                             // caffeine
        "CN1CCC[C@H]1c1cccnc1",                                   // nicotine
        "COc1ccc2[nH]cc(CCN(C)C)c2c1",                            // 5-MeO-DMT
        "CC(C)NCC(O)COc1cccc2ccccc12",                            // propranolol
        "CN(C)C(=N)NC(=N)N",                                      // metformin
        "Clc1ccc(cc1)C(c1ccccc1)N1CCN(CC1)CCOCC(=O)O",            // cetirizine
        "CC1=C(C(=O)OC)C(c2ccccc2[N+](=O)[O-])C(C(=O)OC)=C1C",    // nifedipine
        "O=C(O)c1cc(N=Nc2ccc(O)c(C(=O)O)c2)ccc1O",                // olsalazine
        "CC(C)(C)NCC(O)c1ccc(O)c(CO)c1",                          // salbutamol
        "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc12",                    // diazepam
        "COc1ccc2nc(sc2c1)S(=O)Cc1ncc(C)c(OC)c1C",                // omeprazole-like
        "OC(=O)CCCc1ccc(N(CCCl)CCCl)cc1",                         // chlorambucil
        "Nc1ccc(cc1)S(=O)(=O)Nc1ccccn1",                          // sulfapyridine
        "CCN(CC)CCNC(=O)c1ccc(N)cc1",                             // procainamide
        "C1CCC(CC1)NC(=O)NS(=O)(=O)c1ccc(cc1)C(=O)C",             // acetohexamide-like
        "OCC1OC(O)C(O)C(O)C1O",                                   // glucose
        "CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2O",                    // testosterone
    ];

    public const int LibrarySize = 10_000;

    private Molecule[] _molecules = [];
    private FingerprintIndex _index = new();
    private BitVector _query = new(2048);

    [GlobalSetup]
    public void Setup()
    {
        _molecules = DrugSmiles.Select((s, i) => Molecule.FromSmiles(s, $"m{i}")).ToArray();
        // A 10k-entry library: the 20 drug fingerprints with random extra bits so entries differ.
        var rng = new Random(11);
        var baseFps = _molecules.Select(m => MorganFingerprint.Ecfp4(m)).ToArray();
        _index = new FingerprintIndex(m => MorganFingerprint.Ecfp4(m));
        for (int i = 0; i < LibrarySize; i++)
        {
            var fp = baseFps[i % baseFps.Length].Clone();
            for (int b = 0; b < 8; b++) fp.Set(rng.Next(2048));
            _index.Add($"lib{i}", fp);
        }
        _query = baseFps[0];
    }

    [Benchmark(Description = "Parse 20 drug SMILES")]
    public int ParseSmiles()
    {
        int atoms = 0;
        foreach (var s in DrugSmiles) atoms += Molecule.FromSmiles(s).AtomCount;
        return atoms;
    }

    [Benchmark(Description = "Morgan r=2 2048-bit x 20 molecules")]
    public int Morgan()
    {
        int bits = 0;
        foreach (var m in _molecules) bits += MorganFingerprint.Compute(m, 2, 2048).Count;
        return bits;
    }

    [Benchmark(Description = "Tanimoto top-10 over 10,000 fingerprints")]
    public int TanimotoTopK() => _index.TopK(_query, 10).Count;

    [Benchmark(Description = "Bulk Tanimoto, 10,000 fingerprints")]
    public double BulkTanimoto()
    {
        var fps = new BitVector[_index.Count];
        for (int i = 0; i < fps.Length; i++) fps[i] = _index[i];
        return Similarity.BulkTanimoto(_query, fps)[0];
    }
}
