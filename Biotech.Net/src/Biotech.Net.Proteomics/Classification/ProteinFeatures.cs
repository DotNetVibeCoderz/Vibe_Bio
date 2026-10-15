using System.Text;

namespace Biotech.Net.Proteomics.Classification;

/// <summary>
/// Fixed-length numeric descriptors of protein sequences: amino-acid composition (20), dipeptide composition (400),
/// k-mer spectra and a handful of global physicochemical properties. Only the 20 standard residues are counted.
/// </summary>
public static class ProteinFeatures
{
    /// <summary>Amino-acid composition: 20 fractions in <see cref="AminoAcids.Standard"/> order (sum 1 over standard residues).</summary>
    public static double[] AminoAcidComposition(ReadOnlySpan<byte> residues)
    {
        var v = new double[20];
        int n = 0;
        foreach (byte b in residues)
        {
            int i = AminoAcids.IndexOf(b);
            if (i < 0) continue;
            v[i]++; n++;
        }
        if (n > 0) for (int i = 0; i < 20; i++) v[i] /= n;
        return v;
    }

    /// <summary>Dipeptide composition: 400 fractions, index <c>20·i + j</c> for the pair (residue i, residue j).</summary>
    public static double[] DipeptideComposition(ReadOnlySpan<byte> residues) => KmerSpectrum(residues, 2);

    /// <summary>
    /// Dense k-mer frequency vector of length 20^k (k = 1–3), index = base-20 number of the k-mer in
    /// <see cref="AminoAcids.Standard"/> order. Windows containing a non-standard residue are skipped.
    /// </summary>
    public static double[] KmerSpectrum(ReadOnlySpan<byte> residues, int k)
    {
        if (k is < 1 or > 3) throw new ArgumentOutOfRangeException(nameof(k), "Dense k-mer spectra support k = 1..3 (8000 features).");
        int size = 1;
        for (int i = 0; i < k; i++) size *= 20;
        var v = new double[size];
        int total = 0;
        foreach (var code in EnumerateKmers(residues, k)) { v[code]++; total++; }
        if (total > 0) for (int i = 0; i < size; i++) v[i] /= total;
        return v;
    }

    /// <summary>Sparse k-mer counts (k = 1–7), keyed by the base-20 code.</summary>
    public static Dictionary<int, int> KmerCounts(ReadOnlySpan<byte> residues, int k)
    {
        if (k is < 1 or > 7) throw new ArgumentOutOfRangeException(nameof(k), "Sparse k-mer codes support k = 1..7.");
        var d = new Dictionary<int, int>();
        foreach (var code in EnumerateKmers(residues, k))
            d[code] = d.TryGetValue(code, out int c) ? c + 1 : 1;
        return d;
    }

    /// <summary>The k-mer string for a base-20 code.</summary>
    public static string DecodeKmer(int code, int k)
    {
        var c = new char[k];
        for (int i = k - 1; i >= 0; i--) { c[i] = AminoAcids.Standard[code % 20]; code /= 20; }
        return new string(c);
    }

    /// <summary>
    /// Global descriptors: length (log10), GRAVY, net charge at pH 7 per residue, pI, aromaticity, aliphatic index / 100,
    /// instability index / 100 — seven values on comparable scales.
    /// </summary>
    public static double[] Physicochemical(ReadOnlySpan<byte> residues)
    {
        var a = new ProteinAnalysis(residues);
        return
        [
            Math.Log10(a.Length), a.Gravy(), a.NetCharge(7.0) / a.Length, a.IsoelectricPoint() / 14.0,
            a.Aromaticity(), a.AliphaticIndex() / 100.0, a.InstabilityIndex() / 100.0,
        ];
    }

    /// <summary>Composition (20) + dipeptide composition (400) + physicochemical (7): 427 features.</summary>
    public static double[] Combined(ReadOnlySpan<byte> residues)
    {
        var aac = AminoAcidComposition(residues);
        var dpc = DipeptideComposition(residues);
        var phys = Physicochemical(residues);
        var v = new double[aac.Length + dpc.Length + phys.Length];
        aac.CopyTo(v, 0);
        dpc.CopyTo(v, aac.Length);
        phys.CopyTo(v, aac.Length + dpc.Length);
        return v;
    }

    /// <summary>Applies a span-based extractor to a string.</summary>
    public static double[] Extract(Func<byte[], double[]> extractor, string residues) => extractor(Encoding.ASCII.GetBytes(residues));

    private static List<int> EnumerateKmers(ReadOnlySpan<byte> residues, int k)
    {
        var codes = new List<int>(Math.Max(0, residues.Length - k + 1));
        int mod = 1;
        for (int i = 0; i < k - 1; i++) mod *= 20;
        int code = 0, valid = 0;
        foreach (byte b in residues)
        {
            int idx = AminoAcids.IndexOf(b);
            if (idx < 0) { valid = 0; code = 0; continue; }
            code = (valid >= k - 1 ? code % mod : code) * 20 + idx;
            if (k == 1) code = idx;
            valid++;
            if (valid >= k) codes.Add(code);
        }
        return codes;
    }
}

/// <summary>Alignment-free sequence similarity from k-mer sets and spectra.</summary>
public static class SequenceSimilarity
{
    /// <summary>Jaccard index of the k-mer <em>sets</em>: |A ∩ B| / |A ∪ B| (1 for identical sets; 0 when both are empty).</summary>
    public static double KmerJaccard(ReadOnlySpan<byte> a, ReadOnlySpan<byte> b, int k = 3)
    {
        var ka = ProteinFeatures.KmerCounts(a, k);
        var kb = ProteinFeatures.KmerCounts(b, k);
        if (ka.Count == 0 && kb.Count == 0) return 0;
        int inter = 0;
        foreach (var key in ka.Keys) if (kb.ContainsKey(key)) inter++;
        return (double)inter / (ka.Count + kb.Count - inter);
    }

    /// <summary>Cosine similarity of k-mer count vectors.</summary>
    public static double KmerCosine(ReadOnlySpan<byte> a, ReadOnlySpan<byte> b, int k = 3)
    {
        var ka = ProteinFeatures.KmerCounts(a, k);
        var kb = ProteinFeatures.KmerCounts(b, k);
        double dot = 0, na = 0, nb = 0;
        foreach (var (key, c) in ka)
        {
            na += (double)c * c;
            if (kb.TryGetValue(key, out int d)) dot += (double)c * d;
        }
        foreach (var c in kb.Values) nb += (double)c * c;
        return na == 0 || nb == 0 ? 0 : dot / Math.Sqrt(na * nb);
    }

    /// <summary>Jaccard similarity of two strings.</summary>
    public static double KmerJaccard(string a, string b, int k = 3) => KmerJaccard(Encoding.ASCII.GetBytes(a), Encoding.ASCII.GetBytes(b), k);

    /// <summary>Cosine similarity of two strings.</summary>
    public static double KmerCosine(string a, string b, int k = 3) => KmerCosine(Encoding.ASCII.GetBytes(a), Encoding.ASCII.GetBytes(b), k);

    /// <summary>Symmetric all-against-all k-mer cosine similarity matrix (row-major n×n), computed in parallel.</summary>
    public static double[] CosineMatrix(IReadOnlyList<string> sequences, int k = 3)
    {
        int n = sequences.Count;
        var counts = new Dictionary<int, int>[n];
        var norms = new double[n];
        for (int i = 0; i < n; i++)
        {
            counts[i] = ProteinFeatures.KmerCounts(Encoding.ASCII.GetBytes(sequences[i]), k);
            double s = 0;
            foreach (var c in counts[i].Values) s += (double)c * c;
            norms[i] = Math.Sqrt(s);
        }
        var m = new double[n * n];
        Parallel.For(0, n, i =>
        {
            m[i * n + i] = norms[i] > 0 ? 1 : 0;
            for (int j = i + 1; j < n; j++)
            {
                var (small, large) = counts[i].Count <= counts[j].Count ? (counts[i], counts[j]) : (counts[j], counts[i]);
                double dot = 0;
                foreach (var (key, c) in small) if (large.TryGetValue(key, out int d)) dot += (double)c * d;
                double v = norms[i] == 0 || norms[j] == 0 ? 0 : dot / (norms[i] * norms[j]);
                m[i * n + j] = v;
                m[j * n + i] = v;
            }
        });
        return m;
    }
}
