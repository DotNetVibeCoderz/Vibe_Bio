using System.Diagnostics;
using System.Text;
using Biotech.Net.Core;

namespace Biotech.Net.WebGallery.Services;

/// <summary>Bundled demo data and small helpers shared by pages (mirrors the desktop Gallery's GalleryData).</summary>
public static class GalleryData
{
    /// <summary>EGFP coding sequence (720 bp), the default DNA Explorer and CRISPR target.</summary>
    public const string Egfp =
        "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA";

    public static string EgfpFasta => ">EGFP enhanced green fluorescent protein CDS\n" + Wrap(Egfp, 60);

    /// <summary>Text of a resource embedded from the desktop Gallery's Data folder (region.fa, region.gff3, region.vcf).</summary>
    public static string Resource(string name)
    {
        using var s = typeof(GalleryData).Assembly.GetManifestResourceStream("WebGallery.Data." + name)
            ?? throw new InvalidOperationException($"Missing embedded resource {name}");
        using var r = new StreamReader(s);
        return r.ReadToEnd();
    }

    public static string Wrap(string s, int width)
    {
        var sb = new StringBuilder();
        for (int i = 0; i < s.Length; i += width) sb.Append(s.AsSpan(i, Math.Min(width, s.Length - i))).Append('\n');
        return sb.ToString();
    }

    /// <summary>Back-translates a protein with the first codon of each residue (drives a protein page's trace).</summary>
    public static string BackTranslate(string protein)
    {
        var sb = new StringBuilder(protein.Length * 3);
        foreach (char aa in protein)
        {
            var codons = GeneticCode.Standard.CodonsFor(char.ToUpperInvariant(aa));
            sb.Append(codons.Count > 0 ? codons[0] : "NNN");
        }
        return sb.ToString();
    }

    /// <summary>Maps arbitrary text (SMILES, code) onto bases so every page has a trace.</summary>
    public static string TextToBases(string text)
    {
        var sb = new StringBuilder(text.Length * 2);
        foreach (char c in text) { sb.Append("ACGT"[c & 3]); sb.Append("ACGT"[(c >> 2) & 3]); }
        return sb.ToString();
    }

    public static string Ms(Stopwatch sw) => $"{sw.Elapsed.TotalMilliseconds:0.0} ms";

    public static string Preview(string s, int n) => s.Length <= n ? s : s[..n] + "…";

    /// <summary>C# string literal of a (possibly long) sequence.</summary>
    public static string Literal(string s) => "\"" + s.Replace("\\", "\\\\").Replace("\"", "\\\"") + "\"";
}
