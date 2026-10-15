using System.Text;
using Biotech.Net.Core;

namespace Biotech.Net.Genomics.Tests;

/// <summary>Deterministic genome and read simulators for the tests.</summary>
internal static class Simulation
{
    public static string RandomDna(Random rnd, int length, double gc = 0.5)
    {
        var c = new char[length];
        for (int i = 0; i < length; i++)
        {
            bool strong = rnd.NextDouble() < gc;
            c[i] = strong ? (rnd.Next(2) == 0 ? 'G' : 'C') : (rnd.Next(2) == 0 ? 'A' : 'T');
        }
        return new string(c);
    }

    public static string ReverseComplement(string s) => Sequence.Dna(s).ReverseComplement().ToString();

    public sealed record SimRead(SequenceRead Read, int Origin, bool Reverse, int Errors);

    /// <summary>
    /// Reads of <paramref name="length"/> bases from <paramref name="genome"/> (uniform starts, random strand) with
    /// substitution rate <paramref name="errorRate"/> and quality <paramref name="quality"/>.
    /// </summary>
    public static List<SimRead> Reads(Random rnd, string genome, int count, int length, double errorRate, byte quality = 30, string prefix = "r")
    {
        var reads = new List<SimRead>(count);
        for (int i = 0; i < count; i++)
        {
            int start = rnd.Next(genome.Length - length + 1);
            var s = genome.Substring(start, length).ToCharArray();
            int errors = 0;
            for (int j = 0; j < s.Length; j++)
            {
                if (rnd.NextDouble() >= errorRate) continue;
                char b;
                do b = "ACGT"[rnd.Next(4)]; while (b == s[j]);
                s[j] = b;
                errors++;
            }
            bool rev = rnd.Next(2) == 1;
            string bases = rev ? ReverseComplement(new string(s)) : new string(s);
            var q = new byte[length];
            Array.Fill(q, quality);
            reads.Add(new SimRead(SequenceRead.FromBuffer(Sequence.Dna(bases, $"{prefix}{i}"), q), start, rev, errors));
        }
        return reads;
    }

    /// <summary>Applies substitutions (0-based position → base).</summary>
    public static string Mutate(string genome, params (int Pos, char Base)[] changes)
    {
        var sb = new StringBuilder(genome);
        foreach (var (p, b) in changes) sb[p] = b;
        return sb.ToString();
    }
}
