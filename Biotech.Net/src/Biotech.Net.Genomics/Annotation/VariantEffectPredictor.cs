using System.Globalization;
using System.Text;
using Biotech.Net.Core;
using Biotech.Net.Genomics.Vcf;

namespace Biotech.Net.Genomics.Annotation;

/// <summary>Options for <see cref="VariantEffectPredictor"/>.</summary>
public sealed record VariantEffectOptions
{
    /// <summary>Window 5' of a transcript reported as <c>upstream_gene_variant</c> (default 5,000 bp, as VEP).</summary>
    public int UpstreamDistance { get; init; } = 5000;

    /// <summary>Window 3' of a transcript reported as <c>downstream_gene_variant</c> (default 5,000 bp).</summary>
    public int DownstreamDistance { get; init; } = 5000;

    /// <summary>Genetic code used for translation (standard by default).</summary>
    public GeneticCode Code { get; init; } = GeneticCode.Standard;

    /// <summary>Exonic bases next to an intron that form the splice region (default 3).</summary>
    public int SpliceRegionExonBases { get; init; } = 3;

    /// <summary>Intronic bases from each intron end that form the splice region (default 8; bases 1–2 are the donor/acceptor sites).</summary>
    public int SpliceRegionIntronBases { get; init; } = 8;
}

/// <summary>
/// Predicts the consequences of variants on annotated transcripts, using Sequence Ontology terms (VEP-like) and
/// producing HGVS <c>c.</c>/<c>p.</c> notation for the common cases (substitutions, deletions, duplications,
/// insertions, delins, frameshifts, nonsense, stop-loss extensions and start loss).
/// </summary>
/// <remarks>
/// Gene models come from the <see cref="Genome"/> (genes → transcripts with exons and an optional
/// <see cref="Transcript.CodingRegion"/>); coding consequences need the chromosome sequence to be loaded.
/// Alleles are trimmed to their minimal representation, and indels in transcripts are shifted 3' (HGVS rule)
/// before the <c>c.</c> notation is written.
/// </remarks>
public sealed class VariantEffectPredictor
{
    private readonly Genome _genome;
    private readonly Dictionary<string, (Gene[] Genes, long MaxLength)> _index = new(StringComparer.Ordinal);
    private readonly Dictionary<Transcript, TxModel> _models = new(ReferenceEqualityComparer.Instance);

    /// <summary>Creates a predictor over the genes of <paramref name="genome"/>.</summary>
    public VariantEffectPredictor(Genome genome, VariantEffectOptions? options = null)
    {
        _genome = genome ?? throw new ArgumentNullException(nameof(genome));
        Options = options ?? new VariantEffectOptions();
        foreach (var g in genome.Genes.GroupBy(g => g.Location.Chromosome))
        {
            var arr = g.OrderBy(x => x.Location.Start).ToArray();
            _index[g.Key] = (arr, arr.Length == 0 ? 0 : arr.Max(x => x.Location.Length));
        }
    }

    /// <summary>The options in use.</summary>
    public VariantEffectOptions Options { get; }

    /// <summary>Effects of every alternate allele of <paramref name="variant"/> on every nearby transcript (intergenic when none).</summary>
    public IReadOnlyList<VariantEffect> Predict(Variant variant)
    {
        var result = new List<VariantEffect>();
        for (int a = 0; a < variant.Alternates.Count; a++)
            result.AddRange(Predict(variant.Chromosome, variant.Position, variant.Reference, variant.Alternates[a]));
        return result;
    }

    /// <summary>Effects of one allele given in VCF style (1-based <paramref name="position"/>, anchored REF/ALT).</summary>
    public IReadOnlyList<VariantEffect> Predict(string chromosome, long position, string reference, string alternate)
    {
        var result = new List<VariantEffect>();
        if (alternate is "*" or "." || alternate.StartsWith('<') || alternate.Contains('[') || alternate.Contains(']'))
            return result;
        string refU = reference.ToUpperInvariant(), altU = alternate.ToUpperInvariant();
        long s = position - 1;
        // Minimal representation: trim the common suffix, then the common prefix.
        int suffix = 0;
        while (suffix < refU.Length && suffix < altU.Length && refU[^(suffix + 1)] == altU[^(suffix + 1)]) suffix++;
        refU = refU[..^suffix]; altU = altU[..^suffix];
        int prefix = 0;
        while (prefix < refU.Length && prefix < altU.Length && refU[prefix] == altU[prefix]) prefix++;
        refU = refU[prefix..]; altU = altU[prefix..];
        s += prefix;
        long e = s + refU.Length;
        if (refU.Length == 0 && altU.Length == 0) return result; // no change

        if (_index.TryGetValue(chromosome, out var idx))
        {
            long window = Math.Max(Options.UpstreamDistance, Options.DownstreamDistance);
            long lo = s - window - idx.MaxLength, hi = e + window;
            int i = LowerBound(idx.Genes, lo);
            for (; i < idx.Genes.Length && idx.Genes[i].Location.Start <= hi; i++)
            {
                var gene = idx.Genes[i];
                if (gene.Location.End < s - window) continue;
                if (gene.Transcripts.Count == 0)
                {
                    var tx = new Transcript(gene.Id, gene.Id, [gene.Location]);
                    if (Evaluate(chromosome, s, e, refU, altU, alternate, gene, tx) is { } eff0) result.Add(eff0);
                    continue;
                }
                foreach (var tx in gene.Transcripts)
                    if (Evaluate(chromosome, s, e, refU, altU, alternate, gene, tx) is { } eff) result.Add(eff);
            }
        }
        if (result.Count == 0)
            result.Add(new VariantEffect { Allele = alternate, Consequences = [Consequence.IntergenicVariant] });
        result.Sort((x, y) => x.MostSevere.CompareTo(y.MostSevere));
        return result;
    }

    /// <summary>Effect of one allele on one transcript; null when the variant is farther than the up/downstream window.</summary>
    public VariantEffect? PredictTranscript(Variant variant, int alternateIndex, Gene gene, Transcript transcript)
    {
        string refU = variant.Reference.ToUpperInvariant(), altU = variant.Alternates[alternateIndex].ToUpperInvariant();
        long s = variant.Position - 1;
        int suffix = 0;
        while (suffix < refU.Length && suffix < altU.Length && refU[^(suffix + 1)] == altU[^(suffix + 1)]) suffix++;
        refU = refU[..^suffix]; altU = altU[..^suffix];
        int prefix = 0;
        while (prefix < refU.Length && prefix < altU.Length && refU[prefix] == altU[prefix]) prefix++;
        refU = refU[prefix..]; altU = altU[prefix..];
        s += prefix;
        return Evaluate(variant.Chromosome, s, s + refU.Length, refU, altU, variant.Alternates[alternateIndex], gene, transcript);
    }

    /// <summary>
    /// Adds a VEP-style <c>CSQ</c> INFO field (see <see cref="VariantEffect.ToCsq"/>) to the record and returns the effects.
    /// Call <see cref="AddCsqHeader"/> once on the output header.
    /// </summary>
    public IReadOnlyList<VariantEffect> Annotate(VcfRecord record)
    {
        var effects = Predict(record.Variant);
        if (effects.Count > 0) record.Info["CSQ"] = string.Join(',', effects.Select(e => e.ToCsq()));
        return effects;
    }

    /// <summary>Declares the <c>CSQ</c> INFO field.</summary>
    public static void AddCsqHeader(VcfHeader header) => header.AddInfo("CSQ", ".", "String",
        "Consequence annotations. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature|HGVSc|HGVSp|Protein_position|Amino_acids|Codons");

    private static int LowerBound(Gene[] genes, long start)
    {
        int lo = 0, hi = genes.Length;
        while (lo < hi)
        {
            int mid = (lo + hi) >>> 1;
            if (genes[mid].Location.Start < start) lo = mid + 1; else hi = mid;
        }
        return lo;
    }

    // ------------------------------------------------------------------ transcript model

    private sealed class TxModel
    {
        public required bool Forward;
        public required GenomicInterval[] Exons;   // transcript order
        public required int[] Offsets;             // transcript index of each exon's first base
        public required int Length;
        public byte[]? Seq;                         // upper-case spliced sequence, transcript orientation
        public int Cs = -1, Ce = -1;                // CDS [Cs, Ce) in transcript coordinates
        public bool Coding => Cs >= 0 && Ce > Cs;
        public required GenomicInterval Span;

        public int TxIndex(long g)
        {
            for (int i = 0; i < Exons.Length; i++)
            {
                var e = Exons[i];
                if (g >= e.Start && g < e.End)
                    return Offsets[i] + (int)(Forward ? g - e.Start : e.End - 1 - g);
            }
            return -1;
        }
    }

    private TxModel Model(Gene gene, Transcript t)
    {
        if (_models.TryGetValue(t, out var m)) return m;
        var strand = t.Exons[0].Strand != Strand.Unknown ? t.Exons[0].Strand : gene.Strand;
        bool fwd = strand != Strand.Reverse;
        var exons = fwd ? t.Exons.ToArray() : t.Exons.Reverse().ToArray();
        var offsets = new int[exons.Length];
        int len = 0;
        for (int i = 0; i < exons.Length; i++) { offsets[i] = len; len += (int)exons[i].Length; }
        m = new TxModel { Forward = fwd, Exons = exons, Offsets = offsets, Length = len, Span = t.Span };
        if (_genome.TryGetChromosome(t.Span.Chromosome, out var chrom) && chrom.Sequence is { } seq && t.Span.End <= seq.Length)
        {
            var buf = new byte[len];
            var src = seq.Span;
            foreach (var (e, o) in exons.Zip(offsets))
            {
                var part = src[(int)e.Start..(int)e.End];
                if (fwd) part.CopyTo(buf.AsSpan(o));
                else for (int k = 0; k < part.Length; k++) buf[o + k] = Complement(part[part.Length - 1 - k]);
            }
            for (int k = 0; k < buf.Length; k++) buf[k] = (byte)char.ToUpperInvariant((char)buf[k]);
            m.Seq = buf;
        }
        if (t.CodingRegion is { } cr && cr.Length >= 3)
        {
            int a = m.TxIndex(fwd ? cr.Start : cr.End - 1);
            int b = m.TxIndex(fwd ? cr.End - 1 : cr.Start);
            if (a >= 0 && b >= a) { m.Cs = a; m.Ce = b + 1; }
        }
        _models[t] = m;
        return m;
    }

    private static byte Complement(byte b) => (char)(b & 0xDF) switch
    {
        'A' => (byte)'T', 'C' => (byte)'G', 'G' => (byte)'C', 'T' => (byte)'A', 'U' => (byte)'A', _ => (byte)'N',
    };

    private static string ReverseComplement(string s)
    {
        var c = new char[s.Length];
        for (int i = 0; i < s.Length; i++) c[s.Length - 1 - i] = (char)Complement((byte)s[i]);
        return new string(c);
    }

    // ------------------------------------------------------------------ evaluation

    private VariantEffect? Evaluate(string chrom, long s, long e, string refG, string altG, string allele, Gene gene, Transcript t)
    {
        var tx = Model(gene, t);
        bool ins = s == e;
        long spanS = tx.Span.Start, spanE = tx.Span.End;
        bool inside = ins ? s > spanS && s < spanE : s < spanE && e > spanS;
        if (!inside)
        {
            bool left = ins ? s <= spanS : e <= spanS;
            long dist = left ? spanS - (ins ? s : e) : s - spanE;
            bool upstream = left == tx.Forward;
            if (dist > (upstream ? Options.UpstreamDistance : Options.DownstreamDistance)) return null;
            return new VariantEffect
            {
                Allele = allele,
                Consequences = [upstream ? Consequence.UpstreamGeneVariant : Consequence.DownstreamGeneVariant],
                GeneId = gene.Id, GeneName = gene.Name, TranscriptId = t.Id, Distance = dist,
            };
        }

        var cons = new SortedSet<Consequence>();
        bool Hit(long a, long b) => ins ? s > a && s < b : s < b && e > a;
        bool HitInclusive(long a, long b) => ins ? s >= a && s <= b : s < b && e > a;

        // Splice sites and introns (genomic order).
        var genomicExons = t.Exons;
        int exonBases = Options.SpliceRegionExonBases, intronBases = Options.SpliceRegionIntronBases;
        for (int i = 1; i < genomicExons.Count; i++)
        {
            long iS = genomicExons[i - 1].End, iE = genomicExons[i].Start;
            if (iE - iS < 1) continue;
            var (donorA, donorB) = tx.Forward ? (iS, iS + 2) : (iE - 2, iE);
            var (accA, accB) = tx.Forward ? (iE - 2, iE) : (iS, iS + 2);
            bool site = false;
            if (Hit(Math.Max(iS, donorA), Math.Min(iE, donorB))) { cons.Add(Consequence.SpliceDonorVariant); site = true; }
            if (Hit(Math.Max(iS, accA), Math.Min(iE, accB))) { cons.Add(Consequence.SpliceAcceptorVariant); site = true; }
            if (!site && (HitInclusive(iS - exonBases, iS) || HitInclusive(iS + 2, Math.Min(iE, iS + intronBases))
                || HitInclusive(Math.Max(iS, iE - intronBases), iE - 2) || HitInclusive(iE, iE + exonBases)))
                cons.Add(Consequence.SpliceRegionVariant);
            if (iE - iS > 4 && Hit(iS + 2, iE - 2)) cons.Add(Consequence.IntronVariant);
        }

        bool exonic = genomicExons.Any(x => Hit(x.Start, x.End));
        string refT = tx.Forward ? refG : ReverseComplement(refG);
        string altT = tx.Forward ? altG : ReverseComplement(altG);
        var (contiguous, t0, t1) = MapToTx(tx, s, e);

        int? cdsPos = null, protPos = null;
        string? hgvsP = null, aminoAcids = null, codons = null;
        if (exonic)
        {
            if (!tx.Coding) cons.Add(Consequence.NonCodingTranscriptExonVariant);
            else
            {
                bool hit5 = false, hit3 = false, hitCds = false;
                if (ins)
                {
                    if (contiguous)
                    {
                        if (t0 <= tx.Cs) hit5 = true; else if (t0 >= tx.Ce) hit3 = true; else hitCds = true;
                    }
                }
                else
                {
                    foreach (var x in genomicExons)
                    {
                        long a = Math.Max(s, x.Start), b = Math.Min(e, x.End);
                        if (a >= b) continue;
                        int ta = tx.TxIndex(tx.Forward ? a : b - 1), tb = tx.TxIndex(tx.Forward ? b - 1 : a) + 1;
                        if (ta < tx.Cs) hit5 = true;
                        if (tb > tx.Ce) hit3 = true;
                        if (ta < tx.Ce && tb > tx.Cs) hitCds = true;
                    }
                }
                if (hit5) cons.Add(Consequence.FivePrimeUtrVariant);
                if (hit3) cons.Add(Consequence.ThreePrimeUtrVariant);
                if (hitCds)
                {
                    bool fullyCoding = contiguous && !hit5 && !hit3 && tx.Seq is not null;
                    if (fullyCoding)
                        Coding(tx, t0, t1, refT, altT, cons, out cdsPos, out protPos, out hgvsP, out aminoAcids, out codons);
                    else cons.Add(Consequence.CodingSequenceVariant);
                }
            }
        }
        if (cons.Count == 0) cons.Add(tx.Coding ? Consequence.CodingSequenceVariant : Consequence.NonCodingTranscriptExonVariant);

        string? hgvsC = contiguous && tx.Seq is not null
            ? HgvsFromTx(tx, t0, t1, refT, altT)
            : HgvsFromGenomic(tx, s, e, refT, altT);

        return new VariantEffect
        {
            Allele = allele,
            Consequences = cons.ToArray(),
            GeneId = gene.Id,
            GeneName = gene.Name,
            TranscriptId = t.Id,
            HgvsC = hgvsC,
            HgvsP = hgvsP,
            CdsPosition = cdsPos,
            ProteinPosition = protPos,
            AminoAcids = aminoAcids,
            Codons = codons,
        };
    }

    /// <summary>Maps a genomic change to transcript coordinates; contiguous when every affected base is exonic and adjacent in the transcript.</summary>
    private static (bool Contiguous, int T0, int T1) MapToTx(TxModel tx, long s, long e)
    {
        if (s == e)
        {
            int left = tx.TxIndex(tx.Forward ? s - 1 : s), right = tx.TxIndex(tx.Forward ? s : s - 1);
            if (left >= 0 && right == left + 1) return (true, right, right);
            return (false, -1, -1);
        }
        int a = tx.TxIndex(tx.Forward ? s : e - 1), b = tx.TxIndex(tx.Forward ? e - 1 : s);
        if (a >= 0 && b >= a && b - a + 1 == e - s) return (true, a, b + 1);
        return (false, -1, -1);
    }

    private void Coding(TxModel tx, int t0, int t1, string refT, string altT, SortedSet<Consequence> cons,
        out int? cdsPos, out int? protPos, out string? hgvsP, out string? aminoAcids, out string? codons)
    {
        var code = Options.Code;
        var seq = tx.Seq!;
        int cs = tx.Cs, cdsLen = tx.Ce - tx.Cs;
        int c0 = t0 - cs, refLen = t1 - t0, diff = altT.Length - refLen;
        bool ins = refLen == 0;
        var cdsRef = seq.AsSpan(cs);
        var mut = new byte[cdsRef.Length + diff];
        cdsRef[..c0].CopyTo(mut);
        Encoding.ASCII.GetBytes(altT, mut.AsSpan(c0));
        cdsRef[(c0 + refLen)..].CopyTo(mut.AsSpan(c0 + altT.Length));

        string refProt = Translate(cdsRef[..cdsLen], code, toStop: true);
        string mutProt = Translate(mut, code, toStop: true);
        cdsPos = c0 + 1;
        protPos = c0 / 3 + 1;

        // Codons (VEP style: changed bases upper case).
        codons = null;
        int ca = c0 / 3 * 3;
        int cbRef = Math.Max(ca + 3, (c0 + refLen + 2) / 3 * 3);
        int cbMut = cbRef + diff;
        if (cbRef <= cdsRef.Length && cbMut <= mut.Length && cbMut >= ca)
            codons = Case(cdsRef[ca..cbRef], c0 - ca, refLen) + "/" + Case(mut.AsSpan(ca, cbMut - ca), c0 - ca, altT.Length);

        hgvsP = null;
        aminoAcids = null;
        int refStop = refProt.IndexOf('*');
        int mutStop = mutProt.IndexOf('*');

        // Initiator codon.
        bool touchesStart = ins ? c0 is > 0 and < 3 : c0 < 3;
        if (touchesStart && cdsLen >= 3 && !mut.AsSpan(0, 3).SequenceEqual(cdsRef[..3]))
        {
            int idx = GeneticCode.Index(mut.AsSpan(0, 3));
            if (code.IsStartIndex(idx) && diff % 3 == 0)
            {
                cons.Add(Consequence.StartRetainedVariant);
                hgvsP = $"p.{Three(refProt[0])}1=";
                aminoAcids = refProt[0].ToString();
                return;
            }
            cons.Add(Consequence.StartLost);
            hgvsP = $"p.{Three(refProt[0])}1?";
            aminoAcids = $"{refProt[0]}/{(idx < 0 ? 'X' : code.AminoAcidAt(idx))}";
            return;
        }

        int p = 0;
        while (p < refProt.Length && p < mutProt.Length && refProt[p] == mutProt[p]) p++;

        if (diff % 3 != 0)
        {
            cons.Add(Consequence.FrameshiftVariant);
            if (p >= refProt.Length) { hgvsP = "p.?"; return; }
            protPos = p + 1;
            if (p >= mutProt.Length) { hgvsP = "p.?"; return; }
            char altAa = mutProt[p];
            aminoAcids = $"{refProt[p]}/{altAa}";
            if (refProt[p] == '*')
            {
                hgvsP = $"p.Ter{p + 1}{Three(altAa)}extTer{(mutStop < 0 ? "?" : (mutStop - p + 1).ToString(CultureInfo.InvariantCulture))}";
                return;
            }
            if (altAa == '*') { hgvsP = $"p.{Three(refProt[p])}{p + 1}Ter"; return; }
            int stop = mutProt.IndexOf('*', p);
            hgvsP = $"p.{Three(refProt[p])}{p + 1}{Three(altAa)}fsTer{(stop < 0 ? "?" : (stop - p + 1).ToString(CultureInfo.InvariantCulture))}";
            return;
        }

        if (p == refProt.Length && p == mutProt.Length)
        {
            // Protein unchanged.
            int aa = Math.Min(c0 / 3, refProt.Length - 1);
            bool inStop = refStop >= 0 && aa == refStop;
            cons.Add(inStop ? Consequence.StopRetainedVariant : Consequence.SynonymousVariant);
            hgvsP = $"p.{Three(refProt[aa])}{aa + 1}=";
            aminoAcids = refProt[aa].ToString();
            protPos = aa + 1;
            return;
        }

        protPos = p + 1;
        // Stop lost: the reference stop is changed and no earlier stop appears.
        if (refStop >= 0 && p == refStop)
        {
            if (mutStop == refStop) { cons.Add(Consequence.StopRetainedVariant); hgvsP = $"p.Ter{refStop + 1}="; return; }
            cons.Add(Consequence.StopLost);
            if (diff > 0) cons.Add(Consequence.InframeInsertion);
            else if (diff < 0) cons.Add(Consequence.InframeDeletion);
            char altAa = p < mutProt.Length ? mutProt[p] : 'X';
            aminoAcids = $"*/{altAa}";
            hgvsP = $"p.Ter{refStop + 1}{Three(altAa)}extTer{(mutStop < 0 ? "?" : (mutStop - p + 1).ToString(CultureInfo.InvariantCulture))}";
            return;
        }

        // Premature stop.
        if (mutStop >= 0 && (refStop < 0 || mutStop < refStop + diff / 3))
        {
            cons.Add(Consequence.StopGained);
            if (diff > 0) cons.Add(Consequence.InframeInsertion);
            else if (diff < 0) cons.Add(Consequence.InframeDeletion);
            aminoAcids = $"{refProt[p]}/{mutProt[p]}";
            hgvsP = mutStop == p
                ? $"p.{Three(refProt[p])}{p + 1}Ter"
                : $"p.{Three(refProt[p])}{p + 1}delins{ThreeAll(mutProt.AsSpan(p, mutStop - p + 1))}";
            return;
        }

        // In-frame change: compare the proteins from both ends.
        int q = 0;
        while (q < refProt.Length - p && q < mutProt.Length - p && refProt[^(q + 1)] == mutProt[^(q + 1)]) q++;
        var r = refProt.AsSpan(p, refProt.Length - p - q);
        var i = mutProt.AsSpan(p, mutProt.Length - p - q);
        aminoAcids = $"{(r.IsEmpty ? "-" : r.ToString())}/{(i.IsEmpty ? "-" : i.ToString())}";
        if (diff > 0) cons.Add(Consequence.InframeInsertion);
        else if (diff < 0) cons.Add(Consequence.InframeDeletion);
        else cons.Add(r.Length == i.Length ? Consequence.MissenseVariant : Consequence.ProteinAlteringVariant);

        if (r.Length == 1 && i.Length == 1) hgvsP = $"p.{Three(r[0])}{p + 1}{Three(i[0])}";
        else if (r.IsEmpty)
        {
            int n = i.Length;
            if (p - n >= 0 && refProt.AsSpan(p - n, n).SequenceEqual(i))
                hgvsP = n == 1 ? $"p.{Three(refProt[p - 1])}{p}dup" : $"p.{Three(refProt[p - n])}{p - n + 1}_{Three(refProt[p - 1])}{p}dup";
            else if (p >= 1)
                hgvsP = $"p.{Three(refProt[p - 1])}{p}_{Three(refProt[p])}{p + 1}ins{ThreeAll(i)}";
        }
        else if (i.IsEmpty)
            hgvsP = r.Length == 1 ? $"p.{Three(r[0])}{p + 1}del" : $"p.{Three(r[0])}{p + 1}_{Three(r[^1])}{p + r.Length}del";
        else
            hgvsP = r.Length == 1 ? $"p.{Three(r[0])}{p + 1}delins{ThreeAll(i)}" : $"p.{Three(r[0])}{p + 1}_{Three(r[^1])}{p + r.Length}delins{ThreeAll(i)}";
    }

    private static string Case(ReadOnlySpan<byte> codon, int changedStart, int changedLength)
    {
        var c = new char[codon.Length];
        for (int k = 0; k < c.Length; k++)
        {
            char ch = (char)codon[k];
            c[k] = k >= changedStart && k < changedStart + changedLength ? char.ToUpperInvariant(ch) : char.ToLowerInvariant(ch);
        }
        return c.Length == 0 ? "-" : new string(c);
    }

    internal static string Translate(ReadOnlySpan<byte> nt, GeneticCode code, bool toStop)
    {
        var sb = new StringBuilder(nt.Length / 3 + 1);
        for (int k = 0; k + 3 <= nt.Length; k += 3)
        {
            int idx = GeneticCode.Index(nt.Slice(k, 3));
            char aa = idx < 0 ? 'X' : code.AminoAcidAt(idx);
            sb.Append(aa);
            if (toStop && aa == '*') break;
        }
        return sb.ToString();
    }

    private static readonly Dictionary<char, string> ThreeLetter = new()
    {
        ['A'] = "Ala", ['R'] = "Arg", ['N'] = "Asn", ['D'] = "Asp", ['C'] = "Cys", ['Q'] = "Gln", ['E'] = "Glu",
        ['G'] = "Gly", ['H'] = "His", ['I'] = "Ile", ['L'] = "Leu", ['K'] = "Lys", ['M'] = "Met", ['F'] = "Phe",
        ['P'] = "Pro", ['S'] = "Ser", ['T'] = "Thr", ['W'] = "Trp", ['Y'] = "Tyr", ['V'] = "Val", ['*'] = "Ter",
        ['U'] = "Sec", ['O'] = "Pyl", ['X'] = "Xaa",
    };

    /// <summary>Three-letter amino-acid code (<c>Ter</c> for stop).</summary>
    public static string Three(char aa) => ThreeLetter.TryGetValue(char.ToUpperInvariant(aa), out var s) ? s : "Xaa";

    private static string ThreeAll(ReadOnlySpan<char> aas)
    {
        var sb = new StringBuilder(aas.Length * 3);
        foreach (char c in aas) sb.Append(Three(c));
        return sb.ToString();
    }

    // ------------------------------------------------------------------ HGVS c./n.

    private static string CFromTx(TxModel tx, int t)
    {
        if (!tx.Coding) return (t + 1).ToString(CultureInfo.InvariantCulture);
        if (t < tx.Cs) return "-" + (tx.Cs - t).ToString(CultureInfo.InvariantCulture);
        if (t >= tx.Ce) return "*" + (t - tx.Ce + 1).ToString(CultureInfo.InvariantCulture);
        return (t - tx.Cs + 1).ToString(CultureInfo.InvariantCulture);
    }

    /// <summary>HGVS position of one genomic base: exonic, or intronic as <c>N+d</c>/<c>N-d</c>; null outside the transcript.</summary>
    private static string? CPos(TxModel tx, long g)
    {
        int t = tx.TxIndex(g);
        if (t >= 0) return CFromTx(tx, t);
        for (int i = 0; i + 1 < tx.Exons.Length; i++)
        {
            var a = tx.Exons[i];
            var b = tx.Exons[i + 1];
            long dPrev, dNext;
            if (tx.Forward)
            {
                if (g < a.End || g >= b.Start) continue;
                dPrev = g - (a.End - 1); dNext = b.Start - g;
            }
            else
            {
                if (g >= a.Start || g < b.End) continue;
                dPrev = a.Start - g; dNext = g - (b.End - 1);
            }
            int lastOfPrev = tx.Offsets[i] + (int)a.Length - 1, firstOfNext = tx.Offsets[i + 1];
            return dPrev <= dNext
                ? CFromTx(tx, lastOfPrev) + "+" + dPrev.ToString(CultureInfo.InvariantCulture)
                : CFromTx(tx, firstOfNext) + "-" + dNext.ToString(CultureInfo.InvariantCulture);
        }
        return null;
    }

    private static string Prefix(TxModel tx) => tx.Coding ? "c." : "n.";

    private static string HgvsFromTx(TxModel tx, int t0, int t1, string refT, string altT)
    {
        var seq = tx.Seq!;
        string pre = Prefix(tx);
        if (refT.Length == 1 && altT.Length == 1) return $"{pre}{CFromTx(tx, t0)}{refT}>{altT}";
        if (refT.Length == 0)
        {
            // Shift the insertion 3' in the transcript, then check for a duplication.
            var ins = altT.ToCharArray();
            while (t0 < seq.Length && seq[t0] == ins[0])
            {
                char first = ins[0];
                Array.Copy(ins, 1, ins, 0, ins.Length - 1);
                ins[^1] = first;
                t0++;
            }
            int n = ins.Length;
            string insStr = new(ins);
            if (t0 - n >= 0 && Encoding.ASCII.GetString(seq, t0 - n, n) == insStr)
                return n == 1 ? $"{pre}{CFromTx(tx, t0 - 1)}dup" : $"{pre}{CFromTx(tx, t0 - n)}_{CFromTx(tx, t0 - 1)}dup";
            if (t0 == 0 || t0 >= seq.Length) return $"{pre}{CFromTx(tx, Math.Max(0, t0 - 1))}ins{insStr}";
            return $"{pre}{CFromTx(tx, t0 - 1)}_{CFromTx(tx, t0)}ins{insStr}";
        }
        if (altT.Length == 0)
        {
            while (t1 < seq.Length && seq[t0] == seq[t1]) { t0++; t1++; }
            return t1 - t0 == 1 ? $"{pre}{CFromTx(tx, t0)}del" : $"{pre}{CFromTx(tx, t0)}_{CFromTx(tx, t1 - 1)}del";
        }
        return t1 - t0 == 1 ? $"{pre}{CFromTx(tx, t0)}delins{altT}" : $"{pre}{CFromTx(tx, t0)}_{CFromTx(tx, t1 - 1)}delins{altT}";
    }

    private static string? HgvsFromGenomic(TxModel tx, long s, long e, string refT, string altT)
    {
        string pre = Prefix(tx);
        if (s == e)
        {
            string? l = CPos(tx, tx.Forward ? s - 1 : s), r = CPos(tx, tx.Forward ? s : s - 1);
            return l is null || r is null ? null : $"{pre}{l}_{r}ins{altT}";
        }
        string? first = CPos(tx, tx.Forward ? s : e - 1), last = CPos(tx, tx.Forward ? e - 1 : s);
        if (first is null || last is null) return null;
        string range = e - s == 1 ? first : $"{first}_{last}";
        if (refT.Length == 1 && altT.Length == 1) return $"{pre}{first}{refT}>{altT}";
        if (altT.Length == 0) return $"{pre}{range}del";
        return $"{pre}{range}delins{altT}";
    }
}
