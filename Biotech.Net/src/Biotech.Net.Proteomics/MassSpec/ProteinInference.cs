using Biotech.Net.Core;

namespace Biotech.Net.Proteomics.MassSpec;

/// <summary>A protein group reported by parsimonious inference.</summary>
/// <param name="ProteinIds">Indistinguishable proteins (identical peptide evidence), sorted; the first is the representative.</param>
/// <param name="Peptides">All peptides (by unmodified sequence) assigned to the group.</param>
/// <param name="UniquePeptides">Peptides that map to no protein outside this group.</param>
public sealed record ProteinGroup(IReadOnlyList<string> ProteinIds, IReadOnlyList<string> Peptides, IReadOnlyList<string> UniquePeptides)
{
    /// <summary>Representative accession.</summary>
    public string Representative => ProteinIds[0];
}

/// <summary>
/// Protein inference by parsimony (Occam's razor, as in IDPicker / Nesvizhskii &amp; Aebersold 2005): proteins with
/// identical peptide sets are merged into groups, then a greedy minimum set cover picks the fewest groups that
/// explain every peptide (ties: more peptides, then accession order). Subsumed proteins are not reported.
/// </summary>
public static class ProteinInference
{
    /// <summary>Infers protein groups from a peptide → proteins map.</summary>
    public static List<ProteinGroup> Parsimony(IReadOnlyDictionary<string, IReadOnlyCollection<string>> peptideToProteins)
    {
        ArgumentNullException.ThrowIfNull(peptideToProteins);
        var proteinPeptides = new Dictionary<string, SortedSet<string>>(StringComparer.Ordinal);
        foreach (var (pep, prots) in peptideToProteins)
            foreach (var prot in prots)
            {
                if (!proteinPeptides.TryGetValue(prot, out var set)) proteinPeptides[prot] = set = new(StringComparer.Ordinal);
                set.Add(pep);
            }

        // Merge indistinguishable proteins.
        var byEvidence = new Dictionary<string, (List<string> Proteins, SortedSet<string> Peptides)>(StringComparer.Ordinal);
        foreach (var (prot, peps) in proteinPeptides)
        {
            string key = string.Join('\u0001', peps);
            if (!byEvidence.TryGetValue(key, out var g)) byEvidence[key] = g = ([], peps);
            g.Proteins.Add(prot);
        }
        var candidates = byEvidence.Values.Select(g => { g.Proteins.Sort(StringComparer.Ordinal); return g; })
            .OrderBy(g => g.Proteins[0], StringComparer.Ordinal).ToList();

        var uncovered = new HashSet<string>(peptideToProteins.Keys, StringComparer.Ordinal);
        var chosen = new List<(List<string> Proteins, SortedSet<string> Peptides)>();
        while (uncovered.Count > 0)
        {
            int bestIdx = -1, bestGain = 0;
            for (int i = 0; i < candidates.Count; i++)
            {
                int gain = 0;
                foreach (var p in candidates[i].Peptides) if (uncovered.Contains(p)) gain++;
                if (gain > bestGain || (gain == bestGain && gain > 0 && candidates[i].Peptides.Count > candidates[bestIdx].Peptides.Count))
                {
                    bestGain = gain; bestIdx = i;
                }
            }
            if (bestIdx < 0) break;
            var g = candidates[bestIdx];
            chosen.Add(g);
            candidates.RemoveAt(bestIdx);
            foreach (var p in g.Peptides) uncovered.Remove(p);
        }

        var result = new List<ProteinGroup>(chosen.Count);
        foreach (var (prots, peps) in chosen)
        {
            var groupSet = new HashSet<string>(prots, StringComparer.Ordinal);
            var unique = new List<string>();
            foreach (var p in peps)
            {
                bool only = true;
                foreach (var owner in peptideToProteins[p]) if (!groupSet.Contains(owner)) { only = false; break; }
                if (only) unique.Add(p);
            }
            result.Add(new(prots, peps.ToList(), unique));
        }
        return result;
    }

    /// <summary>Infers protein groups from accepted target PSMs (peptides keyed by unmodified sequence).</summary>
    public static List<ProteinGroup> Parsimony(IEnumerable<PeptideSpectrumMatch> psms)
    {
        ArgumentNullException.ThrowIfNull(psms);
        var map = new Dictionary<string, HashSet<string>>(StringComparer.Ordinal);
        foreach (var psm in psms)
        {
            if (psm.IsDecoy) continue;
            if (!map.TryGetValue(psm.Peptide.Sequence, out var set)) map[psm.Peptide.Sequence] = set = new(StringComparer.Ordinal);
            foreach (var p in psm.ProteinIds) set.Add(p);
        }
        return Parsimony(map.ToDictionary(kv => kv.Key, kv => (IReadOnlyCollection<string>)kv.Value, StringComparer.Ordinal));
    }
}

/// <summary>Label-free abundance estimates for one protein.</summary>
/// <param name="ProteinId">Accession.</param>
/// <param name="Length">Residue count.</param>
/// <param name="SpectralCount">Number of PSMs/observations assigned.</param>
/// <param name="Nsaf">Normalised spectral abundance factor (SpC/L divided by Σ SpC/L), Zybailov et al. 2006.</param>
/// <param name="Intensity">Summed precursor intensity.</param>
/// <param name="ObservablePeptides">Theoretical fully-cleaved peptides of 6–30 residues.</param>
/// <param name="Ibaq">Intensity divided by observable peptides (Schwanhäusser et al. 2011); NaN when none are observable.</param>
public readonly record struct ProteinQuantity(string ProteinId, int Length, int SpectralCount, double Nsaf, double Intensity, int ObservablePeptides, double Ibaq);

/// <summary>Label-free quantification: spectral counting, NSAF and iBAQ.</summary>
public static class LabelFreeQuantification
{
    /// <summary>
    /// Quantifies proteins from observations (protein accession, precursor intensity). Each observation adds one
    /// spectral count and its intensity (NaN intensities count but add nothing) to the named protein.
    /// </summary>
    public static List<ProteinQuantity> Quantify(IEnumerable<(string ProteinId, double Intensity)> observations, IEnumerable<Protein> proteins, Enzyme? enzyme = null)
    {
        ArgumentNullException.ThrowIfNull(observations);
        ArgumentNullException.ThrowIfNull(proteins);
        var byId = new Dictionary<string, Protein>(StringComparer.Ordinal);
        foreach (var p in proteins) byId[p.Id] = p;
        var counts = new Dictionary<string, (int Count, double Intensity)>(StringComparer.Ordinal);
        foreach (var (id, intensity) in observations)
        {
            if (!byId.ContainsKey(id)) continue;
            counts.TryGetValue(id, out var c);
            counts[id] = (c.Count + 1, c.Intensity + (double.IsNaN(intensity) ? 0 : intensity));
        }

        double saf = 0;
        foreach (var (id, c) in counts) saf += (double)c.Count / byId[id].Length;
        var result = new List<ProteinQuantity>(counts.Count);
        foreach (var (id, c) in counts.OrderBy(kv => kv.Key, StringComparer.Ordinal))
        {
            var prot = byId[id];
            int obs = Digestion.ObservablePeptideCount(prot.Sequence.Span, enzyme);
            result.Add(new(id, prot.Length, c.Count, (double)c.Count / prot.Length / saf, c.Intensity, obs, obs > 0 ? c.Intensity / obs : double.NaN));
        }
        return result;
    }

    /// <summary>
    /// Quantifies from accepted target PSMs. A PSM shared by several proteins is credited to each of them
    /// (use <see cref="ProteinInference"/> groups' representatives to avoid double counting).
    /// </summary>
    public static List<ProteinQuantity> Quantify(IEnumerable<PeptideSpectrumMatch> psms, IEnumerable<Protein> proteins, Enzyme? enzyme = null)
    {
        ArgumentNullException.ThrowIfNull(psms);
        var obs = new List<(string, double)>();
        foreach (var psm in psms)
        {
            if (psm.IsDecoy) continue;
            foreach (var id in psm.ProteinIds) obs.Add((id, psm.Spectrum.PrecursorIntensity));
        }
        return Quantify(obs, proteins, enzyme);
    }

    /// <summary>NSAF for parallel arrays of spectral counts and protein lengths.</summary>
    public static double[] Nsaf(ReadOnlySpan<double> spectralCounts, ReadOnlySpan<int> lengths)
    {
        if (spectralCounts.Length != lengths.Length) throw new ArgumentException("Counts and lengths differ in length.");
        var r = new double[spectralCounts.Length];
        double total = 0;
        for (int i = 0; i < r.Length; i++) { r[i] = spectralCounts[i] / lengths[i]; total += r[i]; }
        if (total > 0) for (int i = 0; i < r.Length; i++) r[i] /= total;
        return r;
    }
}
