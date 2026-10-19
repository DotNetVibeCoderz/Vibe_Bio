namespace Biotech.Net.WebGallery.Playground;

/// <summary>A Playground example: taken from the docs tutorials (docs/en/tutorials), trimmed to run in the browser.</summary>
public sealed record Preset(string Id, string Title, string Module, string Tutorial, string Code);

public static class Presets
{
    public static readonly Preset[] All =
    [
        new("sequences", "Sequences: translate and reverse-complement", "Core · Sequences", "tutorials/sequences.md", """
            using Biotech.Net.Core;

            var dna = Sequence.Dna("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", id: "demo");
            Console.WriteLine(dna.Length);                  // 39
            Console.WriteLine(dna.GcContent());             // 0.564...
            Console.WriteLine(dna.ReverseComplement());     // CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT
            Console.WriteLine(dna.Transcribe());            // AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG
            Console.WriteLine(dna.Translate());             // MAIVMGR*KGAR*
            Console.WriteLine(dna.Translate(toStop: true)); // MAIVMGR

            // Other genetic codes: NCBI tables 1, 2, 4 and 11
            Console.WriteLine(dna.Translate(GeneticCode.VertebrateMitochondrial));

            // Alphabet detection and subsequences (0-based start, length)
            var unknown = Sequence.Parse("MKTAYIAKQR");
            Console.WriteLine(unknown.Type);                // Protein
            var cds = dna.Subsequence(0, 21);
            Console.WriteLine(cds.IndexOf("GGG"));          // 14
            Console.WriteLine($"kernels: {BiotechInfo.Backend}");
            """),

        new("alignment", "Pairwise alignment (BLOSUM62, local)", "Alignments", "tutorials/alignments.md", """
            using Biotech.Net.Alignments;
            using Biotech.Net.Native;

            var aligner = new PairwiseAligner
            {
                Mode = AlignmentMode.Local,
                SubstitutionMatrix = SubstitutionMatrix.Blosum62,
                OpenGapScore = -11,
                ExtendGapScore = -1,
            };
            var aln = aligner.Align("HEAGAWGHEE", "PAWHEAE");
            Console.WriteLine(aln.Score);           // 17
            Console.WriteLine(aln.Format());        // three-line view with match markers
            Console.WriteLine(aln.Cigar());
            Console.WriteLine($"identity {aln.Identity:P1}, gaps {aln.Gaps}");

            // Nucleotide global alignment with a simple match/mismatch scheme
            var dnaAligner = new PairwiseAligner { MatchScore = 2, MismatchScore = -3, OpenGapScore = -5, ExtendGapScore = -2 };
            Console.WriteLine(dnaAligner.Score("GATTACAGATTACA", "GATCACAGATACA"));
            """),

        new("msa", "Multiple alignment and a neighbour-joining tree", "Alignments", "tutorials/alignments.md", """
            using Biotech.Net.Alignments;
            using Biotech.Net.Alignments.Phylogenetics;

            Sequence[] globins =
            [
                Sequence.Protein("VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF", "HBA_HUMAN"),
                Sequence.Protein("VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESF", "HBB_HUMAN"),
                Sequence.Protein("VLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRF", "MYG_HUMAN"),
                Sequence.Protein("VLSAADKTNVKAAWSKVGGHAGEYGAEALERMFLGFPTTKTYFPHF", "HBA_MOUSE"),
            ];

            var msa = new ProgressiveAligner { SubstitutionMatrix = SubstitutionMatrix.Blosum62 }.Align(globins);
            Console.WriteLine(msa.ToClustal());
            Console.WriteLine($"consensus: {msa.Consensus()}");

            var distances = DistanceMatrix.FromAlignment(msa, DistanceModel.PDistance);
            var tree = TreeBuilder.NeighborJoining(distances).RootAtMidpoint();
            Console.WriteLine(tree.ToNewick(decimals: 3));
            """),

        new("protparam", "ProtParam: insulin B chain", "Proteomics", "tutorials/proteomics.md", """
            using Biotech.Net.Proteomics;

            // Human insulin B chain
            var analysis = new ProteinAnalysis("FVNQHLCGSHLVEALYLVCGERGFFYTPKT");
            Console.WriteLine($"MW        {analysis.MolecularWeight():F2} Da");
            Console.WriteLine($"pI        {analysis.IsoelectricPoint():F2}");
            Console.WriteLine($"GRAVY     {analysis.Gravy():F3}");
            Console.WriteLine($"Instab.   {analysis.InstabilityIndex():F2} (stable: {analysis.IsStable})");
            Console.WriteLine($"Aliphatic {analysis.AliphaticIndex():F2}");
            Console.WriteLine($"Ext. coef {analysis.ExtinctionCoefficient().Reduced} M^-1 cm^-1");
            Console.WriteLine($"Half-life {analysis.HalfLife().Text}");

            // Sequence-based predictors
            Console.WriteLine(ChouFasman.Predict("FVNQHLCGSHLVEALYLVCGERGFFYTPKT"));
            var tm = TransmembranePredictor.Predict("MKTLLILAVVAAALAGSLLLVVFLLAGIAALLGYKRRSS");
            var signal = SignalPeptidePredictor.Predict("MKKTAIAIAVALAGFATVAQAAPKDNTWYTGAKLGWSQ");
            Console.WriteLine($"{tm.Count} TM segment(s); signal peptide: {signal.HasSignalPeptide} (cleavage after {signal.CleavageSite})");
            """),

        new("structure", "Crambin: DSSP, SASA, disulfides, 3-D view", "StructuralBiology", "tutorials/structural-biology.md", """
            using Biotech.Net.StructuralBiology;

            var structure = SampleStructures.Crambin();      // PDB 1CRN, embedded in the library
            var chain = structure["A"];
            Console.WriteLine(structure);
            Console.WriteLine($"Rg = {Geometry.RadiusOfGyration(structure.Atoms):F2} Å");

            foreach (var t in Geometry.BackboneTorsions(chain).Take(3))
                Console.WriteLine($"{t.Residue.Id} phi {t.Phi:F1} psi {t.Psi:F1} {t.Region}");

            DsspResult dssp = Dssp.Assign(structure);
            Console.WriteLine(dssp.ToSecondaryStructureString());

            var sasa = ShrakeRupley.PerResidue(structure.FirstModel);
            Console.WriteLine($"buried residues: {sasa.Count(r => r.IsBuried)}/{sasa.Count}");

            foreach (var ss in Disulfides.Find(structure))
                Console.WriteLine($"{ss.CysteineA.Id} - {ss.CysteineB.Id} {ss.Distance:F2} Å");

            Show(structure);                                 // interactive 3Dmol.js viewer
            """),

        new("chemistry", "Aspirin: descriptors, Lipinski, depiction", "Chemistry", "tutorials/chemistry.md", """
            using Biotech.Net.Chemistry;
            using Biotech.Net.Visualization;

            var aspirin = Molecule.FromSmiles("CC(=O)Oc1ccccc1C(=O)O", "aspirin");
            Console.WriteLine($"{aspirin.Name}: {aspirin.Formula}, {aspirin.AtomCount} heavy atoms");
            Console.WriteLine($"canonical SMILES: {aspirin.ToSmiles()}");

            DescriptorSet d = MolecularDescriptors.Compute(aspirin);
            Console.WriteLine($"MW {d.MolecularWeight:F2}  logP {d.LogP:F2}  TPSA {d.Tpsa:F1}  HBD {d.HBondDonors}  HBA {d.HBondAcceptors}  RotB {d.RotatableBonds}");
            Console.WriteLine(DrugLikeness.Lipinski(d));
            Console.WriteLine(DrugLikeness.Veber(d));
            Console.WriteLine($"QED {Qed.Compute(aspirin):F3}");
            Console.WriteLine($"Murcko scaffold: {MurckoScaffold.GetScaffoldSmiles(aspirin)}");

            Show(new MoleculeViewer(aspirin, 420, 300) { Theme = BioTheme.Dark });
            """),

        new("admet", "ADMET profile and PAINS alerts", "DrugDiscovery", "tutorials/drug-discovery.md", """
            using Biotech.Net.Chemistry;
            using Biotech.Net.DrugDiscovery;

            var mol = Molecule.FromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "ibuprofen");
            AdmetProfile admet = Admet.Profile(mol);
            Console.WriteLine(admet);                       // summary text
            Console.WriteLine($"ESOL logS {admet.EsolLogS:F2} ({admet.SolubilityClass}), BBB {admet.BloodBrainBarrier}, GI {admet.GiAbsorption}");
            Console.WriteLine($"hERG risk {admet.HergRisk:F2}, QED {admet.Qed:F2}");

            var rhodanine = Molecule.FromSmiles("O=C1NC(=S)SC1=Cc1ccccc1", "rhodanine-like");
            foreach (var alert in StructuralAlerts.Match(rhodanine)) Console.WriteLine($"{alert.Family}: {alert.Name}");

            var filter = new CompoundFilter().RequireLipinski().ExcludeAlerts(StructuralAlerts.Pains);
            Console.WriteLine(filter.Evaluate(rhodanine));
            Console.WriteLine(filter.Evaluate(mol));
            """),

        new("codons", "Codon optimisation for E. coli", "SyntheticBiology", "tutorials/synthetic-biology.md", """
            using Biotech.Net.SyntheticBiology.Codons;

            var options = new CodonOptimizationOptions
            {
                Target = CodonUsageTable.EColiK12,
                Method = CodonOptimizationMethod.MostFrequent,
                Constraints = SequenceConstraints.ForGeneSynthesis,
            };
            var result = CodonOptimizer.Optimize("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ", options);
            Console.WriteLine($"CAI {result.Cai:F3}, GC {result.GcContent:P1}, {result.RemainingViolations.Count} violations left");
            Console.WriteLine(result.Dna);

            // Re-code the gene for a human host, keeping the protein identical
            var humanOptions = new CodonOptimizationOptions
            {
                Target = CodonUsageTable.HomoSapiens,
                Constraints = SequenceConstraints.ForGeneSynthesis,
            };
            var recoded = CodonOptimizer.OptimizeCoding(result.Dna, humanOptions);
            Console.WriteLine($"codons changed: {recoded.CodonsChanged}");
            """),

        new("vep", "Variant effect on a toy gene", "Genomics", "tutorials/genomics.md", """
            using Biotech.Net.Genomics.Annotation;

            // A one-exon gene: ATG AAA CGC ... TAA on chr1 positions 10..39 (0-based)
            var chrSeq = Sequence.Dna("GGGGGGGGGG" + "ATGAAACGCATTAGCACCACCATTACCTAA" + "GGGGGGGGGG", "chr1");
            var genome = new Genome("toy").Add(Chromosome.FromSequence(chrSeq));
            var gene = new Gene("G1", new GenomicInterval("chr1", 10, 40, Strand.Forward), "TOY1");
            gene.Transcripts.Add(new Transcript("TOY1-201", "G1", [new GenomicInterval("chr1", 10, 40, Strand.Forward)])
                { CodingRegion = new GenomicInterval("chr1", 10, 40, Strand.Forward) });
            genome.Add(gene);

            var vep = new VariantEffectPredictor(genome);
            // Position 17 (1-based) is the first base of codon 3 (CGC, Arg); C>T gives TGC (Cys)
            foreach (var effect in vep.Predict("chr1", 17, "C", "T"))
                Console.WriteLine($"{effect.ConsequenceTerms} {effect.Impact} {effect.HgvsC} {effect.HgvsP}");
            """),

        new("volcano", "DESeq2-style test and a volcano plot", "Transcriptomics · Visualization", "tutorials/visualization.md", """
            using Biotech.Net.Transcriptomics;
            using Biotech.Net.Visualization;

            var rng = new Random(42);
            int genes = 400, samples = 6;
            var values = new double[genes * samples];
            for (int g = 0; g < genes; g++)
            {
                double mean = Math.Exp(2 + rng.NextDouble() * 6);
                double fold = g < 25 ? 4 : g < 50 ? 0.25 : 1;
                for (int s = 0; s < samples; s++)
                    values[g * samples + s] = Math.Round(mean * (s >= 3 ? fold : 1)
                        * Math.Exp(0.25 * (rng.NextDouble() + rng.NextDouble() + rng.NextDouble() - 1.5)));
            }
            var counts = new Dataset("counts",
                Enumerable.Range(1, genes).Select(i => $"gene{i:D4}").ToArray(),
                ["ctl1", "ctl2", "ctl3", "trt1", "trt2", "trt3"], values);
            counts.Annotate("condition", "ctl", "ctl", "ctl", "trt", "trt", "trt");

            var result = DESeq.Run(counts, "condition", "ctl", "trt");
            foreach (var row in result.Top(5))
                Console.WriteLine($"{row.Gene}  log2FC {row.Log2FoldChange:F2}  padj {row.PAdj:E2}");

            var volcano = new VolcanoPlot(result) { Title = "Treated vs control", LabelTop = 8, Theme = BioTheme.Dark };
            Console.WriteLine(volcano.Counts());   // (Up, Down, NotSignificant)
            Show(volcano);
            """),
    ];

    public static Preset Default => All[0];
}
