namespace Biotech.Net.StructuralBiology.Tests;

public class PdbIoTests
{
    [Fact]
    public void Crambin_HasExpectedComposition()
    {
        var s = SampleStructures.Crambin();
        Assert.Equal("1CRN", s.Id);
        Assert.Single(s.Models);
        Assert.Single(s.Chains);
        Assert.Equal(327, s.AtomCount);
        Assert.Equal(46, s["A"].Residues.Count);
        // Sequence from the PDB SEQRES / UniProt P01542.
        Assert.Equal("TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN", s["A"].ToSequence().ToString());
    }

    [Fact]
    public void Crambin_HeaderRecordsAreParsed()
    {
        var s = SampleStructures.Crambin();
        Assert.Equal("PLANT PROTEIN", s.Classification);
        Assert.StartsWith("WATER STRUCTURE OF A HYDROPHOBIC PROTEIN", s.Title);
        Assert.Contains("PENTAGON RINGS", s.Title);
        Assert.NotNull(s.Cell);
        Assert.Equal(40.960, s.Cell!.A, 3);
        Assert.Equal(90.77, s.Cell.Beta, 2);
        Assert.Equal("P 1 21 1", s.Cell.SpaceGroup);
        Assert.Equal(2, s.Cell.Z);
        Assert.Equal(2, s.Helices.Count);
        Assert.Equal((7, 19), (s.Helices[0].StartNumber, s.Helices[0].EndNumber));
        Assert.Equal("ILE", s.Helices[0].StartResidueName);
        Assert.Equal(2, s.Sheets.Count);
        Assert.Equal(-1, s.Sheets[1].Sense);
        Assert.Equal(3, s.Disulfides.Count);
        Assert.Equal((3, 40), (s.Disulfides[0].NumberA, s.Disulfides[0].NumberB));
        Assert.Equal(2.00, s.Disulfides[0].Length, 2);
    }

    [Fact]
    public void Crambin_FirstAtomFieldsAreExact()
    {
        var a = SampleStructures.Crambin().Atoms.First();
        Assert.Equal("N", a.Name);
        Assert.Equal("N", a.Element);
        Assert.Equal(new Vec3(17.047, 14.099, 3.625), a.Position);
        Assert.Equal(13.79, a.BFactor, 2);
        Assert.Equal(1.0, a.Occupancy);
        Assert.Equal("THR", a.ResidueName);
        Assert.Equal(1, a.ResidueNumber);
        Assert.Equal("A", a.ChainId);
        Assert.False(a.IsHetero);
        Assert.True(a.IsBackbone);
    }

    [Fact]
    public void Crambin_UnitCellVolume()
    {
        // Monoclinic cell: V = abc·sinβ.
        var c = SampleStructures.Crambin().Cell!;
        Assert.Equal(40.960 * 18.650 * 22.520 * Math.Sin(90.77 * Math.PI / 180), c.Volume, 6);
    }

    [Fact]
    public void Writer_RoundTripsCoordinatesAndHeader()
    {
        var s = SampleStructures.Crambin();
        string text = PdbWriter.ToPdb(s);
        var back = PdbReader.Parse(text, "1CRN");
        Assert.Equal(s.AtomCount, back.AtomCount);
        foreach (var (a, b) in s.Atoms.Zip(back.Atoms))
        {
            Assert.Equal(a.Name, b.Name);
            Assert.Equal(a.Element, b.Element);
            Assert.Equal(a.Position, b.Position);
            Assert.Equal(a.BFactor, b.BFactor, 2);
            Assert.Equal(a.ResidueNumber, b.ResidueNumber);
        }
        Assert.Equal(s.Title, back.Title);
        Assert.Equal(s.Helices, back.Helices);
        Assert.Equal(s.Sheets, back.Sheets);
        Assert.Equal(s.Disulfides.Select(d => (d.ChainA, d.NumberA, d.ChainB, d.NumberB)),
            back.Disulfides.Select(d => (d.ChainA, d.NumberA, d.ChainB, d.NumberB)));
        Assert.Equal(s.Cell, back.Cell);
    }

    [Fact]
    public void Writer_AtomLineHasPdbColumns()
    {
        var s = SampleStructures.Crambin();
        var ca = s["A"].Residues[0]["CA"];
        string line = PdbWriter.FormatAtom(ca, ca.Residue!, "A", 2);
        Assert.Equal(80, line.Length);
        Assert.Equal("ATOM      2  CA  THR A   1      16.967  12.784   4.338  1.00 10.80           C  ", line);
    }

    private const string AltLocPdb = """
        ATOM      1  N   SER A  10       1.000   2.000   3.000  1.00 10.00           N
        ATOM      2  CA ASER A  10       2.000   2.000   3.000  0.40 10.00           C
        ATOM      3  CA BSER A  10       2.100   2.100   3.100  0.60 11.00           C
        ATOM      4  C   SER A  10       3.000   2.000   3.000  1.00 10.00           C
        HETATM    5 ZN    ZN A 101      10.000  10.000  10.000  1.00 20.00          ZN2+
        HETATM    6  O   HOH A 201      12.000  10.000  10.000  1.00 30.00           O
        """;

    [Fact]
    public void AltLocs_HighestOccupancyIsSelected()
    {
        var s = PdbReader.Parse(AltLocPdb);
        var ser = s["A"].Residues[0];
        Assert.Equal(3, ser.Atoms.Count);
        var ca = ser["CA"];
        Assert.Equal('B', ca.AltLoc);
        Assert.Equal(0.60, ca.Occupancy);
        Assert.Single(ca.Alternates);
        Assert.Equal('A', ca.Alternates[0].AltLoc);
    }

    [Fact]
    public void AltLocs_AreWrittenBack()
    {
        var s = PdbReader.Parse(AltLocPdb);
        string text = PdbWriter.ToPdb(s);
        Assert.Contains(" CA ASER A  10", text);
        Assert.Contains(" CA BSER A  10", text);
        Assert.DoesNotContain(" CA ASER", PdbWriter.ToPdb(s, new PdbWriteOptions { IncludeAlternates = false }));
    }

    [Fact]
    public void Hetero_ChargeElementLigandAndWater()
    {
        var s = PdbReader.Parse(AltLocPdb);
        var zn = s.Atoms.Single(a => a.Element == "ZN");
        Assert.True(zn.IsHetero);
        Assert.Equal(2, zn.Charge);
        Assert.Equal("ZN", s.Ligands.Single().Name);
        Assert.True(s.Residues.Single(r => r.Name == "HOH").IsWater);
        Assert.Contains("ZN2+", PdbWriter.ToPdb(s));
    }

    [Fact]
    public void MissingElementColumns_AreInferredFromNameAlignment()
    {
        const string pdb = """
            ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00
            HETATM    2 CA    CA A   2       5.000   0.000   0.000  1.00  0.00
            ATOM      3 HG21 THR A   3       9.000   0.000   0.000  1.00  0.00
            HETATM    4 FE   HEM A   4       9.000   5.000   0.000  1.00  0.00
            """;
        var atoms = PdbReader.Parse(pdb).Atoms.ToList();
        Assert.Equal(["C", "CA", "H", "FE"], atoms.Select(a => a.Element));
    }

    [Fact]
    public void MultipleModels_AreSeparated()
    {
        const string pdb = """
            MODEL        1
            ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00  0.00           C
            ATOM      2  CA  GLY A   2       3.800   0.000   0.000  1.00  0.00           C
            ENDMDL
            MODEL        2
            ATOM      1  CA  GLY A   1       0.000   1.000   0.000  1.00  0.00           C
            ATOM      2  CA  GLY A   2       3.800   1.000   0.000  1.00  0.00           C
            ENDMDL
            END
            """;
        var s = PdbReader.Parse(pdb);
        Assert.Equal(2, s.Models.Count);
        Assert.Equal([1, 2], s.Models.Select(m => m.Number));
        Assert.Equal(1.0, s.Models[1].Atoms.First().Position.Y);
        string written = PdbWriter.ToPdb(s);
        Assert.Contains("MODEL        2", written);
        Assert.Equal(2, PdbReader.Parse(written).Models.Count);
    }

    [Fact]
    public void Conect_RecordsAreBidirectionalAndRoundTrip()
    {
        var s = TestData.Trypsin();
        Assert.NotEmpty(s.Connections);
        foreach (var (a, partners) in s.Connections)
            foreach (int b in partners) Assert.Contains(a, s.Connections[b]);
        var back = PdbReader.Parse(PdbWriter.ToPdb(s));
        Assert.Equal(s.Connections.Count, back.Connections.Count);
    }

    [Fact]
    public void Trypsin_LigandsAndAtoms()
    {
        var s = TestData.Trypsin();
        var ligands = s.Ligands.Select(r => r.Name).ToList();
        Assert.Contains("BEN", ligands);
        Assert.Contains("CA", ligands);
        Assert.Equal(9, s.Ligands.Single(r => r.Name == "BEN").Atoms.Count);
        Assert.Equal(1629, s.Atoms.Count(a => !a.IsHetero));
    }

    [Fact]
    public void TruncatedAtomLine_Throws()
    {
        Assert.Throws<InvalidDataException>(() => PdbReader.Parse("ATOM      1  CA  ALA A   1       0.000   0.000"));
    }

    [Fact]
    public void ChainToSequence_IsCoreProteinSequence()
    {
        var seq = SampleStructures.Crambin()["A"].ToSequence();
        Assert.Equal(Biotech.Net.Core.SequenceType.Protein, seq.Type);
        Assert.Equal(46, seq.Length);
        Assert.Equal("1CRN_A", seq.Id);
    }
}

public class MmCifTests
{
    [Fact]
    public void Crambin_MmCifMatchesPdb()
    {
        var cif = MmCifReader.Read(TestData.PathOf("1crn.cif"));
        var pdb = SampleStructures.Crambin();
        Assert.Equal(327, cif.AtomCount);
        foreach (var (a, b) in cif.Atoms.Zip(pdb.Atoms))
        {
            Assert.Equal(b.Name, a.Name);
            Assert.Equal(b.Element, a.Element);
            Assert.Equal(b.Position, a.Position);
            Assert.Equal(b.ResidueNumber, a.ResidueNumber);
            Assert.Equal(b.ChainId, a.ChainId);
        }
        Assert.Equal(pdb["A"].ToSequence().ToString(), cif["A"].ToSequence().ToString());
    }

    [Fact]
    public void Crambin_MmCifMetadata()
    {
        var cif = MmCifReader.Read(TestData.PathOf("1crn.cif"));
        Assert.Equal(40.960, cif.Cell!.A, 3);
        Assert.Equal("P 1 21 1", cif.Cell.SpaceGroup);
        Assert.Contains("PENTAGON RINGS", cif.Title, StringComparison.OrdinalIgnoreCase);
        Assert.Equal(3, cif.Disulfides.Count);
        Assert.Equal(2, cif.Helices.Count);
        Assert.Equal(7, cif.Helices[0].StartNumber);
    }

    [Fact]
    public void Tokenizer_HandlesQuotesSemicolonFieldsAndPrimes()
    {
        const string cif = """
            data_test
            _struct.title 'A title with spaces'
            _struct_keywords.pdbx_keywords
            ;Multi-line
            text field
            ;
            loop_
            _atom_site.group_PDB
            _atom_site.id
            _atom_site.type_symbol
            _atom_site.label_atom_id
            _atom_site.label_comp_id
            _atom_site.label_asym_id
            _atom_site.label_seq_id
            _atom_site.Cartn_x
            _atom_site.Cartn_y
            _atom_site.Cartn_z
            _atom_site.occupancy
            ATOM 1 O "O5'" DA B 1 1.5 2.5 3.5 1.00
            ATOM 2 C "C5'" DA B 1 2.5 2.5 3.5 ?
            HETATM 3 O O HOH W . 9.0 9.0 9.0 1.00 # trailing comment
            """;
        var s = MmCifReader.Parse(cif);
        Assert.Equal("test", s.Id);
        Assert.Equal("A title with spaces", s.Title);
        Assert.Equal("Multi-line\ntext field", s.Classification.Replace("\r", ""));
        var atoms = s.Atoms.ToList();
        Assert.Equal(3, atoms.Count);
        Assert.Equal("O5'", atoms[0].Name);
        Assert.Equal("C5'", atoms[1].Name);
        Assert.Equal(1.0, atoms[1].Occupancy);
        Assert.True(atoms[2].IsHetero);
        Assert.Equal(["B", "W"], s.Chains.Select(c => c.Id));
    }

    [Fact]
    public void GenericCifBlock_ExposesLoopsAndItems()
    {
        var blocks = MmCifReader.ParseCif(File.ReadAllText(TestData.PathOf("1crn.cif")));
        var block = Assert.Single(blocks);
        Assert.Equal("1CRN", block.Name);
        var loop = block.GetLoop("_atom_site");
        Assert.NotNull(loop);
        Assert.Equal(327, loop!.Rows.Count);
        Assert.Equal("40.960", block.Get("_cell.length_a"));
        Assert.Null(block.Get("_no.such_tag"));
    }

    [Fact]
    public void ToMmCif_RoundTrips()
    {
        var s = SampleStructures.Crambin();
        var back = MmCifReader.Parse(MmCifReader.ToMmCif(s));
        Assert.Equal(s.AtomCount, back.AtomCount);
        Assert.Equal(s.Atoms.Select(a => a.Position), back.Atoms.Select(a => a.Position));
    }
}
