using Biotech.Net.Chemistry;

namespace Biotech.Net.Chemistry.Tests;

public class MolFileTests
{
    private const string Ethanol = """
        ethanol
          handmade

          3  2  0  0  0  0  0  0  0  0999 V2000
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
            1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
            2.5981    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
          1  2  1  0
          2  3  1  0
        M  END
        """;

    [Fact]
    public void Reads_a_v2000_block()
    {
        var m = MolFile.Parse(Ethanol);
        Assert.Equal("ethanol", m.Name);
        Assert.Equal(3, m.AtomCount);
        Assert.Equal("C2H6O", m.Formula);
        Assert.Equal(2.5981, m.Coordinates![2].X, 4);
    }

    [Fact]
    public void Charges_isotopes_and_explicit_hydrogens()
    {
        // ammonium written with explicit H atoms and M  CHG; 13C via M  ISO
        var src = Molecule.FromSmiles("[NH4+].[13CH3]C(=O)[O-]");
        var block = MolFile.Write(src);
        Assert.Contains("M  CHG", block);
        Assert.Contains("M  ISO", block);
        var back = MolFile.Parse(block);
        Assert.Equal(src.ToSmiles(), back.ToSmiles());

        var withH = MolFile.Write(Molecule.FromSmiles("C[NH3+]").AddHydrogens());
        var readH = MolFile.Parse(withH);
        Assert.Equal(2, readH.AtomCount); // explicit H folded away
        Assert.Equal("CH6N+", readH.Formula);
        var keepH = MolFile.Parse(withH, new MolReaderOptions { RemoveHydrogens = false });
        Assert.Equal(8, keepH.AtomCount);
    }

    [Theory]
    [InlineData("CC(=O)Oc1ccccc1C(=O)O")]
    [InlineData("Cn1cnc2c1c(=O)n(C)c(=O)n2C")]
    [InlineData("c1cc[nH]c1")]
    [InlineData("C[S+](C)[O-]")]
    [InlineData("[Na+].[Cl-]")]
    [InlineData("C#N")]
    public void Write_then_read_preserves_the_molecule(string smiles)
    {
        var m = Molecule.FromSmiles(smiles);
        Assert.Equal(m.ToSmiles(), MolFile.Parse(MolFile.Write(m)).ToSmiles());
    }

    [Fact]
    public void Writes_kekule_bonds_and_valence_for_unusual_hydrogen_counts()
    {
        var block = MolFile.Write(Molecule.FromSmiles("c1ccccc1"));
        var bondLines = block.Split('\n').Where(l => l.Length == 12 && !l.StartsWith('M')).ToList();
        Assert.Equal(3, bondLines.Count(l => l.Substring(6, 3).Trim() == "2"));
        var radical = Molecule.FromSmiles("[CH3]"); // methyl radical: 3 H, not 4
        Assert.Equal("CH3", MolFile.Parse(MolFile.Write(radical)).Formula);
    }

    [Fact]
    public void Rejects_v3000()
    {
        const string v3 = "x\n\n\n  0  0  0     0  0            999 V3000\nM  END\n";
        Assert.Throws<FormatException>(() => MolFile.Parse(v3));
    }

    [Fact]
    public void Sdf_streaming_round_trip_with_properties()
    {
        var mols = new[] { "CCO", "c1ccccc1", "CC(=O)O" }.Select((s, i) =>
        {
            var m = Molecule.FromSmiles(s, $"mol{i}");
            m.Properties["ID"] = $"ID-{i}";
            m.Properties["Note"] = "line one\nline two";
            return m;
        }).ToList();
        string sdf = SdfWriter.ToSdf(mols);
        Assert.Equal(3, sdf.Split("$$$$").Length - 1);
        var read = SdfReader.Read(new StringReader(sdf)).ToList();
        Assert.Equal(3, read.Count);
        for (int i = 0; i < 3; i++)
        {
            Assert.Equal($"mol{i}", read[i].Name);
            Assert.Equal($"ID-{i}", read[i].Properties["ID"]);
            Assert.Equal("line one\nline two", read[i].Properties["Note"]);
            Assert.Equal(mols[i].ToSmiles(), read[i].ToSmiles());
        }
    }

    [Fact]
    public void Sdf_reader_is_lazy_and_can_skip_bad_records()
    {
        string good = SdfWriter.ToSdf([Molecule.FromSmiles("CCO", "a")]);
        string bad = "broken\n\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0 0.0 0.0 Zz\nM  END\n$$$$\n";
        string sdf = good + bad + good;
        Assert.Throws<FormatException>(() => SdfReader.Read(new StringReader(sdf)).ToList());
        var skipped = SdfReader.Read(new StringReader(sdf), new MolReaderOptions { SkipInvalid = true }).ToList();
        Assert.Equal(2, skipped.Count);
        // laziness: taking the first record never touches the broken one
        Assert.Equal("a", SdfReader.Read(new StringReader(sdf)).First().Name);
    }

    [Fact]
    public void Sdf_file_io_including_gzip()
    {
        string dir = Path.Combine(Path.GetTempPath(), "biotechnet-chem-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(dir);
        try
        {
            string path = Path.Combine(dir, "set.sdf");
            Assert.Equal(2, SdfWriter.WriteFile(path, [Molecule.FromSmiles("CCO", "e"), Molecule.FromSmiles("CCN", "a")]));
            Assert.Equal(["e", "a"], SdfReader.ReadFile(path).Select(m => m.Name).ToArray());
            string gz = path + ".gz";
            using (var fs = File.Create(gz))
            using (var z = new System.IO.Compression.GZipStream(fs, System.IO.Compression.CompressionLevel.Fastest))
            using (var src = File.OpenRead(path)) src.CopyTo(z);
            Assert.Equal(2, SdfReader.ReadFile(gz).Count());
        }
        finally { Directory.Delete(dir, true); }
    }
}

public class DepictionTests
{
    public static TheoryData<string> Molecules() =>
    [
        "CCCCCCCCCC", "CC(C)(C)C", "c1ccccc1", "c1ccc2ccccc2c1", "c1ccc2cc3ccccc3cc2c1", "C1CCC2(CC1)CCCC2",
        "CC(=O)Oc1ccccc1C(=O)O", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
        "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O",
        "CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C", "C#CC#CC", "OC(=O)CCC(=O)O.[Na+]",
        "C1CCCCCCCCCCC1", "c1ccc(-c2ccccc2)cc1", "CC(C)(C)c1cc(C(C)(C)C)cc(C(C)(C)C)c1",
    ];

    [Theory]
    [MemberData(nameof(Molecules))]
    public void Layout_has_uniform_bonds_and_no_clashes(string smiles)
    {
        var m = Molecule.FromSmiles(smiles);
        var c = Depiction.Compute2DCoordinates(m);
        Assert.Equal(m.AtomCount, c.Length);
        Assert.All(c, p => Assert.Equal(0, p.Z));
        foreach (var b in m.Bonds)
            Assert.InRange(Point3.Distance(c[b.Begin], c[b.End]), 1.3, 1.7);
        Assert.Equal(0, Depiction.CountClashes(m));
    }

    [Fact]
    public void Benzene_is_a_regular_hexagon()
    {
        var m = Molecule.FromSmiles("c1ccccc1");
        var c = Depiction.Compute2DCoordinates(m);
        double cx = c.Average(p => p.X), cy = c.Average(p => p.Y);
        Assert.All(c, p => Assert.Equal(1.5, Math.Sqrt((p.X - cx) * (p.X - cx) + (p.Y - cy) * (p.Y - cy)), 2));
    }

    [Fact]
    public void Fragments_do_not_overlap()
    {
        var m = Molecule.FromSmiles("c1ccccc1.c1ccccc1");
        var c = Depiction.Compute2DCoordinates(m);
        double maxA = c.Take(6).Max(p => p.X), minB = c.Skip(6).Min(p => p.X);
        Assert.True(minB > maxA + 1);
    }
}

public class ScaffoldTests
{
    [Theory]
    [InlineData("CC(=O)Oc1ccccc1C(=O)O", "c1ccccc1")]
    [InlineData("CCCC", "")]
    [InlineData("O=C1CCCCC1CC", "O=C1CCCCC1")]
    [InlineData("Cn1ccnc1", "c1c[nH]cn1")]
    [InlineData("c1ccccc1CCc1ccncc1C", "c1ccc(CCc2ccncc2)cc1")]
    public void Murcko_scaffolds(string smiles, string scaffold)
    {
        string expected = scaffold.Length == 0 ? "" : Molecule.FromSmiles(scaffold).ToSmiles();
        Assert.Equal(expected, MurckoScaffold.GetScaffoldSmiles(Molecule.FromSmiles(smiles)));
    }

    [Fact]
    public void Exocyclic_double_bonds_can_be_dropped_and_generic_framework()
    {
        var m = Molecule.FromSmiles("O=C1CCCCC1Cc1ccccc1");
        Assert.Equal(Molecule.FromSmiles("C1CCC(CC1)Cc1ccccc1").ToSmiles(), MurckoScaffold.GetScaffoldSmiles(m, keepExocyclicDoubleBonds: false));
        Assert.Equal(Molecule.FromSmiles("C1CCC(CC1)CC1CCCCC1").ToSmiles(), MurckoScaffold.GetGenericFrameworkSmiles(m));
    }
}
