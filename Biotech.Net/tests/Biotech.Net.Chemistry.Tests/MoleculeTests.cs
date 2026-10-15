using Biotech.Net.Chemistry;

namespace Biotech.Net.Chemistry.Tests;

public class SmilesAndGraphTests
{
    [Theory]
    [InlineData("c1ccccc1", 6, 6, 1)]
    [InlineData("c1ccc2ccccc2c1", 10, 11, 2)]
    [InlineData("C1CC2CCC1C2", 7, 8, 2)]
    [InlineData("C1C2CC3CC1CC(C2)C3", 10, 12, 3)]
    [InlineData("C1CCC2(CC1)CCCC2", 10, 11, 2)]
    [InlineData("CCO", 3, 2, 0)]
    public void Ring_perception_gives_sssr(string smiles, int atoms, int bonds, int rings)
    {
        var m = Molecule.FromSmiles(smiles);
        Assert.Equal(atoms, m.AtomCount);
        Assert.Equal(bonds, m.BondCount);
        Assert.Equal(rings, m.Rings.Count);
    }

    [Theory]
    [InlineData("c1ccccc1", true)]
    [InlineData("c1ccncc1", true)]
    [InlineData("c1cc[nH]c1", true)]
    [InlineData("c1ccoc1", true)]
    [InlineData("c1ccsc1", true)]
    [InlineData("C1=CC=CC=C1", true)] // Kekulé input is perceived as aromatic
    [InlineData("C1=CCC=C1", false)]  // cyclopentadiene
    [InlineData("C1=CC=CC=CC=C1", false)] // cyclooctatetraene (8 π)
    [InlineData("O=C1C=CC(=O)C=C1", false)] // benzoquinone
    [InlineData("O=c1cccc[nH]1", true)] // 2-pyridone (RDKit model)
    public void Aromaticity_of_classic_rings(string smiles, bool aromatic)
    {
        var m = Molecule.FromSmiles(smiles);
        int ringAtom = Enumerable.Range(0, m.AtomCount).First(i => m.Rings.IsAtomInRing(i));
        Assert.Equal(aromatic, m.Atoms[ringAtom].IsAromatic);
    }

    [Fact]
    public void Kekule_and_aromatic_inputs_give_the_same_canonical_smiles()
    {
        Assert.Equal(Molecule.FromSmiles("c1ccccc1O").ToSmiles(), Molecule.FromSmiles("OC1=CC=CC=C1").ToSmiles());
        Assert.Equal(Molecule.FromSmiles("Cn1cnc2c1c(=O)n(C)c(=O)n2C").ToSmiles(), Molecule.FromSmiles("CN1C=NC2=C1C(=O)N(C)C(=O)N2C").ToSmiles());
    }

    [Fact]
    public void Kekulize_produces_alternating_bonds()
    {
        var k = Aromaticity.Kekulize(Molecule.FromSmiles("c1ccccc1"));
        Assert.Equal(3, k.Bonds.ToArray().Count(b => b.Order == BondOrder.Double));
        Assert.DoesNotContain(k.Atoms.ToArray(), a => a.IsAromatic);
    }

    [Theory]
    [InlineData("C[C@H](N)C(=O)O")]
    [InlineData("C[C@@H](N)C(=O)O")]
    [InlineData("F/C=C/F")]
    [InlineData("F/C=C\\F")]
    [InlineData("[NH4+].[Cl-]")]
    [InlineData("[13CH4]")]
    [InlineData("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O")]
    public void Canonical_smiles_round_trips(string smiles)
    {
        var a = Molecule.FromSmiles(smiles).ToSmiles();
        Assert.Equal(a, Molecule.FromSmiles(a).ToSmiles());
    }

    [Fact]
    public void Stereo_isomers_have_different_canonical_smiles()
    {
        Assert.NotEqual(Molecule.FromSmiles("C[C@H](N)C(=O)O").ToSmiles(), Molecule.FromSmiles("C[C@@H](N)C(=O)O").ToSmiles());
        Assert.Equal(Molecule.FromSmiles("C[C@H](N)C(=O)O").ToSmiles(), Molecule.FromSmiles("N[C@@H](C)C(=O)O").ToSmiles());
        Assert.NotEqual(Molecule.FromSmiles("F/C=C/F").ToSmiles(), Molecule.FromSmiles("F/C=C\\F").ToSmiles());
        Assert.Equal(Molecule.FromSmiles("F/C=C/F").ToSmiles(), Molecule.FromSmiles("F\\C=C\\F").ToSmiles());
    }

    [Theory]
    [InlineData("C1CC")]
    [InlineData("C(C")]
    [InlineData("C)")]
    [InlineData("[Xx]")]
    public void Malformed_smiles_throw(string smiles)
    {
        Assert.Throws<FormatException>(() => Molecule.FromSmiles(smiles));
        Assert.Null(SmilesParser.TryParse(smiles));
    }

    [Fact]
    public void Hydrogens_can_be_added_and_removed()
    {
        var m = Molecule.FromSmiles("CCO");
        var h = m.AddHydrogens();
        Assert.Equal(9, h.AtomCount);
        Assert.Equal("C2H6O", h.Formula);
        var back = h.RemoveHydrogens();
        Assert.Equal(3, back.AtomCount);
        Assert.Equal(m.ToSmiles(), back.ToSmiles());
    }

    [Fact]
    public void Added_hydrogens_get_coordinates_when_the_molecule_has_them()
    {
        var m = Molecule.FromSmiles("CC(=O)O");
        Depiction.Compute2DCoordinates(m);
        var h = m.AddHydrogens();
        Assert.NotNull(h.Coordinates);
        for (int i = m.AtomCount; i < h.AtomCount; i++)
        {
            int parent = h.Neighbors(i)[0];
            double d = Point3.Distance(h.Coordinates![i], h.Coordinates[parent]);
            Assert.InRange(d, 0.9, 1.2);
        }
    }

    [Fact]
    public void Formula_and_masses()
    {
        var aspirin = Molecule.FromSmiles("CC(=O)Oc1ccccc1C(=O)O");
        Assert.Equal("C9H8O4", aspirin.Formula);
        Assert.Equal(180.159, MolecularDescriptors.MolecularWeight(aspirin), 3);
        Assert.Equal(180.042259, MolecularDescriptors.ExactMass(aspirin), 5);
        Assert.Equal("C2H8N+", Molecule.FromSmiles("CC[NH3+]").Formula);
        Assert.Equal("CHD3O", MolecularDescriptors.Formula(Molecule.FromSmiles("[2H]C([2H])([2H])O"), separateIsotopes: true));
    }

    [Theory]
    [InlineData("C9H8O4", 180.159, 180.042259)]
    [InlineData("H2O", 18.015, 18.010565)]
    [InlineData("Ca(OH)2", 74.092, 73.968070)]
    [InlineData("CuSO4·5H2O", 249.677, 248.934150)]
    public void Molecular_formula_parsing(string formula, double average, double mono)
    {
        var f = MolecularFormula.Parse(formula);
        Assert.Equal(average, f.AverageMass, 2);
        Assert.Equal(mono, f.MonoisotopicMass, 4);
    }

    [Fact]
    public void Molecular_formula_charge_and_mz()
    {
        var f = MolecularFormula.Parse("C2H8N+");
        Assert.Equal(1, f.Charge);
        Assert.Equal("C2H8N+", f.ToString());
        Assert.Equal(f.MonoisotopicMass - 0.00054858, f.MassToCharge, 6);
        Assert.Equal(-2, MolecularFormula.Parse("SO4-2").Charge);
        Assert.Throws<FormatException>(() => MolecularFormula.Parse("C9Xx"));
    }

    [Fact]
    public void Fragments_and_distances()
    {
        var m = Molecule.FromSmiles("CCO.[Na+]");
        Assert.Equal(2, m.Fragments().Count);
        var d = m.TopologicalDistances();
        Assert.Equal(2, d[0, 2]);
        Assert.Equal(-1, d[0, 3]);
    }
}

public class SmartsTests
{
    [Theory]
    [InlineData("CC(=O)Oc1ccccc1C(=O)O", "C(=O)[OH]", 1)]
    [InlineData("CC(=O)Oc1ccccc1C(=O)O", "[CX3](=O)[OX2H0]", 1)] // ester
    [InlineData("CC(=O)Oc1ccccc1C(=O)O", "c1ccccc1", 1)]
    [InlineData("CC(=O)Oc1ccccc1C(=O)O", "[#6]", 9)]
    [InlineData("c1ccncc1", "[nX2]", 1)]
    [InlineData("c1ccncc1", "[$(n1ccccc1)]", 1)]
    [InlineData("C1CCCCC1", "[R]", 6)]
    [InlineData("C1CCCCC1", "[r5]", 0)]
    [InlineData("CCN(CC)CC", "[NX3;H0;!$(NC=O)]", 1)]
    [InlineData("CC(=O)NC", "[NX3;H0;!$(NC=O)]", 0)]
    [InlineData("OCC(O)CO", "[OH]", 3)]
    [InlineData("C[N+](=O)[O-]", "[N+](=O)[O-]", 1)]
    [InlineData("c1ccccc1Cl", "c[Cl,Br,I]", 1)]
    [InlineData("CC=CC", "C=,#C", 1)]
    [InlineData("CCCC", "C!@C", 3)]
    [InlineData("C1CCC1", "C@C", 4)]
    public void Pattern_match_counts(string smiles, string smarts, int expected)
        => Assert.Equal(expected, SmartsPattern.Parse(smarts).CountMatches(Molecule.FromSmiles(smiles)));

    [Fact]
    public void Explicit_hydrogen_query_atoms_are_merged()
    {
        var p = SmartsPattern.Parse("[#6]([#1])([#1])[#1]");
        Assert.Equal(1, p.AtomCount); // H atoms folded into "at least 3 H"
        Assert.Equal(1, p.CountMatches(Molecule.FromSmiles("CCO"))); // the methyl, not the CH2
    }

    [Fact]
    public void Atom_maps_are_recorded()
    {
        var p = SmartsPattern.Parse("[C:1](=O)[O:2]");
        Assert.Equal(1, p.AtomMap(0));
        Assert.Equal(0, p.AtomMap(1));
        Assert.Equal(2, p.AtomMap(2));
    }

    [Fact]
    public void Molecule_as_query()
    {
        Assert.True(Substructure.HasMatch(Molecule.FromSmiles("CC(=O)Oc1ccccc1C(=O)O"), Molecule.FromSmiles("c1ccccc1")));
        Assert.False(Substructure.HasMatch(Molecule.FromSmiles("C1CCCCC1"), Molecule.FromSmiles("c1ccccc1")));
    }

    [Fact]
    public void Malformed_smarts_throws()
    {
        Assert.Throws<FormatException>(() => SmartsPattern.Parse("[C"));
        Assert.Null(SmartsPattern.TryParse("C(("));
    }
}
