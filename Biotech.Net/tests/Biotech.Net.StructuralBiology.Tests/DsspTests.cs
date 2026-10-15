namespace Biotech.Net.StructuralBiology.Tests;

public class DsspTests
{
    [Fact]
    public void Crambin_MatchesReferenceDssp()
    {
        // Reference: DSSP 4.6.1 output for 1CRN (PDB-REDO databank), stored in Data/1crn.dssp.
        var result = Dssp.Assign(SampleStructures.Crambin());
        Assert.Equal(TestData.CrambinDsspReference(), result.ToSecondaryStructureString());
    }
}
