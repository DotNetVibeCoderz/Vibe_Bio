namespace Biotech.Net.StructuralBiology.Tests;

internal static class TestData
{
    public static string PathOf(string name) => Path.Combine(AppContext.BaseDirectory, "Data", name);

    private static Structure? _trypsin;

    /// <summary>
    /// Bovine trypsin with benzamidine (PDB 3PTB), coordinate records only. The instance is parsed once and shared by
    /// all tests, so callers must treat it as read-only (no <see cref="Model.Transform"/> or <see cref="Atom.Position"/> writes).
    /// </summary>
    public static Structure Trypsin() => _trypsin ??= PdbReader.Read(PathOf("3ptb.pdb"));

    /// <summary>
    /// Reference DSSP string for 1CRN from DSSP 4.6.1 (PDB-REDO DSSP databank, pdb-redo.eu/dssp/db/1crn/legacy;
    /// Hekkelman et al. 2025, Kabsch &amp; Sander 1983), column 17 of the legacy output, loops as '-'.
    /// </summary>
    public static string CrambinDsspReference()
    {
        var lines = File.ReadAllLines(PathOf("1crn.dssp")).Skip(1).Where(l => l.Length > 16);
        return string.Concat(lines.Select(l => l[16] == ' ' ? '-' : l[16]));
    }
}
