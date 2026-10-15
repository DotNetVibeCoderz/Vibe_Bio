using Biotech.Net.Core;

namespace Biotech.Net.StructuralBiology;

/// <summary>
/// One atom of a macromolecular structure. The <see cref="Position"/> is mutable so superpositions and docking
/// poses can be applied in place; everything else is fixed after parsing.
/// </summary>
public sealed class Atom
{
    private List<Atom>? _alternates;

    /// <summary>Creates an atom.</summary>
    /// <param name="name">Atom name, e.g. <c>CA</c> (trimmed).</param>
    /// <param name="element">Element symbol, e.g. <c>C</c> or <c>FE</c> (normalised to upper case).</param>
    /// <param name="position">Cartesian coordinates in Å.</param>
    public Atom(string name, string element, Vec3 position)
    {
        Name = name.Trim();
        Element = ElementData.Normalize(element);
        Position = position;
    }

    /// <summary>Atom name (e.g. <c>CA</c>, <c>OG1</c>).</summary>
    public string Name { get; }

    /// <summary>Upper-case element symbol (e.g. <c>C</c>, <c>SE</c>, <c>ZN</c>).</summary>
    public string Element { get; }

    /// <summary>Coordinates in Å.</summary>
    public Vec3 Position { get; set; }

    /// <summary>Atom serial number from the file.</summary>
    public int Serial { get; init; }

    /// <summary>Alternate location indicator (<c>' '</c> when none).</summary>
    public char AltLoc { get; init; } = ' ';

    /// <summary>Occupancy (0–1).</summary>
    public double Occupancy { get; init; } = 1.0;

    /// <summary>Isotropic temperature factor (B-factor, Å²).</summary>
    public double BFactor { get; init; }

    /// <summary>Formal charge (0 when not given).</summary>
    public int Charge { get; init; }

    /// <summary>True for HETATM records (ligands, ions, water, modified residues).</summary>
    public bool IsHetero { get; init; }

    /// <summary>The residue this atom belongs to (null for a free-standing atom).</summary>
    public Residue? Residue { get; internal set; }

    /// <summary>Other alternate-location conformers of this atom (this atom is the selected one, highest occupancy).</summary>
    public IReadOnlyList<Atom> Alternates => (IReadOnlyList<Atom>?)_alternates ?? [];

    /// <summary>Residue name, or empty.</summary>
    public string ResidueName => Residue?.Name ?? string.Empty;

    /// <summary>Residue sequence number, or 0.</summary>
    public int ResidueNumber => Residue?.Number ?? 0;

    /// <summary>Chain identifier, or empty.</summary>
    public string ChainId => Residue?.Chain?.Id ?? string.Empty;

    /// <summary>Average atomic mass of <see cref="Element"/>.</summary>
    public double Mass => ElementData.Mass(Element);

    /// <summary>Van der Waals radius of <see cref="Element"/>.</summary>
    public double VdwRadius => ElementData.VdwRadius(Element);

    /// <summary>True for protein backbone atoms (N, CA, C, O, OXT) of amino-acid residues.</summary>
    public bool IsBackbone => Residue is { IsAminoAcid: true } && Name is "N" or "CA" or "C" or "O" or "OXT";

    internal void AddAlternate(Atom other) => (_alternates ??= []).Add(other);

    /// <summary>A copy at the same position, detached from any residue.</summary>
    public Atom Clone() => new(Name, Element, Position)
    {
        Serial = Serial, AltLoc = AltLoc, Occupancy = Occupancy, BFactor = BFactor, Charge = Charge, IsHetero = IsHetero,
    };

    public override string ToString() => Residue is null ? $"{Name} {Position}" : $"{ChainId}:{ResidueName}{ResidueNumber}:{Name}";
}

/// <summary>A residue (amino acid, nucleotide, ligand, ion or water molecule).</summary>
public sealed class Residue
{
    private readonly List<Atom> _atoms = [];

    /// <summary>Creates an empty residue.</summary>
    public Residue(string name, int number, char insertionCode = ' ', bool isHetero = false)
    {
        Name = name.Trim();
        Number = number;
        InsertionCode = insertionCode;
        IsHetero = isHetero;
    }

    /// <summary>Three-letter residue name (<c>ALA</c>, <c>HOH</c>, <c>BEN</c>, ...).</summary>
    public string Name { get; }

    /// <summary>Author residue sequence number.</summary>
    public int Number { get; }

    /// <summary>Insertion code (<c>' '</c> when none).</summary>
    public char InsertionCode { get; }

    /// <summary>True when the residue came from HETATM records.</summary>
    public bool IsHetero { get; }

    /// <summary>The owning chain.</summary>
    public Chain? Chain { get; internal set; }

    /// <summary>Atoms in file order (selected alternate locations only).</summary>
    public IReadOnlyList<Atom> Atoms => _atoms;

    /// <summary>True for amino acids (standard or common modified names).</summary>
    public bool IsAminoAcid => ResidueNames.IsAminoAcid(Name);

    /// <summary>True for water.</summary>
    public bool IsWater => ResidueNames.IsWater(Name);

    /// <summary>True for HETATM residues that are not water and not a modified amino acid: ligands, cofactors, ions.</summary>
    public bool IsLigand => IsHetero && !IsWater && !IsAminoAcid;

    /// <summary>One-letter amino-acid code (<c>'X'</c> for non-amino acids).</summary>
    public char OneLetterCode => ResidueNames.ToOneLetter(Name);

    /// <summary>The atom called <paramref name="atomName"/>, or null.</summary>
    public Atom? GetAtom(string atomName)
    {
        foreach (var a in _atoms) if (a.Name == atomName) return a;
        return null;
    }

    /// <summary>The atom called <paramref name="atomName"/>.</summary>
    /// <exception cref="KeyNotFoundException">No such atom.</exception>
    public Atom this[string atomName] => GetAtom(atomName) ?? throw new KeyNotFoundException($"{this} has no atom '{atomName}'.");

    /// <summary>Whether the residue has an atom named <paramref name="atomName"/>.</summary>
    public bool HasAtom(string atomName) => GetAtom(atomName) is not null;

    /// <summary>Adds an atom. An atom whose name already exists with a different altloc becomes an alternate
    /// conformer; the conformer with the highest occupancy (first on ties) stays the selected one.</summary>
    public void AddAtom(Atom atom)
    {
        if (atom.AltLoc != ' ')
        {
            for (int i = 0; i < _atoms.Count; i++)
            {
                var existing = _atoms[i];
                if (existing.Name != atom.Name || existing.AltLoc == ' ') continue;
                if (atom.Occupancy > existing.Occupancy)
                {
                    foreach (var alt in existing.Alternates) atom.AddAlternate(alt);
                    atom.AddAlternate(existing);
                    existing.Residue = null;
                    _atoms[i] = atom;
                    atom.Residue = this;
                }
                else
                {
                    existing.AddAlternate(atom);
                    atom.Residue = this;
                }
                return;
            }
        }
        atom.Residue = this;
        _atoms.Add(atom);
    }

    /// <summary>Identifier like <c>A:CYS3</c> (with insertion code when present).</summary>
    public string Id => $"{Chain?.Id}:{Name}{Number}{(InsertionCode == ' ' ? "" : InsertionCode.ToString())}";

    public override string ToString() => Id;
}

/// <summary>A polymer chain (plus the hetero groups and waters assigned to it).</summary>
public sealed class Chain
{
    private readonly List<Residue> _residues = [];

    /// <summary>Creates an empty chain.</summary>
    public Chain(string id) => Id = id;

    /// <summary>Chain identifier (author chain id).</summary>
    public string Id { get; }

    /// <summary>The owning model.</summary>
    public Model? Model { get; internal set; }

    /// <summary>Residues in file order.</summary>
    public IReadOnlyList<Residue> Residues => _residues;

    /// <summary>Adds a residue.</summary>
    public void AddResidue(Residue residue)
    {
        residue.Chain = this;
        _residues.Add(residue);
    }

    /// <summary>Amino-acid residues that have a CA atom, in order (the protein backbone trace).</summary>
    public IReadOnlyList<Residue> AminoAcids() => _residues.Where(r => r.IsAminoAcid && r.HasAtom("CA")).ToList();

    /// <summary>The residue with author number <paramref name="number"/> (and insertion code), or null.</summary>
    public Residue? FindResidue(int number, char insertionCode = ' ')
    {
        foreach (var r in _residues) if (r.Number == number && r.InsertionCode == insertionCode) return r;
        return null;
    }

    /// <summary>All atoms of the chain.</summary>
    public IEnumerable<Atom> Atoms => _residues.SelectMany(r => r.Atoms);

    /// <summary>The amino-acid sequence of the chain as a Core protein <see cref="Sequence"/>.</summary>
    public Sequence ToSequence(string? id = null)
    {
        var chars = AminoAcids().Select(r => r.OneLetterCode).ToArray();
        return Sequence.Protein(new string(chars), id ?? $"{Model?.Structure?.Id}_{Id}".TrimStart('_'));
    }

    public override string ToString() => $"Chain {Id} ({_residues.Count} residues)";
}

/// <summary>One model (conformation) of a structure — NMR ensembles have several.</summary>
public sealed class Model
{
    private readonly List<Chain> _chains = [];

    /// <summary>Creates an empty model.</summary>
    public Model(int number = 1) => Number = number;

    /// <summary>Model serial number (1 for single-model files).</summary>
    public int Number { get; }

    /// <summary>The owning structure.</summary>
    public Structure? Structure { get; internal set; }

    /// <summary>Chains in file order.</summary>
    public IReadOnlyList<Chain> Chains => _chains;

    /// <summary>The chain with identifier <paramref name="chainId"/>.</summary>
    public Chain this[string chainId] => GetChain(chainId) ?? throw new KeyNotFoundException($"No chain '{chainId}'.");

    /// <summary>The chain with identifier <paramref name="chainId"/>, or null.</summary>
    public Chain? GetChain(string chainId) => _chains.FirstOrDefault(c => c.Id == chainId);

    /// <summary>Adds a chain.</summary>
    public void AddChain(Chain chain)
    {
        chain.Model = this;
        _chains.Add(chain);
    }

    /// <summary>Every residue of every chain.</summary>
    public IEnumerable<Residue> Residues => _chains.SelectMany(c => c.Residues);

    /// <summary>Every atom (selected altlocs) in file order.</summary>
    public IEnumerable<Atom> Atoms => _chains.SelectMany(c => c.Atoms);

    /// <summary>The atoms as a list (materialised once per call).</summary>
    public List<Atom> GetAtoms() => Atoms.ToList();

    /// <summary>The coordinates of <see cref="Atoms"/> as a contiguous array.</summary>
    public Vec3[] GetCoordinates() => Atoms.Select(a => a.Position).ToArray();

    /// <summary>HETATM residues that are neither water nor amino acids (ligands, cofactors, ions).</summary>
    public IEnumerable<Residue> Ligands => Residues.Where(r => r.IsLigand);

    /// <summary>Applies <c>x → R·x + t</c> to every atom (including alternates).</summary>
    public void Transform(Matrix3 rotation, Vec3 translation)
    {
        foreach (var a in Atoms)
        {
            a.Position = rotation * a.Position + translation;
            foreach (var alt in a.Alternates) alt.Position = rotation * alt.Position + translation;
        }
    }

    public override string ToString() => $"Model {Number} ({_chains.Count} chains)";
}

/// <summary>Unit cell and space group from a CRYST1 record or <c>_cell</c>/<c>_symmetry</c> categories.</summary>
public sealed record UnitCell(double A, double B, double C, double Alpha, double Beta, double Gamma, string SpaceGroup, int Z)
{
    /// <summary>Cell volume in Å³.</summary>
    public double Volume
    {
        get
        {
            double ca = Math.Cos(Alpha * Math.PI / 180), cb = Math.Cos(Beta * Math.PI / 180), cg = Math.Cos(Gamma * Math.PI / 180);
            return A * B * C * Math.Sqrt(1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg);
        }
    }
}

/// <summary>A HELIX record: a helix from the start to the end residue (inclusive).</summary>
/// <param name="HelixClass">PDB helix class: 1 right-handed α, 5 right-handed 3₁₀, 3 π, ...</param>
public sealed record HelixRecord(int Serial, string HelixId, string StartChain, string StartResidueName, int StartNumber,
    string EndChain, string EndResidueName, int EndNumber, int HelixClass, int Length);

/// <summary>A SHEET record: one strand of a β-sheet.</summary>
/// <param name="Sense">0 for the first strand, 1 parallel, −1 antiparallel to the previous strand.</param>
public sealed record SheetRecord(int Strand, string SheetId, int StrandCount, string StartChain, string StartResidueName, int StartNumber,
    string EndChain, string EndResidueName, int EndNumber, int Sense);

/// <summary>An SSBOND record: a disulfide bridge between two cysteines.</summary>
public sealed record SsBondRecord(int Serial, string ChainA, int NumberA, string ChainB, int NumberB, double Length);

/// <summary>
/// A macromolecular structure: models → chains → residues → atoms, plus the header records that were read
/// (unit cell, HELIX/SHEET/SSBOND annotation and CONECT bonds).
/// </summary>
public sealed class Structure
{
    private readonly List<Model> _models = [];

    /// <summary>Creates an empty structure.</summary>
    public Structure(string id) => Id = id;

    /// <summary>Identifier (PDB id or file name).</summary>
    public string Id { get; }

    /// <summary>Title from TITLE / <c>_struct.title</c>.</summary>
    public string Title { get; set; } = string.Empty;

    /// <summary>Header classification from HEADER.</summary>
    public string Classification { get; set; } = string.Empty;

    /// <summary>Unit cell, when present.</summary>
    public UnitCell? Cell { get; set; }

    /// <summary>Models in file order.</summary>
    public IReadOnlyList<Model> Models => _models;

    /// <summary>The first model.</summary>
    public Model FirstModel => _models.Count > 0 ? _models[0] : throw new InvalidOperationException($"Structure {Id} has no models.");

    /// <summary>HELIX records.</summary>
    public List<HelixRecord> Helices { get; } = [];

    /// <summary>SHEET records.</summary>
    public List<SheetRecord> Sheets { get; } = [];

    /// <summary>SSBOND records.</summary>
    public List<SsBondRecord> Disulfides { get; } = [];

    /// <summary>CONECT bonds: atom serial → bonded atom serials.</summary>
    public Dictionary<int, List<int>> Connections { get; } = [];

    /// <summary>Chains of the first model.</summary>
    public IReadOnlyList<Chain> Chains => FirstModel.Chains;

    /// <summary>Atoms of the first model.</summary>
    public IEnumerable<Atom> Atoms => FirstModel.Atoms;

    /// <summary>Residues of the first model.</summary>
    public IEnumerable<Residue> Residues => FirstModel.Residues;

    /// <summary>Ligands (non-water, non-amino-acid HETATM residues) of the first model.</summary>
    public IEnumerable<Residue> Ligands => FirstModel.Ligands;

    /// <summary>The first-model chain called <paramref name="chainId"/>.</summary>
    public Chain this[string chainId] => FirstModel[chainId];

    /// <summary>Adds a model.</summary>
    public void AddModel(Model model)
    {
        model.Structure = this;
        _models.Add(model);
    }

    /// <summary>Records a CONECT bond in both directions (duplicates ignored).</summary>
    public void AddConnection(int serialA, int serialB)
    {
        Add(serialA, serialB);
        Add(serialB, serialA);

        void Add(int a, int b)
        {
            if (!Connections.TryGetValue(a, out var list)) Connections[a] = list = [];
            if (!list.Contains(b)) list.Add(b);
        }
    }

    /// <summary>Number of atoms in the first model.</summary>
    public int AtomCount => _models.Count == 0 ? 0 : Atoms.Count();

    public override string ToString() =>
        $"{Id}: {_models.Count} model(s), {(_models.Count > 0 ? Chains.Count : 0)} chain(s), {AtomCount} atoms";
}
