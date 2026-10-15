using System.Reflection;
using System.Runtime.InteropServices;

namespace Biotech.Net.Native;

/// <summary>
/// P/Invoke bridge to the Rust library <c>biotech_native</c> (native/biotech-native). Every method mirrors a
/// <c>bn_*</c> export one to one; argument validation happens in <see cref="BioKernels"/> before a call gets here.
/// </summary>
/// <remarks>
/// The library is resolved through <see cref="NativeLibrary"/>: the default probing paths first, then
/// <c>runtimes/{rid}/native</c> and the assembly directory. <see cref="IsAvailable"/> is false (never an exception)
/// when the library is missing, fails to load, reports a different <see cref="ExpectedAbiVersion"/>, or when the
/// environment variable <c>BIOTECH_NET_DISABLE_NATIVE=1</c> is set.
/// </remarks>
internal static unsafe partial class NativeBridge
{
    internal const string LibraryName = "biotech_native";

    /// <summary>The C ABI version this bridge was written against (<c>ABI_VERSION</c> in lib.rs).</summary>
    internal const uint ExpectedAbiVersion = 1;

    /// <summary>Environment variable that forces the managed kernels when set to <c>1</c> or <c>true</c>.</summary>
    internal const string DisableEnvironmentVariable = "BIOTECH_NET_DISABLE_NATIVE";

    private static readonly Lazy<(bool Available, string Version, string? Reason)> State = new(Probe);

    /// <summary>True when the library was loaded and its ABI version matches.</summary>
    public static bool IsAvailable => State.Value.Available;

    /// <summary>Version string reported by the library, or <c>"managed"</c>.</summary>
    public static string Version => State.Value.Version;

    /// <summary>Why the native library is not used (null when it is).</summary>
    public static string? UnavailableReason => State.Value.Reason;

    static NativeBridge()
    {
        try { NativeLibrary.SetDllImportResolver(typeof(NativeBridge).Assembly, Resolve); }
        catch (InvalidOperationException) { /* a resolver was already registered for this assembly */ }
    }

    private static (bool, string, string?) Probe()
    {
        string? flag = Environment.GetEnvironmentVariable(DisableEnvironmentVariable);
        if (flag is not null && (flag.Trim() == "1" || flag.Trim().Equals("true", StringComparison.OrdinalIgnoreCase)))
            return (false, "managed", $"{DisableEnvironmentVariable} is set");
        try
        {
            uint abi = bn_abi_version();
            if (abi != ExpectedAbiVersion)
                return (false, "managed", $"ABI version mismatch: library {abi}, bridge {ExpectedAbiVersion}");
            string version = Marshal.PtrToStringUTF8((nint)bn_version()) ?? "biotech-native";
            return (true, version, null);
        }
        catch (Exception ex) when (ex is DllNotFoundException or EntryPointNotFoundException or BadImageFormatException or TypeInitializationException)
        {
            return (false, "managed", ex.Message);
        }
    }

    private static nint Resolve(string name, Assembly assembly, DllImportSearchPath? searchPath)
    {
        if (name != LibraryName) return 0;
        if (NativeLibrary.TryLoad(name, assembly, searchPath, out nint handle)) return handle;

        string file = OperatingSystem.IsWindows() ? "biotech_native.dll"
            : OperatingSystem.IsMacOS() ? "libbiotech_native.dylib"
            : "libbiotech_native.so";
        foreach (string dir in CandidateDirectories(assembly))
        {
            string path = Path.Combine(dir, file);
            if (File.Exists(path) && NativeLibrary.TryLoad(path, out handle)) return handle;
        }
        return 0; // fall through to the runtime's default (which throws DllNotFoundException)
    }

    private static IEnumerable<string> CandidateDirectories(Assembly assembly)
    {
        string arch = RuntimeInformation.ProcessArchitecture.ToString().ToLowerInvariant();
        string os = OperatingSystem.IsWindows() ? "win" : OperatingSystem.IsMacOS() ? "osx" : "linux";
        string[] rids = [RuntimeInformation.RuntimeIdentifier, $"{os}-{arch}"];
        var roots = new List<string>();
        if (!string.IsNullOrEmpty(assembly.Location)) roots.Add(Path.GetDirectoryName(assembly.Location)!);
        roots.Add(AppContext.BaseDirectory);
        foreach (string root in roots.Distinct())
        {
            foreach (string rid in rids.Distinct()) yield return Path.Combine(root, "runtimes", rid, "native");
            yield return root;
        }
    }

    // ---- Raw imports (signatures match lib.rs exactly; usize ↔ nuint) -------------------------------

    [LibraryImport(LibraryName)] private static partial uint bn_abi_version();
    [LibraryImport(LibraryName)] private static partial byte* bn_version();
    [LibraryImport(LibraryName)] private static partial void bn_reverse_complement(byte* src, nuint len, byte* dst);
    [LibraryImport(LibraryName)] private static partial void bn_complement(byte* src, nuint len, byte* dst);
    [LibraryImport(LibraryName)] private static partial nuint bn_translate(byte* nuc, nuint len, byte* table, byte* protein);
    [LibraryImport(LibraryName)] private static partial void bn_count_nucleotides(byte* s, nuint len, long* output);
    [LibraryImport(LibraryName)] private static partial void bn_byte_histogram(byte* data, nuint len, long* counts);
    [LibraryImport(LibraryName)] private static partial long bn_count_kmers(byte* s, nuint len, uint k, int* counts);
    [LibraryImport(LibraryName)] private static partial nuint bn_hamming(byte* a, byte* b, nuint len);

    [LibraryImport(LibraryName)]
    private static partial void bn_welch_ttest(double* matrix, nuint rows, nuint cols, int* ga, nuint na, int* gb, nuint nb,
        double* t, double* p, double* diff);

    [LibraryImport(LibraryName)] private static partial void bn_benjamini_hochberg(double* p, nuint n, double* q);
    [LibraryImport(LibraryName)] private static partial double bn_protein_mass(byte* s, nuint len, int monoisotopic);
    [LibraryImport(LibraryName)] private static partial double bn_net_charge(byte* s, nuint len, double ph);
    [LibraryImport(LibraryName)] private static partial double bn_isoelectric_point(byte* s, nuint len);
    [LibraryImport(LibraryName)] private static partial nuint bn_hydropathy(byte* s, nuint len, nuint window, double* output);

    [LibraryImport(LibraryName)]
    private static partial int bn_align(byte* a, nuint alen, byte* b, nuint blen, int* matrix, int gapOpen, int gapExtend,
        int mode, long band, byte* ops, nuint opsCap, int* result);

    [LibraryImport(LibraryName)]
    private static partial int bn_local_score(byte* q, nuint qlen, byte* t, nuint tlen, int* matrix, int gapOpen, int gapExtend);

    [LibraryImport(LibraryName)]
    private static partial int bn_local_score_batch(byte* q, nuint qlen, byte* targets, nuint targetsLen, int* offsets, nuint count,
        int* matrix, int gapOpen, int gapExtend, int* scores);

    [LibraryImport(LibraryName)] private static partial double bn_kabsch_rmsd(double* a, double* b, nuint n, double* rotation);

    // ---- Span wrappers ------------------------------------------------------------------------------

    public static void ReverseComplement(ReadOnlySpan<byte> s, Span<byte> d)
    {
        fixed (byte* ps = s) fixed (byte* pd = d) bn_reverse_complement(ps, (nuint)s.Length, pd);
    }

    public static void Complement(ReadOnlySpan<byte> s, Span<byte> d)
    {
        fixed (byte* ps = s) fixed (byte* pd = d) bn_complement(ps, (nuint)s.Length, pd);
    }

    public static int Translate(ReadOnlySpan<byte> n, ReadOnlySpan<byte> t, Span<byte> p)
    {
        fixed (byte* pn = n) fixed (byte* pt = t) fixed (byte* pp = p) return (int)bn_translate(pn, (nuint)n.Length, pt, pp);
    }

    public static NucleotideCounts CountNucleotides(ReadOnlySpan<byte> s)
    {
        long* c = stackalloc long[6];
        fixed (byte* ps = s) bn_count_nucleotides(ps, (nuint)s.Length, c);
        return new NucleotideCounts(c[0], c[1], c[2], c[3], c[4], c[5]);
    }

    public static void ByteHistogram(ReadOnlySpan<byte> s, Span<long> c)
    {
        fixed (byte* ps = s) fixed (long* pc = c) bn_byte_histogram(ps, (nuint)s.Length, pc);
    }

    public static long CountKmers(ReadOnlySpan<byte> s, int k, Span<int> c)
    {
        fixed (byte* ps = s) fixed (int* pc = c) return bn_count_kmers(ps, (nuint)s.Length, (uint)k, pc);
    }

    public static int Hamming(ReadOnlySpan<byte> a, ReadOnlySpan<byte> b)
    {
        fixed (byte* pa = a) fixed (byte* pb = b) return (int)bn_hamming(pa, pb, (nuint)a.Length);
    }

    public static void WelchTTest(ReadOnlySpan<double> m, int r, int c, ReadOnlySpan<int> a, ReadOnlySpan<int> b,
        Span<double> t, Span<double> p, Span<double> d)
    {
        fixed (double* pm = m) fixed (int* pa = a) fixed (int* pb = b)
        fixed (double* pt = t) fixed (double* pp = p) fixed (double* pd = d)
            bn_welch_ttest(pm, (nuint)r, (nuint)c, pa, (nuint)a.Length, pb, (nuint)b.Length, pt, pp, pd);
    }

    public static void BenjaminiHochberg(ReadOnlySpan<double> p, Span<double> q)
    {
        fixed (double* pp = p) fixed (double* pq = q) bn_benjamini_hochberg(pp, (nuint)p.Length, pq);
    }

    public static double ProteinMass(ReadOnlySpan<byte> s, bool mono)
    {
        fixed (byte* ps = s) return bn_protein_mass(ps, (nuint)s.Length, mono ? 1 : 0);
    }

    public static double NetCharge(ReadOnlySpan<byte> s, double pH)
    {
        fixed (byte* ps = s) return bn_net_charge(ps, (nuint)s.Length, pH);
    }

    public static double IsoelectricPoint(ReadOnlySpan<byte> s)
    {
        fixed (byte* ps = s) return bn_isoelectric_point(ps, (nuint)s.Length);
    }

    public static int Hydropathy(ReadOnlySpan<byte> s, int w, Span<double> o)
    {
        fixed (byte* ps = s) fixed (double* po = o) return (int)bn_hydropathy(ps, (nuint)s.Length, (nuint)w, po);
    }

    /// <summary>Calls <c>bn_align</c>; returns the status code (0 or negative) and fills <paramref name="result"/> (6 ints).</summary>
    public static int Align(ReadOnlySpan<byte> a, ReadOnlySpan<byte> b, ReadOnlySpan<int> matrix, int gapOpen, int gapExtend,
        int mode, long band, Span<byte> ops, Span<int> result)
    {
        fixed (byte* pa = a) fixed (byte* pb = b) fixed (int* pm = matrix) fixed (byte* po = ops) fixed (int* pr = result)
            return bn_align(pa, (nuint)a.Length, pb, (nuint)b.Length, pm, gapOpen, gapExtend, mode, band, po, (nuint)ops.Length, pr);
    }

    public static int LocalScore(ReadOnlySpan<byte> q, ReadOnlySpan<byte> t, ReadOnlySpan<int> matrix, int gapOpen, int gapExtend)
    {
        fixed (byte* pq = q) fixed (byte* pt = t) fixed (int* pm = matrix)
            return bn_local_score(pq, (nuint)q.Length, pt, (nuint)t.Length, pm, gapOpen, gapExtend);
    }

    public static int LocalScoreBatch(ReadOnlySpan<byte> q, ReadOnlySpan<byte> targets, ReadOnlySpan<int> offsets,
        ReadOnlySpan<int> matrix, int gapOpen, int gapExtend, Span<int> scores)
    {
        fixed (byte* pq = q) fixed (byte* pt = targets) fixed (int* po = offsets) fixed (int* pm = matrix) fixed (int* ps = scores)
            return bn_local_score_batch(pq, (nuint)q.Length, pt, (nuint)targets.Length, po, (nuint)(offsets.Length - 1),
                pm, gapOpen, gapExtend, ps);
    }

    public static double KabschRmsd(ReadOnlySpan<double> a, ReadOnlySpan<double> b, int n, Span<double> rotation)
    {
        fixed (double* pa = a) fixed (double* pb = b) fixed (double* pr = rotation)
            return bn_kabsch_rmsd(pa, pb, (nuint)n, rotation.IsEmpty ? null : pr);
    }
}
