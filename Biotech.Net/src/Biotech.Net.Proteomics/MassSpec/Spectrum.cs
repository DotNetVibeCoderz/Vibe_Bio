using System.Globalization;
using System.IO.Compression;
using System.Xml;

namespace Biotech.Net.Proteomics.MassSpec;

/// <summary>A centroided mass spectrum: parallel m/z and intensity arrays sorted by m/z, plus precursor information.</summary>
public sealed class Spectrum
{
    /// <summary>Creates a spectrum; peaks are sorted by m/z if they are not already.</summary>
    public Spectrum(string id, double[] mz, double[] intensity)
    {
        ArgumentNullException.ThrowIfNull(id);
        ArgumentNullException.ThrowIfNull(mz);
        ArgumentNullException.ThrowIfNull(intensity);
        if (mz.Length != intensity.Length) throw new ArgumentException("m/z and intensity arrays differ in length.");
        for (int i = 1; i < mz.Length; i++)
            if (mz[i] < mz[i - 1]) { Array.Sort(mz, intensity); break; }
        Id = id;
        Mz = mz;
        Intensity = intensity;
    }

    /// <summary>Spectrum identifier (MGF TITLE or mzML id).</summary>
    public string Id { get; }

    /// <summary>Peak m/z values, ascending.</summary>
    public double[] Mz { get; }

    /// <summary>Peak intensities aligned with <see cref="Mz"/>.</summary>
    public double[] Intensity { get; }

    /// <summary>Number of peaks.</summary>
    public int Count => Mz.Length;

    /// <summary>MS level (1 = survey, 2 = MS/MS).</summary>
    public int MsLevel { get; init; } = 2;

    /// <summary>Retention time in seconds (NaN when unknown).</summary>
    public double RetentionTime { get; init; } = double.NaN;

    /// <summary>Precursor m/z (NaN for MS1).</summary>
    public double PrecursorMz { get; init; } = double.NaN;

    /// <summary>Precursor charge (0 when unknown).</summary>
    public int PrecursorCharge { get; init; }

    /// <summary>Precursor intensity (NaN when unknown); used for intensity-based quantification.</summary>
    public double PrecursorIntensity { get; init; } = double.NaN;

    /// <summary>Total ion current.</summary>
    public double TotalIonCurrent
    {
        get
        {
            double s = 0;
            foreach (double x in Intensity) s += x;
            return s;
        }
    }

    /// <summary>Index of the most intense peak within ±<paramref name="tolerance"/> of <paramref name="mz"/>, or −1.</summary>
    public int FindPeak(double mz, double tolerance)
    {
        int i = Array.BinarySearch(Mz, mz - tolerance);
        if (i < 0) i = ~i;
        int best = -1;
        for (; i < Mz.Length && Mz[i] <= mz + tolerance; i++)
            if (best < 0 || Intensity[i] > Intensity[best]) best = i;
        return best;
    }

    /// <summary>A copy keeping the <paramref name="n"/> most intense peaks (m/z order preserved).</summary>
    public Spectrum TopPeaks(int n)
    {
        if (n >= Count) return this;
        var idx = new int[Count];
        for (int i = 0; i < idx.Length; i++) idx[i] = i;
        var keys = (double[])Intensity.Clone();
        Array.Sort(keys, idx);
        var keep = idx.AsSpan(Count - n).ToArray();
        Array.Sort(keep);
        var mz = new double[n];
        var it = new double[n];
        for (int i = 0; i < n; i++) { mz[i] = Mz[keep[i]]; it[i] = Intensity[keep[i]]; }
        return new(Id, mz, it)
        {
            MsLevel = MsLevel, RetentionTime = RetentionTime, PrecursorMz = PrecursorMz,
            PrecursorCharge = PrecursorCharge, PrecursorIntensity = PrecursorIntensity,
        };
    }

    /// <inheritdoc/>
    public override string ToString() => $"{Id}: MS{MsLevel}, {Count} peaks" + (double.IsNaN(PrecursorMz) ? "" : $", precursor {PrecursorMz:F4} ({PrecursorCharge}+)");
}

/// <summary>Reader for Mascot Generic Format (MGF) peak lists.</summary>
public static class MgfReader
{
    /// <summary>Streams the spectra of an MGF file.</summary>
    public static IEnumerable<Spectrum> Read(TextReader reader)
    {
        ArgumentNullException.ThrowIfNull(reader);
        var ci = CultureInfo.InvariantCulture;
        string? line;
        bool inside = false;
        string title = "";
        double pepMz = double.NaN, pepInt = double.NaN, rt = double.NaN;
        int charge = 0, count = 0, lineNo = 0;
        var mz = new List<double>();
        var it = new List<double>();
        while ((line = reader.ReadLine()) is not null)
        {
            lineNo++;
            var t = line.Trim();
            if (t.Length == 0 || t[0] is '#' or ';' or '!') continue;
            if (t.Equals("BEGIN IONS", StringComparison.OrdinalIgnoreCase))
            {
                inside = true; title = ""; pepMz = pepInt = rt = double.NaN; charge = 0; mz.Clear(); it.Clear();
                continue;
            }
            if (!inside) continue;
            if (t.Equals("END IONS", StringComparison.OrdinalIgnoreCase))
            {
                inside = false;
                count++;
                yield return new Spectrum(title.Length > 0 ? title : $"spectrum{count}", mz.ToArray(), it.ToArray())
                {
                    MsLevel = 2, PrecursorMz = pepMz, PrecursorIntensity = pepInt, PrecursorCharge = charge, RetentionTime = rt,
                };
                continue;
            }
            int eq = t.IndexOf('=');
            if (eq > 0 && char.IsLetter(t[0]))
            {
                string key = t[..eq].ToUpperInvariant(), value = t[(eq + 1)..].Trim();
                switch (key)
                {
                    case "TITLE": title = value; break;
                    case "PEPMASS":
                        var pm = value.Split([' ', '\t'], StringSplitOptions.RemoveEmptyEntries);
                        pepMz = double.Parse(pm[0], ci);
                        if (pm.Length > 1) pepInt = double.Parse(pm[1], ci);
                        break;
                    case "CHARGE": charge = ParseCharge(value); break;
                    case "RTINSECONDS": rt = double.Parse(value.Split('-')[0], ci); break;
                }
                continue;
            }
            var parts = t.Split([' ', '\t'], StringSplitOptions.RemoveEmptyEntries);
            if (parts.Length < 2 || !double.TryParse(parts[0], NumberStyles.Float, ci, out double m) || !double.TryParse(parts[1], NumberStyles.Float, ci, out double i))
                throw new InvalidDataException($"MGF line {lineNo}: cannot parse peak '{t}'.");
            mz.Add(m);
            it.Add(i);
        }
        if (inside) throw new InvalidDataException("MGF ended inside a BEGIN IONS block.");
    }

    /// <summary>Reads all spectra of an MGF file.</summary>
    public static List<Spectrum> ReadFile(string path)
    {
        using var r = new StreamReader(path);
        return Read(r).ToList();
    }

    private static int ParseCharge(string value)
    {
        // "2+", "3", "2+ and 3+" (first wins), "1-".
        var first = value.Split([' ', ','], StringSplitOptions.RemoveEmptyEntries)[0];
        int sign = first.EndsWith('-') ? -1 : 1;
        return sign * int.Parse(first.TrimEnd('+', '-'), CultureInfo.InvariantCulture);
    }
}

/// <summary>
/// A basic streaming mzML 1.1 reader: spectrum id, MS level, scan start time, precursor m/z / charge / intensity and the
/// m/z and intensity binary arrays (base64, 32- or 64-bit float, zlib or no compression). Numpress and indexed
/// offsets are not used; chromatograms are skipped.
/// </summary>
public static class MzMLReader
{
    /// <summary>Streams the spectra of an mzML document.</summary>
    public static IEnumerable<Spectrum> Read(Stream stream)
    {
        ArgumentNullException.ThrowIfNull(stream);
        var settings = new XmlReaderSettings { IgnoreWhitespace = true, IgnoreComments = true, DtdProcessing = DtdProcessing.Ignore };
        using var xml = XmlReader.Create(stream, settings);
        while (xml.Read())
        {
            if (xml.NodeType != XmlNodeType.Element || xml.LocalName != "spectrum") continue;
            using var sub = xml.ReadSubtree();
            yield return ReadSpectrum(sub);
        }
    }

    /// <summary>Reads all spectra of an mzML file.</summary>
    public static List<Spectrum> ReadFile(string path)
    {
        using var s = File.OpenRead(path);
        return Read(s).ToList();
    }

    private static Spectrum ReadSpectrum(XmlReader x)
    {
        var ci = CultureInfo.InvariantCulture;
        x.Read();
        string id = x.GetAttribute("id") ?? "";
        int msLevel = 1, charge = 0;
        double rt = double.NaN, precMz = double.NaN, precInt = double.NaN;
        double[]? mz = null, intensity = null;

        // Binary array state.
        bool inBinaryArray = false, is64 = true, zlib = false;
        int arrayKind = 0; // 1 = m/z, 2 = intensity
        bool inPrecursor = false;

        while (x.Read())
        {
            if (x.NodeType == XmlNodeType.EndElement)
            {
                if (x.LocalName == "binaryDataArray") inBinaryArray = false;
                else if (x.LocalName == "precursor") inPrecursor = false;
                continue;
            }
            if (x.NodeType != XmlNodeType.Element) continue;
            switch (x.LocalName)
            {
                case "precursor": inPrecursor = true; break;
                case "binaryDataArray": inBinaryArray = true; is64 = true; zlib = false; arrayKind = 0; break;
                case "cvParam":
                    string acc = x.GetAttribute("accession") ?? "";
                    string val = x.GetAttribute("value") ?? "";
                    if (inBinaryArray)
                    {
                        switch (acc)
                        {
                            case "MS:1000521": is64 = false; break;         // 32-bit float
                            case "MS:1000523": is64 = true; break;          // 64-bit float
                            case "MS:1000574": zlib = true; break;          // zlib compression
                            case "MS:1000576": zlib = false; break;         // no compression
                            case "MS:1000514": arrayKind = 1; break;        // m/z array
                            case "MS:1000515": arrayKind = 2; break;        // intensity array
                            default:
                                if (acc.StartsWith("MS:100261", StringComparison.Ordinal) || acc == "MS:1002312" || acc == "MS:1002313" || acc == "MS:1002314")
                                    throw new NotSupportedException("MS-Numpress compressed arrays are not supported.");
                                break;
                        }
                    }
                    else
                    {
                        switch (acc)
                        {
                            case "MS:1000511": msLevel = int.Parse(val, ci); break;                 // ms level
                            case "MS:1000016":                                                       // scan start time
                                rt = double.Parse(val, ci);
                                if ((x.GetAttribute("unitName") ?? "").StartsWith("min", StringComparison.OrdinalIgnoreCase)
                                    || x.GetAttribute("unitAccession") == "UO:0000031") rt *= 60;
                                break;
                            case "MS:1000744": if (inPrecursor) precMz = double.Parse(val, ci); break; // selected ion m/z
                            case "MS:1000041": if (inPrecursor) charge = int.Parse(val, ci); break;    // charge state
                            case "MS:1000042": if (inPrecursor) precInt = double.Parse(val, ci); break; // peak intensity
                        }
                    }
                    break;
                case "binary":
                    if (!inBinaryArray) break;
                    string b64 = x.ReadElementContentAsString();
                    var values = Decode(b64, is64, zlib);
                    if (arrayKind == 1) mz = values; else if (arrayKind == 2) intensity = values;
                    break;
            }
        }
        mz ??= [];
        intensity ??= new double[mz.Length];
        return new Spectrum(id, mz, intensity)
        {
            MsLevel = msLevel, RetentionTime = rt, PrecursorMz = precMz, PrecursorCharge = charge, PrecursorIntensity = precInt,
        };
    }

    /// <summary>Decodes an mzML binary array (base64, optionally zlib-compressed, little-endian floats).</summary>
    public static double[] Decode(string base64, bool is64Bit, bool zlibCompressed)
    {
        var bytes = Convert.FromBase64String(base64.Trim());
        if (zlibCompressed)
        {
            using var input = new MemoryStream(bytes);
            using var z = new ZLibStream(input, CompressionMode.Decompress);
            using var output = new MemoryStream();
            z.CopyTo(output);
            bytes = output.ToArray();
        }
        int width = is64Bit ? 8 : 4;
        if (bytes.Length % width != 0) throw new InvalidDataException("Binary array length is not a multiple of the value width.");
        var values = new double[bytes.Length / width];
        for (int i = 0; i < values.Length; i++)
            values[i] = is64Bit
                ? BitConverter.ToDouble(bytes, i * 8)
                : BitConverter.ToSingle(bytes, i * 4);
        if (!BitConverter.IsLittleEndian) throw new PlatformNotSupportedException("Big-endian hosts are not supported.");
        return values;
    }

    /// <summary>Encodes values the way mzML stores them (used by writers and tests).</summary>
    public static string Encode(ReadOnlySpan<double> values, bool is64Bit, bool zlibCompressed)
    {
        var bytes = new byte[values.Length * (is64Bit ? 8 : 4)];
        for (int i = 0; i < values.Length; i++)
        {
            if (is64Bit) BitConverter.TryWriteBytes(bytes.AsSpan(i * 8), values[i]);
            else BitConverter.TryWriteBytes(bytes.AsSpan(i * 4), (float)values[i]);
        }
        if (zlibCompressed)
        {
            using var output = new MemoryStream();
            using (var z = new ZLibStream(output, CompressionLevel.Optimal, leaveOpen: true)) z.Write(bytes);
            bytes = output.ToArray();
        }
        return Convert.ToBase64String(bytes);
    }
}
