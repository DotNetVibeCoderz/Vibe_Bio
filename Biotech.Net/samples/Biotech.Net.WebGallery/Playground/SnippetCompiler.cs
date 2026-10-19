using System.Diagnostics;
using System.Reflection;
using System.Reflection.Metadata;
using Microsoft.CodeAnalysis;
using Microsoft.CodeAnalysis.CSharp;

namespace Biotech.Net.WebGallery.Playground;

public sealed record SnippetDiagnostic(int Line, int Column, string Severity, string Id, string Message);

public sealed record SnippetRun(
    bool Compiled,
    IReadOnlyList<SnippetDiagnostic> Diagnostics,
    string Console,
    IReadOnlyList<SnippetOutput> Outputs,
    string? Exception,
    TimeSpan CompileTime,
    TimeSpan RunTime);

/// <summary>
/// Compiles a C# snippet (top-level statements) with Roslyn inside the browser and runs it. References are the
/// assemblies already loaded into the WebAssembly runtime, read straight from memory with
/// <see cref="AssemblyExtensions.TryGetRawMetadata"/>: no reference assemblies are downloaded.
/// </summary>
/// <remarks>
/// Only touched after the Playground has lazy-loaded Microsoft.CodeAnalysis, so no Roslyn type appears in any
/// signature that the rest of the app compiles against.
/// </remarks>
public static class SnippetCompiler
{
    /// <summary>Assemblies a snippet can use; loaded on demand so their metadata is available.</summary>
    private static readonly string[] Preload =
    [
        "System.Runtime", "System.Private.CoreLib", "System.Console", "System.Linq", "System.Collections",
        "System.Collections.Immutable", "System.Runtime.Extensions", "System.Text.RegularExpressions", "System.Memory",
        "netstandard",
        "Biotech.Net.Native", "Biotech.Net.Core", "Biotech.Net.Sequences", "Biotech.Net.Alignments", "Biotech.Net.Genomics",
        "Biotech.Net.Transcriptomics", "Biotech.Net.Proteomics", "Biotech.Net.StructuralBiology", "Biotech.Net.Chemistry",
        "Biotech.Net.DrugDiscovery", "Biotech.Net.SyntheticBiology", "Biotech.Net.Visualization",
        "Visual.Net.Core", "Matplotlib.Net", "Biotech.Net.WebGallery",
    ];

    /// <summary>Usings every snippet gets, like the SDK's implicit usings plus the Biotech.Net core and Show().</summary>
    public const string GlobalUsings = """
        global using System;
        global using System.Collections.Generic;
        global using System.IO;
        global using System.Linq;
        global using System.Text;
        global using System.Threading.Tasks;
        global using Biotech.Net.Core;
        global using static Biotech.Net.WebGallery.Playground.Display;
        """;

    private static List<MetadataReference>? _references;
    private static readonly CSharpParseOptions ParseOptions = new(LanguageVersion.Latest);

    /// <summary>How many assemblies the compiler references (after <see cref="Warm"/>).</summary>
    public static int ReferenceCount => _references?.Count ?? 0;

    /// <summary>Builds the reference set and compiles a tiny program once, so the first real run is quick.</summary>
    public static TimeSpan Warm()
    {
        var sw = Stopwatch.StartNew();
        References();
        Compile("Console.WriteLine(1);", out _);
        return sw.Elapsed;
    }

    private static List<MetadataReference> References()
    {
        if (_references is not null) return _references;
        foreach (var name in Preload)
        {
            try { Assembly.Load(name); } catch { /* not part of this build */ }
        }
        var refs = new List<MetadataReference>();
        var seen = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
        foreach (var asm in AppDomain.CurrentDomain.GetAssemblies())
        {
            if (asm.IsDynamic || asm.GetName().Name is not { } n || !seen.Add(n)) continue;
            if (n.StartsWith("Snippet", StringComparison.Ordinal)) continue;
            unsafe
            {
                if (asm.TryGetRawMetadata(out byte* blob, out int length))
                {
                    var module = ModuleMetadata.CreateFromMetadata((nint)blob, length);
                    refs.Add(AssemblyMetadata.Create(module).GetReference(display: n));
                }
            }
        }
        return _references = refs;
    }

    private static CSharpCompilation Compilation(string code)
    {
        var trees = new[]
        {
            CSharpSyntaxTree.ParseText(GlobalUsings, ParseOptions, path: "GlobalUsings.cs"),
            CSharpSyntaxTree.ParseText(code, ParseOptions, path: "Snippet.cs"),
        };
        var options = new CSharpCompilationOptions(OutputKind.ConsoleApplication,
            optimizationLevel: OptimizationLevel.Release,
            allowUnsafe: true,
            concurrentBuild: false,            // the browser runtime is single-threaded
            nullableContextOptions: NullableContextOptions.Annotations);   // `string?` allowed, no flow warnings
        return CSharpCompilation.Create("Snippet" + Guid.NewGuid().ToString("N")[..8], trees, References(), options);
    }

    private static byte[]? Compile(string code, out List<SnippetDiagnostic> diagnostics)
    {
        var compilation = Compilation(code);
        using var pe = new MemoryStream();
        var result = compilation.Emit(pe);
        diagnostics = result.Diagnostics
            .Where(d => d.Severity >= DiagnosticSeverity.Warning && d.Location.SourceTree?.FilePath == "Snippet.cs")
            .Concat(result.Diagnostics.Where(d => d.Severity == DiagnosticSeverity.Error && d.Location.SourceTree?.FilePath != "Snippet.cs"))
            .OrderByDescending(d => d.Severity)
            .Select(d =>
            {
                var span = d.Location.GetLineSpan();
                return new SnippetDiagnostic(span.StartLinePosition.Line + 1, span.StartLinePosition.Character + 1,
                    d.Severity.ToString().ToLowerInvariant(), d.Id, d.GetMessage(System.Globalization.CultureInfo.InvariantCulture));
            })
            .Where(d => d.Id != "CS8019")
            .ToList();
        return result.Success ? pe.ToArray() : null;
    }

    /// <summary>Compiles and runs <paramref name="code"/>, capturing Console output and <c>Show(...)</c> outputs.</summary>
    public static async Task<SnippetRun> RunAsync(string code)
    {
        var compileWatch = Stopwatch.StartNew();
        byte[]? image;
        List<SnippetDiagnostic> diagnostics;
        try
        {
            image = Compile(code, out diagnostics);
        }
        catch (Exception ex)
        {
            return new SnippetRun(false, [new(0, 0, "error", ex.GetType().Name, ex.Message)], "", [], null, compileWatch.Elapsed, TimeSpan.Zero);
        }
        compileWatch.Stop();
        if (image is null) return new SnippetRun(false, diagnostics, "", [], null, compileWatch.Elapsed, TimeSpan.Zero);

        var console = new StringWriter();
        var outputs = new List<SnippetOutput>();
        var original = Console.Out;
        string? exception = null;
        var runWatch = Stopwatch.StartNew();
        try
        {
            Console.SetOut(console);
            Display.Sink = outputs;
            var assembly = Assembly.Load(image);
            var entry = assembly.EntryPoint ?? throw new InvalidOperationException("The snippet has no entry point.");
            object?[]? args = entry.GetParameters().Length == 0 ? null : [Array.Empty<string>()];
            var result = entry.Invoke(null, args);
            if (result is Task task) await task;
        }
        catch (Exception ex)
        {
            var inner = ex is TargetInvocationException { InnerException: { } ie } ? ie : ex;
            exception = $"{inner.GetType().FullName}: {inner.Message}";
            var frame = inner.StackTrace?.Split('\n').FirstOrDefault(l => l.Contains("Program.<Main>", StringComparison.Ordinal) || l.Contains("<Main>$", StringComparison.Ordinal));
            if (frame is not null) exception += "\n   " + frame.Trim();
        }
        finally
        {
            Console.SetOut(original);
            Display.Sink = null;
        }
        runWatch.Stop();
        return new SnippetRun(true, diagnostics, console.ToString(), outputs, exception, compileWatch.Elapsed, runWatch.Elapsed);
    }
}
