using System.Text.Json;
using Biotech.Net.Data.Http;

namespace Biotech.Net.Data.Ncbi;

/// <summary>Record formats that <see cref="EntrezClient.FetchSequencesAsync(string, IEnumerable{string}, EntrezSequenceFormat, CancellationToken)"/> parses.</summary>
public enum EntrezSequenceFormat
{
    /// <summary><c>rettype=fasta</c>: identifier, description and residues.</summary>
    Fasta,

    /// <summary><c>rettype=gb</c> (GenBank flat file; GenPept for proteins): annotations and the feature table as well.</summary>
    GenBank,
}

/// <summary>A result set kept on the NCBI history server (<c>WebEnv</c> + <c>query_key</c>).</summary>
public sealed record EntrezHistory(string WebEnv, string QueryKey);

/// <summary>The answer to an ESearch call.</summary>
/// <param name="Database">The Entrez database searched.</param>
/// <param name="Count">Total number of matching records.</param>
/// <param name="RetStart">Offset of the first returned id.</param>
/// <param name="RetMax">Number of ids requested.</param>
/// <param name="Ids">The returned UIDs (at most <paramref name="RetMax"/>).</param>
/// <param name="QueryTranslation">How Entrez interpreted the term.</param>
/// <param name="History">The history-server handle when <c>usehistory</c> was requested.</param>
public sealed record ESearchResult(string Database, int Count, int RetStart, int RetMax, IReadOnlyList<string> Ids,
    string? QueryTranslation, EntrezHistory? History);

/// <summary>One ESummary document (version 0.3 JSON): the UID plus every field as JSON.</summary>
public sealed class ESummaryDocument
{
    internal ESummaryDocument(string uid, JsonElement fields)
    {
        Uid = uid;
        Fields = fields;
    }

    /// <summary>The record UID.</summary>
    public string Uid { get; }

    /// <summary>All fields of the summary (a detached clone; safe to keep).</summary>
    public JsonElement Fields { get; }

    /// <summary>A field as text (numbers are returned as their JSON text), or null.</summary>
    public string? GetString(string field) => Fields.Str(field);

    public override string ToString() => $"{Uid}: {GetString("title") ?? GetString("caption") ?? GetString("name")}";
}

/// <summary>A PubMed citation summary.</summary>
public sealed record PubMedSummary(string Pmid, string Title, string Journal, string PubDate, IReadOnlyList<string> Authors,
    string? Doi, string? Volume, string? Issue, string? Pages)
{
    public override string ToString() =>
        $"{(Authors.Count > 0 ? Authors[0] + (Authors.Count > 1 ? " et al." : "") + " " : "")}{Title} {Journal} {PubDate}. PMID {Pmid}";
}

/// <summary>The links found for one ELink input set.</summary>
public sealed record ELinkSet(string DbFrom, IReadOnlyList<string> Ids, IReadOnlyList<ELinkResult> Links);

/// <summary>Linked UIDs in one target database under one link name (e.g. <c>protein_gene</c>).</summary>
public sealed record ELinkResult(string DbTo, string LinkName, IReadOnlyList<string> Ids);

/// <summary>A searchable field of an Entrez database.</summary>
public sealed record EInfoField(string Name, string FullName, string Description);

/// <summary>A link from an Entrez database to another.</summary>
public sealed record EInfoLink(string Name, string Menu, string Description, string DbTo);

/// <summary>EInfo statistics for one Entrez database.</summary>
public sealed record EInfoDatabase(string Name, string MenuName, string Description, long Count, string? LastUpdate,
    IReadOnlyList<EInfoField> Fields, IReadOnlyList<EInfoLink> Links);
