using Biotech.Net.Data.Ensembl;
using Biotech.Net.Data.Http;
using Biotech.Net.Data.Ncbi;
using Biotech.Net.Data.Pdb;
using Biotech.Net.Data.UniProt;

namespace Biotech.Net.Data;

/// <summary>
/// One entry point to every connector, sharing one <see cref="HttpClient"/>, one cache and one set of options.
/// </summary>
/// <example>
/// <code>
/// var data = new BioDataClient();
/// var hbb = await data.Entrez.FetchSequenceAsync("nuccore", "NM_000518.5", EntrezSequenceFormat.GenBank);
/// var crambin = await data.UniProt.GetProteinAsync("P01542");
/// var structure = await data.Pdb.GetStructureAsync("1CRN");
/// </code>
/// </example>
public sealed class BioDataClient
{
    /// <param name="httpClient">An <see cref="HttpClient"/> to share (e.g. from <c>IHttpClientFactory</c>); null uses a shared one.</param>
    /// <param name="options">Cache, retry and rate-limit settings; null uses <see cref="DataClientOptions.Default"/>.</param>
    public BioDataClient(HttpClient? httpClient = null, DataClientOptions? options = null)
    {
        Options = options ?? DataClientOptions.Default;
        Entrez = new EntrezClient(httpClient, Options);
        UniProt = new UniProtClient(httpClient, Options);
        Pdb = new PdbClient(httpClient, Options);
        AlphaFold = new AlphaFoldClient(httpClient, Options);
        Ensembl = new EnsemblClient(httpClient, Options);
    }

    /// <summary>The options every connector uses.</summary>
    public DataClientOptions Options { get; }

    /// <summary>NCBI E-utilities (credentials from <c>NCBI_API_KEY</c> / <c>NCBI_EMAIL</c>).</summary>
    public EntrezClient Entrez { get; }

    /// <summary>UniProt REST.</summary>
    public UniProtClient UniProt { get; }

    /// <summary>RCSB PDB files, Data API and Search API.</summary>
    public PdbClient Pdb { get; }

    /// <summary>AlphaFold DB predictions.</summary>
    public AlphaFoldClient AlphaFold { get; }

    /// <summary>Ensembl REST (current release).</summary>
    public EnsemblClient Ensembl { get; }
}
