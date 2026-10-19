using System.Net;
using System.Text;

namespace Biotech.Net.WebGallery.Shared;

/// <summary>A tiny line-at-a-time C# highlighter: keywords, types, strings, comments and numbers as spans.</summary>
public static class CSharpHighlighter
{
    private static readonly HashSet<string> Keywords =
    [
        "using", "var", "new", "foreach", "for", "in", "if", "else", "return", "true", "false", "null", "string", "int",
        "double", "bool", "long", "char", "byte", "void", "class", "record", "static", "public", "private", "async",
        "await", "is", "not", "and", "or", "switch", "case", "break", "const", "readonly", "while", "do", "out", "ref",
        "this", "typeof", "nameof", "default", "throw", "try", "catch", "finally", "struct", "object", "decimal", "float",
        "sealed", "namespace", "with", "when", "select", "from", "where", "orderby", "let", "uint", "ulong",
    ];

    public static string Line(string line)
    {
        var sb = new StringBuilder(line.Length + 32);
        int i = 0;
        while (i < line.Length)
        {
            char c = line[i];
            if (c == '/' && i + 1 < line.Length && line[i + 1] == '/')
            {
                Span(sb, "cm", line[i..]);
                break;
            }
            if (c == '"' || (c == '$' && i + 1 < line.Length && line[i + 1] == '"') || (c == '@' && i + 1 < line.Length && line[i + 1] == '"'))
            {
                int j = c == '"' ? i + 1 : i + 2;
                while (j < line.Length && line[j] != '"') { if (line[j] == '\\') j++; j++; }
                j = Math.Min(line.Length, j + 1);
                Span(sb, "st", line[i..j]);
                i = j;
                continue;
            }
            if (c == '\'' && i + 2 < line.Length)
            {
                int j = line.IndexOf('\'', i + 1 + (line[i + 1] == '\\' ? 2 : 1));
                if (j > i && j - i <= 4) { Span(sb, "st", line[i..(j + 1)]); i = j + 1; continue; }
            }
            if (char.IsDigit(c) && (i == 0 || !char.IsLetterOrDigit(line[i - 1]) && line[i - 1] != '_'))
            {
                int j = i;
                while (j < line.Length && (char.IsLetterOrDigit(line[j]) || line[j] is '.' or '_')) j++;
                Span(sb, "nu", line[i..j]);
                i = j;
                continue;
            }
            if (char.IsLetter(c) || c == '_')
            {
                int j = i;
                while (j < line.Length && (char.IsLetterOrDigit(line[j]) || line[j] == '_')) j++;
                string word = line[i..j];
                if (Keywords.Contains(word)) Span(sb, "kw", word);
                else if (char.IsUpper(word[0])) Span(sb, "ty", word);
                else sb.Append(WebUtility.HtmlEncode(word));
                i = j;
                continue;
            }
            sb.Append(WebUtility.HtmlEncode(c.ToString()));
            i++;
        }
        return sb.ToString();
    }

    private static void Span(StringBuilder sb, string cls, string text)
        => sb.Append("<span class=\"").Append(cls).Append("\">").Append(WebUtility.HtmlEncode(text)).Append("</span>");
}
