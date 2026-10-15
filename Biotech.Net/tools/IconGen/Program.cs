// Renders the Biotech.Net icons: assets/icon.png (NuGet, 128 px), assets/icon-512.png and assets/app.ico.
// Usage: dotnet run --project tools/IconGen -- <output-directory>
// Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil.
using SkiaSharp;

#pragma warning disable CS0618 // SKPath mutators: SKPathBuilder is not worth it for a build-time tool.

string outDir = args.Length > 0 ? args[0] : "assets";
Directory.CreateDirectory(outDir);
File.WriteAllBytes(Path.Combine(outDir, "icon-512.png"), Render(512));
File.WriteAllBytes(Path.Combine(outDir, "icon.png"), Render(128));
WriteIco(Path.Combine(outDir, "app.ico"), [16, 24, 32, 48, 64, 128, 256]);
Console.WriteLine($"Icons written to {Path.GetFullPath(outDir)}");

static byte[] Render(int size)
{
    using var surface = SKSurface.Create(new SKImageInfo(size, size, SKColorType.Rgba8888, SKAlphaType.Premul));
    var c = surface.Canvas;
    c.Clear(SKColors.Transparent);
    c.Scale(size / 512f);

    // Tile: deep ocean gradient, the "lab at night" backdrop the helix glows against.
    var tile = new SKRoundRect(new SKRect(16, 16, 496, 496), 112, 112);
    using (var bg = new SKPaint { IsAntialias = true })
    {
        bg.Shader = SKShader.CreateLinearGradient(new SKPoint(16, 16), new SKPoint(496, 496),
            [SKColor.Parse("#04131F"), SKColor.Parse("#0A2E3F"), SKColor.Parse("#1A1B4B")],
            [0f, 0.55f, 1f], SKShaderTileMode.Clamp);
        c.DrawRoundRect(tile, bg);
    }

    c.Save();
    c.ClipRoundRect(tile, antialias: true);

    using (var glow = new SKPaint { IsAntialias = true })
    {
        glow.Shader = SKShader.CreateRadialGradient(new SKPoint(256, 256), 250,
            [SKColor.Parse("#2DD4BF").WithAlpha(60), SKColors.Transparent], SKShaderTileMode.Clamp);
        c.DrawRect(0, 0, 512, 512, glow);
    }

    // Faint hexagonal lattice: chemistry under the biology.
    using (var hex = new SKPaint { IsAntialias = true, Style = SKPaintStyle.Stroke, StrokeWidth = 2, Color = SKColors.White.WithAlpha(14) })
    {
        const float r = 34;
        float w = r * MathF.Sqrt(3);
        for (int row = -1; row < 11; row++)
            for (int col = -1; col < 10; col++)
            {
                float cx = col * w + (row % 2 == 0 ? 0 : w / 2), cy = row * r * 1.5f;
                c.DrawPath(Hexagon(cx, cy, r), hex);
            }
    }

    // The double helix, drawn along a diagonal so it reads at 16 px too.
    c.Save();
    c.RotateDegrees(-38, 256, 256);
    const float x0 = 50, x1 = 462, amp = 92, cy0 = 256;
    const float period = 206;
    float Y(float x, float phase) => cy0 + amp * MathF.Sin((x - x0) / period * 2 * MathF.PI + phase);

    var baseColors = new[] { "#34D399", "#F87171", "#FBBF24", "#60A5FA" }; // A T G C
    // Depth of strand A at x: +1 facing the viewer, -1 behind. Strand B is always opposite.
    float Depth(float x) => MathF.Cos((x - x0) / period * 2 * MathF.PI);

    // Back halves first, then base pairs, then front halves: the strands weave over and under.
    DrawStrand(c, x0, x1, x => Y(x, 0), x => Depth(x) < 0, SKColor.Parse("#5EEAD4"), SKColor.Parse("#22D3EE"), back: true);
    DrawStrand(c, x0, x1, x => Y(x, MathF.PI), x => Depth(x) > 0, SKColor.Parse("#A78BFA"), SKColor.Parse("#F472B6"), back: true);

    int rung = 0;
    for (float x = x0 + 14; x < x1 - 4; x += 24.5f, rung++)
    {
        float ya = Y(x, 0), yb = Y(x, MathF.PI);
        byte alpha = (byte)(150 + 90 * MathF.Abs(Depth(x)));
        float mid = (ya + yb) / 2;
        int b = rung % 4, pair = b ^ 1; // A-T, G-C
        using var p1 = new SKPaint { IsAntialias = true, StrokeWidth = 10, StrokeCap = SKStrokeCap.Round, Color = SKColor.Parse(baseColors[b]).WithAlpha(alpha) };
        using var p2 = new SKPaint { IsAntialias = true, StrokeWidth = 10, StrokeCap = SKStrokeCap.Round, Color = SKColor.Parse(baseColors[pair]).WithAlpha(alpha) };
        c.DrawLine(x, ya, x, mid, p1);
        c.DrawLine(x, mid, x, yb, p2);
    }

    DrawStrand(c, x0, x1, x => Y(x, 0), x => Depth(x) >= 0, SKColor.Parse("#5EEAD4"), SKColor.Parse("#22D3EE"), back: false);
    DrawStrand(c, x0, x1, x => Y(x, MathF.PI), x => Depth(x) <= 0, SKColor.Parse("#A78BFA"), SKColor.Parse("#F472B6"), back: false);
    c.Restore();

    // A bright node: the "Net" — data flowing out of the molecule.
    using (var node = new SKPaint { IsAntialias = true })
    {
        node.Shader = SKShader.CreateRadialGradient(new SKPoint(392, 120), 30,
            [SKColors.White, SKColor.Parse("#5EEAD4"), SKColor.Parse("#5EEAD4").WithAlpha(0)], [0f, 0.45f, 1f], SKShaderTileMode.Clamp);
        c.DrawCircle(392, 120, 30, node);
    }

    c.Restore();

    using (var border = new SKPaint { IsAntialias = true, Style = SKPaintStyle.Stroke, StrokeWidth = 4 })
    {
        border.Shader = SKShader.CreateLinearGradient(new SKPoint(16, 16), new SKPoint(496, 496),
            [SKColors.White.WithAlpha(100), SKColors.White.WithAlpha(12)], SKShaderTileMode.Clamp);
        c.DrawRoundRect(tile, border);
    }

    using var image = surface.Snapshot();
    using var data = image.Encode(SKEncodedImageFormat.Png, 100);
    return data.ToArray();
}

static void DrawStrand(SKCanvas c, float x0, float x1, Func<float, float> y, Func<float, bool> visible, SKColor from, SKColor to, bool back)
{
    // One path per visible run, so a strand can pass behind its partner and come out in front again.
    using var path = new SKPath();
    bool open = false;
    for (float x = x0; x <= x1; x += 1.5f)
    {
        if (!visible(x)) { open = false; continue; }
        if (!open) { path.MoveTo(x, y(x)); open = true; }
        else path.LineTo(x, y(x));
    }

    if (!back)
    {
        using var shadow = new SKPaint { IsAntialias = true, Style = SKPaintStyle.Stroke, StrokeWidth = 34, StrokeCap = SKStrokeCap.Round, Color = from.WithAlpha(70) };
        shadow.MaskFilter = SKMaskFilter.CreateBlur(SKBlurStyle.Normal, 12);
        c.DrawPath(path, shadow);
    }

    using var line = new SKPaint { IsAntialias = true, Style = SKPaintStyle.Stroke, StrokeWidth = back ? 18 : 22, StrokeCap = SKStrokeCap.Round };
    var a = back ? (byte)150 : (byte)255;
    line.Shader = SKShader.CreateLinearGradient(new SKPoint(x0, 0), new SKPoint(x1, 0), [from.WithAlpha(a), to.WithAlpha(a)], SKShaderTileMode.Clamp);
    c.DrawPath(path, line);

    if (back) return;
    using var shine = new SKPaint { IsAntialias = true, Style = SKPaintStyle.Stroke, StrokeWidth = 6, StrokeCap = SKStrokeCap.Round, Color = SKColors.White.WithAlpha(110) };
    c.Save();
    c.Translate(0, -5);
    c.DrawPath(path, shine);
    c.Restore();
}

static SKPath Hexagon(float cx, float cy, float r)
{
    var p = new SKPath();
    for (int i = 0; i < 6; i++)
    {
        float a = MathF.PI / 6 + i * MathF.PI / 3;
        float x = cx + r * MathF.Cos(a), y = cy + r * MathF.Sin(a);
        if (i == 0) p.MoveTo(x, y); else p.LineTo(x, y);
    }
    p.Close();
    return p;
}

// ICO with PNG-compressed entries (supported since Windows Vista).
static void WriteIco(string file, int[] sizes)
{
    var images = sizes.Select(Render).ToArray();
    using var fs = File.Create(file);
    using var w = new BinaryWriter(fs);
    w.Write((short)0); w.Write((short)1); w.Write((short)sizes.Length);
    int offset = 6 + 16 * sizes.Length;
    for (int i = 0; i < sizes.Length; i++)
    {
        byte dim = (byte)(sizes[i] >= 256 ? 0 : sizes[i]);
        w.Write(dim); w.Write(dim); w.Write((byte)0); w.Write((byte)0);
        w.Write((short)1); w.Write((short)32);
        w.Write(images[i].Length); w.Write(offset);
        offset += images[i].Length;
    }
    foreach (var img in images) w.Write(img);
}
