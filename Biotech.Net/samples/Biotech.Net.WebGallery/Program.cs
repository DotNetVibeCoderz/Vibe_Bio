using Biotech.Net.WebGallery;
using Biotech.Net.WebGallery.Services;
using Microsoft.AspNetCore.Components.Web;
using Microsoft.AspNetCore.Components.WebAssembly.Hosting;
using Microsoft.AspNetCore.Components.WebAssembly.Services;

// Biotech.Net Web Gallery. Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil.
var builder = WebAssemblyHostBuilder.CreateDefault(args);
builder.RootComponents.Add<App>("#app");
builder.RootComponents.Add<HeadOutlet>("head::after");

builder.Services.AddScoped(_ => new HttpClient { BaseAddress = new Uri(builder.HostEnvironment.BaseAddress) });
builder.Services.AddScoped<LazyAssemblyLoader>();
builder.Services.AddSingleton<PlaygroundState>();

await builder.Build().RunAsync();
