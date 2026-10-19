// Biotech.Net Web Gallery: small browser helpers called through JS interop.
// Dibuat oleh Gravicode Studios dipimpin oleh Kang Fadhil.
window.gallery = {
  copy: async function (text) {
    try { await navigator.clipboard.writeText(text); return true; } catch { return false; }
  },

  download: function (fileName, mime, text) {
    const blob = new Blob([text], { type: mime });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url; a.download = fileName;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  },

  // Code editor niceties for a plain <textarea>: Tab indents, Ctrl/Cmd+Enter runs.
  attachEditor: function (el, dotnet) {
    if (!el || el.dataset.attached) return;
    el.dataset.attached = '1';
    el.addEventListener('keydown', (e) => {
      if (e.key === 'Tab' && !e.shiftKey && !e.ctrlKey && !e.metaKey) {
        e.preventDefault();
        const s = el.selectionStart, t = el.selectionEnd;
        el.setRangeText('    ', s, t, 'end');
        el.dispatchEvent(new Event('input', { bubbles: true }));
      } else if (e.key === 'Enter' && (e.ctrlKey || e.metaKey)) {
        e.preventDefault();
        el.dispatchEvent(new Event('change', { bubbles: true }));
        dotnet.invokeMethodAsync('RunFromKeyboard', el.value);
      }
    });
    const sync = () => {
      const gutter = el.parentElement && el.parentElement.querySelector('.gutter');
      if (gutter) gutter.scrollTop = el.scrollTop;
    };
    el.addEventListener('scroll', sync);
  },

  setValue: function (el, value) { if (el) { el.value = value; el.scrollTop = 0; } },

  // Lets the browser paint "computing…" before a long synchronous WebAssembly call.
  nextFrame: function () { return new Promise((r) => requestAnimationFrame(() => setTimeout(r, 0))); },

  scrollTop: function () { const m = document.querySelector('.stage'); if (m) m.scrollTop = 0; window.scrollTo(0, 0); }
};
