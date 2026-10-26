/**
 * Fetch a gzip'd binary data file from public/data/ and return the decoded bytes.
 * Hosts that cannot serve binary files get the same gzip bytes as base64 text (*.b64.txt;
 * see scripts/make-artifact.py). Servers that already decoded the gzip are handled too.
 */
export async function fetchBin(url: string): Promise<ArrayBuffer> {
  const res = await fetch(url);
  if (!res.ok) throw new Error(`${res.status} ${url}`);
  const buf = url.endsWith('.b64.txt')
    ? Uint8Array.from(atob((await res.text()).trim()), (c) => c.charCodeAt(0)).buffer
    : await res.arrayBuffer();
  const b = new Uint8Array(buf);
  if (b[0] !== 0x1f || b[1] !== 0x8b) return buf; // server already decoded it
  const stream = new Blob([buf]).stream().pipeThrough(new DecompressionStream('gzip'));
  return new Response(stream).arrayBuffer();
}
