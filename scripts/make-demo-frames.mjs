// Generates placeholder SVG frames in public/animations/demo/ so the
// FrameAnimation component can be tested before real frames exist.
// Delete public/animations/demo/ once real animations are in place.
import { mkdirSync, writeFileSync } from 'node:fs';

const dir = new URL('../public/animations/demo/', import.meta.url);
mkdirSync(dir, { recursive: true });

const N = 24;
for (let i = 0; i < N; i++) {
  const t = i / N;
  const x = 80 + 440 * t;
  const y = 230 - 120 * t;
  const svg = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 640 320" width="640" height="320">
  <rect width="640" height="320" fill="#dfe9f3"/>
  <ellipse cx="${x}" cy="${y}" rx="150" ry="28" transform="rotate(-18 ${x} ${y})" fill="#2a6fb0" opacity="0.75"/>
  <text x="20" y="300" font-family="sans-serif" font-size="16" fill="#333">Placeholder frame ${i + 1}/${N}</text>
</svg>`;
  writeFileSync(new URL(`frame_${String(i).padStart(3, '0')}.svg`, dir), svg);
}
console.log(`Wrote ${N} frames to ${dir.pathname}`);
