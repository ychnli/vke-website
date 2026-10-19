// @ts-check
import { defineConfig } from 'astro/config';
import mdx from '@astrojs/mdx';
import remarkMath from 'remark-math';
import rehypeKatex from 'rehype-katex';
import { rehypeHeadingIds } from '@astrojs/markdown-remark';
import rehypeAutolinkHeadings from 'rehype-autolink-headings';

// The site can be served from a sub-path (e.g. https://web.stanford.edu/~user/ar/).
// Set SITE_URL and SITE_BASE at build time when the final host is known, e.g.
//   SITE_URL=https://web.stanford.edu SITE_BASE=/~user/ar npm run build
export default defineConfig({
  site: process.env.SITE_URL,
  base: process.env.SITE_BASE || '/',  // GitHub Pages passes '' for a user or custom-domain site
  // PALETTE=mono builds the monochrome palette into its own folder.
  outDir: process.env.PALETTE === 'mono' ? './dist-mono' : './dist',
  // Emit ./page/index.html so the build works on plain static file hosts.
  trailingSlash: 'ignore',
  build: { format: 'directory' },
  integrations: [mdx()],
  markdown: {
    // LaTeX: $inline$ and $$display$$ are rendered to HTML by KaTeX at build time.
    remarkPlugins: [remarkMath],
    // Section headings get ids, and a clickable "#" anchor before the title.
    rehypePlugins: [
      // Budget-term macros. <kind> is one of the term kinds defined in src/styles/global.css:
      //   \term{pe}{...}  colored and clickable (opens the matching <Term kind="pe"> panel)
      //   \tc{pe}{...}    colored only
      // Colors are set in CSS, never in the math, so they can be changed in one place.
      [
        rehypeKatex,
        {
          macros: {
            '\\term': '\\htmlData{term=#1}{#2}',
            '\\tc': '\\htmlData{color=#1}{#2}',
          },
          trust: (ctx) => ctx.command === '\\htmlData',
          strict: (code) => (code === 'htmlExtension' ? 'ignore' : 'warn'),
        },
      ],
      rehypeHeadingIds,
      [
        rehypeAutolinkHeadings,
        {
          behavior: 'prepend',
          test: ['h2', 'h3'],
          properties: { className: ['heading-anchor'], ariaHidden: 'true', tabIndex: -1 },
          content: { type: 'text', value: '#' },
        },
      ],
    ],
  },
});
