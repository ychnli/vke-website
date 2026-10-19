/**
 * Resolve a path to a file in `public/` against the site's base path, so
 * links keep working when the site is deployed to a sub-directory.
 *   asset('figures/fig1.png') -> '/~user/ar/figures/fig1.png'
 * Absolute URLs (http:, https:, data:) are returned unchanged.
 */
export function asset(path: string): string {
  if (/^(https?:|data:)/.test(path)) return path;
  const base = import.meta.env.BASE_URL.replace(/\/$/, '');
  return `${base}/${path.replace(/^\//, '')}`;
}
