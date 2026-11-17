"""Use PyPI-compatible README links while checking their local MkDocs targets."""

from pathlib import Path
import re


root = Path(__file__).resolve().parents[1]


def local_link(match):
    route, anchor = match.groups()
    target = route.rstrip('/') + '.md' if route else 'index.md'
    if target != 'index.md' and not (root / 'docs' / target).is_file():
        raise ValueError('README documentation target does not exist: ' + target)
    return '](' + target + (anchor or '') + ')'


readme = (root / 'README.md').read_text()
readme = re.sub(r'\]\(https://openvax\.github\.io/vaxrank/([^#)]*)(#[^)]*)?\)', local_link, readme)
(root / 'docs' / 'index.md').write_text(readme.replace('](docs/', ']('))
