"""Documented allele lists survive shell globbing and reach MHCtools intact."""

import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

import pytest


@pytest.mark.parametrize('shell', ['bash', 'zsh'])
@pytest.mark.parametrize('separator', [',', ' '])
def test_quoted_allele_lists_in_shells(shell, separator, tmp_path):
    executable = shutil.which(shell)
    if executable is None:
        pytest.skip(f'{shell} is not installed')
    # These would match the unquoted examples and corrupt their allele values.
    for name in ('HLA-Ax02:01', 'HLA-By07:02', 'HLA-Ax02:01,HLA-By07:02'):
        (tmp_path / name).touch()
    code = (
        'import json, sys; '
        'from vaxrank.cli import make_vaxrank_arg_parser; '
        'from mhctools.cli import mhc_alleles_from_args; '
        'args = make_vaxrank_arg_parser().parse_args(sys.argv[1:]); '
        'print(json.dumps(mhc_alleles_from_args(args)))'
    )
    alleles = ['HLA-A*02:01', 'HLA-B*07:02']
    command = shlex.join([
        sys.executable, '-c', code, '--vcf', 'tumor.vcf', '--bam', 'tumor.bam',
        '--mhc-predictor', 'random', '--mhc-alleles', separator.join(alleles),
        '--output-csv', 'output.csv',
    ])
    env = os.environ.copy()
    env.pop('BASH_ENV', None)
    env.pop('ENV', None)
    repo = str(Path(__file__).resolve().parents[1])
    env['PYTHONPATH'] = os.pathsep.join(filter(None, [repo, env.get('PYTHONPATH')]))
    flags = ['--noprofile', '--norc'] if shell == 'bash' else ['-f']
    result = subprocess.run(
        [executable, *flags, '-c', command], cwd=tmp_path, env=env,
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout) == alleles
