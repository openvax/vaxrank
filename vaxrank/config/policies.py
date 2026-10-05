"""Discover and inspect bundled selection definitions without biological inputs."""

from importlib.resources import files

import msgspec

from .loader import load_vaxrank_config, extract_epitope_config_kwargs
from ..epitope_config import EpitopeConfig


def list_policies():
    """Return bundled policy names and descriptions in stable name order."""
    records = []
    for path in files('vaxrank.config').iterdir():
        if path.name.endswith('.yaml') and path.name != 'default.yaml':
            definition = msgspec.yaml.decode(path.read_bytes())
            if definition.get('name'):
                records.append(dict(name=definition['name'], description=(
                    definition.get('policy_metadata', {}).get('description',
                        'OpenVax selection bundle; see --show-policy'))))
    return sorted(records, key=lambda record: record['name'])


def show_policy(name):
    """Resolve a bundled policy over the frozen baseline and current defaults.

    Returns the same configuration used by ``--config builtin:openvax-v1
    --config builtin:NAME``, with the complete canonical Topiary definition.
    This inspection does not evaluate evidence or execute predictors.
    """
    if name not in {record['name'] for record in list_policies()}:
        raise ValueError('Unknown policy %r; use --list-policies' % name)
    configuration = load_vaxrank_config(
        config_path=['builtin:openvax-v1', 'builtin:' + name])
    epitope_config = msgspec.convert(extract_epitope_config_kwargs(configuration), EpitopeConfig)
    configuration['epitopes']['selection_policy'] = epitope_config.selection_policy
    return configuration
