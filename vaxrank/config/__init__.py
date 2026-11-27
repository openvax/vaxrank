"""Configuration loading and defaults, imported lazily to keep CLI startup light."""

from importlib import import_module

_LAZY_EXPORTS = {
    "extract_epitope_config_kwargs": (".loader", "extract_epitope_config_kwargs"),
    "extract_vaccine_config_kwargs": (".loader", "extract_vaccine_config_kwargs"),
    "load_vaxrank_config": (".loader", "load_vaxrank_config"),
}

__all__ = list(_LAZY_EXPORTS)


def __getattr__(name):
    if name in _LAZY_EXPORTS:
        module_name, attr_name = _LAZY_EXPORTS[name]
        value = getattr(import_module(module_name, package=__name__), attr_name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | set(__all__))
