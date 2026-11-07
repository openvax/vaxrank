"""Varcode effect classes that native serialization accepts."""

import pytest
from varcode.effects import effect_classes

from vaxrank import native_serialization


@pytest.fixture
def fresh_registry():
    native_serialization._varcode_graph_classes.cache_clear()
    yield native_serialization._varcode_graph_classes
    native_serialization._varcode_graph_classes.cache_clear()


def test_unannotated_contig_translocations_are_allowed(fresh_registry, monkeypatch):
    # Varcode 11.9 reports a breakend mate on a contig the annotation lacks
    # (a decoy, EBV, hg19's chrM) with this subclass (openvax/varcode#709).
    cls = getattr(effect_classes, "TranslocationToUnannotatedContig", None)
    if cls is None:  # An older Varcode: stand in for it.
        cls = type("TranslocationToUnannotatedContig", (effect_classes.TranslocationToIntergenic,),
                   {"__module__": effect_classes.__name__})
        monkeypatch.setattr(effect_classes, "TranslocationToUnannotatedContig", cls, raising=False)
    assert fresh_registry()[(cls.__module__, cls.__name__)] is cls


def test_the_registry_builds_without_the_newer_class(fresh_registry, monkeypatch):
    monkeypatch.delattr(effect_classes, "TranslocationToUnannotatedContig", raising=False)
    registry = fresh_registry()
    intergenic = effect_classes.TranslocationToIntergenic
    assert registry[(intergenic.__module__, intergenic.__name__)] is intergenic
    assert not any(name == "TranslocationToUnannotatedContig" for _, name in registry)
