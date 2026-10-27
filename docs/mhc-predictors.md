# MHC binding predictors

Vaxrank runs MHC binding predictors through
[MHCtools](https://github.com/openvax/mhctools). Select one with
`--mhc-predictor <name>`. Names marked "requires local install" need the
vendor's own executable on your `PATH`. The historical IEDB names also run
local tools; they do not call a web API or fall back to HTTP.

MHCflurry also needs model weights (`mhcflurry-downloads fetch`). Check the
[MHCtools predictor guide](https://openvax.github.io/mhctools/) for model
installation and allele/length support.

| `--mhc-predictor` | Tool | MHC Class | Notes |
|--------------------|------|-----------|-------|
| `mhcflurry` | [MHCflurry](https://github.com/openvax/mhcflurry) | I | Open-source neural network; Python package installed with MHCtools; fetch model weights separately |
| `bigmhc` | [BigMHC](https://github.com/KarchinLab/bigmhc) | I | Auto-detects EL or IM model |
| `bigmhc-el` | [BigMHC](https://github.com/KarchinLab/bigmhc) EL | I | Presentation (eluted ligand) model |
| `bigmhc-im` | [BigMHC](https://github.com/KarchinLab/bigmhc) IM | I | Immunogenicity model |
| `pepsickle` | [Pepsickle](https://github.com/pdxgx/pepsickle) | -- | Proteasomal cleavage predictor |
| `netmhc` | [NetMHC](https://services.healthtech.dtu.dk/services/NetMHC-4.0/) | I | Auto-detects NetMHC3 or NetMHC4 |
| `netmhc3` | NetMHC 3.x | I | Requires local install |
| `netmhc4` | [NetMHC 4.0](https://services.healthtech.dtu.dk/services/NetMHC-4.0/) | I | Requires local install |
| `netmhcpan` | [NetMHCpan](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) | I | Auto-detects installed version |
| `netmhcpan28` | NetMHCpan 2.8 | I | Requires local install |
| `netmhcpan3` | NetMHCpan 3.x | I | Requires local install |
| `netmhcpan4` | [NetMHCpan 4.0](https://services.healthtech.dtu.dk/services/NetMHCpan-4.0/) | I | Default mode (EL + BA) |
| `netmhcpan4-ba` | NetMHCpan 4.0 | I | Binding affinity mode only |
| `netmhcpan4-el` | NetMHCpan 4.0 | I | Eluted ligand mode only |
| `netmhcpan41` | [NetMHCpan 4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) | I | Default mode (EL + BA) |
| `netmhcpan41-ba` | NetMHCpan 4.1 | I | Binding affinity mode only |
| `netmhcpan41-el` | NetMHCpan 4.1 | I | Eluted ligand mode only |
| `netmhcpan42` | NetMHCpan 4.2 | I | Default mode (EL + BA) |
| `netmhcpan42-ba` | NetMHCpan 4.2 | I | Binding affinity mode only |
| `netmhcpan42-el` | NetMHCpan 4.2 | I | Eluted ligand mode only |
| `netmhccons` | [NetMHCcons](https://services.healthtech.dtu.dk/services/NetMHCcons-1.1/) | I | Requires local install |
| `netmhcstabpan` | [NetMHCstabpan](https://services.healthtech.dtu.dk/services/NetMHCstabpan-1.0/) | I | Stability predictor; requires local install |
| `netchop` | [NetChop](https://services.healthtech.dtu.dk/services/NetChop-3.1/) | -- | Proteasomal cleavage predictor |
| `netmhciipan` | [NetMHCIIpan](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/) | II | Auto-detects installed version |
| `netmhciipan3` | NetMHCIIpan 3.x | II | Requires local install |
| `netmhciipan4` | [NetMHCIIpan 4.0](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.0/) | II | Default mode (EL + BA) |
| `netmhciipan4-ba` | NetMHCIIpan 4.0 | II | Binding affinity mode only |
| `netmhciipan4-el` | NetMHCIIpan 4.0 | II | Eluted ligand mode only |
| `netmhciipan43` | [NetMHCIIpan 4.3](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/) | II | Default mode (EL + BA) |
| `netmhciipan43-ba` | NetMHCIIpan 4.3 | II | Binding affinity mode only |
| `netmhciipan43-el` | NetMHCIIpan 4.3 | II | Eluted ligand mode only |
| `mixmhcpred` | [MixMHCpred](https://github.com/GfellerLab/MixMHCpred) | I | Requires local install |
| `netmhcpan-iedb` | NetMHCpan compatibility wrapper | I | Requires local NetMHCpan 4.1 BA |
| `netmhccons-iedb` | NetMHCcons compatibility wrapper | I | Requires local NetMHCcons |
| `netmhciipan-iedb` | NetMHCIIpan compatibility wrapper | II | Requires local NetMHCIIpan 4.3 BA |
| `smm-iedb` | SMM compatibility wrapper | I | Requires local SMM |
| `smm-pmbec-iedb` | SMM-PMBEC compatibility wrapper | I | Requires local SMM-PMBEC |
| `random` | Random | -- | Returns random scores; for testing only |

## Optional Weirdo category scores

Topiary 5.95 supports the released Weirdo model through
`pip install 'topiary[weirdo]>=5.95,<6'`. Its `foreignness` predictions are
allele-independent and retain their exact model/artifact version and score
policy through Topiary tables, native `CandidateEpitope` records, and saved
Vaxrank datasets. Imported matched-comparator rows also retain their own
measurement context.

For a dataset that already contains Weirdo predictions, opt into the score
with `epitopes.score_expr=peptide_view(foreignness[weirdo].score)`. Ordinary
binding-based scoring is unchanged when these observations are merely present.
Keep the complete `predictor_version` and `measurement_context` columns when
transporting the table; native dataset save/load preserves them automatically.

The released policy is `human-self-complement-v1`: one minus the model's
human-category output. It accepts canonical 8–11-residue peptides and only a
human reference host (taxon 9606). Do not apply this policy as a mouse/B16
self reference; a rodent training category is not a mouse-specific model.
This is an experimental category score, with no confirmed immune-response
ranking benefit. It is not a self-safety verdict or a substitute for binding,
presentation, or self-proteome evidence. It is also distinct from the published
IEDB similarity feature used by `tesla-recognition-v1`.
