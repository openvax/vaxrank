# Shared fixture mechanics

Acquisition, original-record selection, transport, and integrity mechanics
come from Osteosarc (migrated at 0.2.3). Since 3.26, Vaxrank's Sid reads come
from openvax-v1, the OpenVax libraries' shared Sid test data published by
osteosarc 0.11 ([iskandr/osteosarc#56](https://github.com/iskandr/osteosarc/issues/56)),
and Vaxrank no longer acquires or packages reads itself. Keep the reviewed
recipe format, source/correction/reference pins, data bytes and scientific
expectations unchanged.

`vaxrank.sid_test_data` selects each cohort's historical records by exact
digest and multiplicity from its openvax-v1 member, then writes them with the
reviewed header and order, so the regressions run on unchanged inputs. This
consumer uses published Osteosarc 0.12.x and Topiary releases. See
[TEST_DATA.md](TEST_DATA.md) and osteosarc's
[shared test data](https://iskandr.github.io/osteosarc/test-data/#shared-test-data-openvax-v1).
