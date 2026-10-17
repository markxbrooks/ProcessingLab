# CBF fixtures

| File | Notes |
| ---- | ----- |
| `sample_minicbf.cbf` | Symlink to a local miniCBF (`KIT-02AT05_w1_1_0001.cbf` from mxflask test data when present) |

## Smoke test

From `mosflm-cpp/tools`:

```bash
cc -O0 -g -I../../cbf/include -o cbf_smoke cbf_smoke.c ../../cbf/lib/libcbf.a -lm
./cbf_smoke ../testdata/cbf/sample_minicbf.cbf
```

Expect a line starting with `CBF_SMOKE_OK`.
