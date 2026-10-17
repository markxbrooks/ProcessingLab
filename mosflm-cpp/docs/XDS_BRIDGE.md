# XDS as an optional external bridge

ipmosflm and mosflm-cpp **do not** require [XDS](https://wiki.uni-konstanz.de/xds/index.php/Installation) to build or run. XDS is useful when you need to **convert** an unsupported detector format into miniCBF (or process HDF5 with vendor plugins), then feed the result into MOSFLM’s existing CBFlib path ([CBF_STATUS.md](CBF_STATUS.md)).

## When to use it

| Situation | Suggested workflow |
| --------- | ------------------ |
| Native CBF / miniCBF | Read directly (`DETECTOR CBF` / Pilatus / `cbfwrap_`) |
| Format XDS understands but MOSFLM does not | Convert with XDS `2cbf` (or related tools) → miniCBF → MOSFLM |
| Eiger HDF5 | Prefer MOSFLM `fromhdf5`, **or** convert frames externally (`H5ToXds` / `hdf2mini-cbf`) if needed |
| `*.cbf.gz` | XDS `xds-zcbf` plugin can speed XDS; for MOSFLM, decompress or convert first |

## Install notes (operators only)

Follow the [XDSwiki Installation](https://wiki.uni-konstanz.de/xds/index.php/Installation) page:

- Academic download via Dectris form; commercial users need a subscription.
- macOS Apple Silicon: `XDS-Apple_M1.zip` → `/usr/local/bin` or `$HOME/bin`.
- Quarantine: `sudo xattr -dr com.apple.quarantine XDS-Apple_M1/*` if Gatekeeper blocks binaries.
- HDF5 plugins (Durin / Neggia) and helpers (`generate_XDS.INP`, XDSGUI) are optional extras from the wiki `get_folder_silicon.sh` scripts.
- CCP4 helps XDSSTAT / XDSCONV / POINTLESS from XDSGUI; it is **not** required to compile ipmosflm’s CBF reader.

## Explicit non-goals

- Do **not** vendor XDS binaries into this repository.
- Do **not** add Neggia/Durin as link-time dependencies of `bin/ipmosflm` or `mosflm-cpp`.
- Do **not** treat `2cbf` as a runtime requirement for every image open.
